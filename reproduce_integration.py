"""Run integrity checks, result analysis, controls and model calculations."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path, PurePosixPath
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
RUNTIME_DIRECTORIES = {'.git', '.venv', '.lake', '.deps', '__pycache__', 'runs'}


def file_identity(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            digest.update(chunk)
    return {'name': path.name, 'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def external_gubpi(args, output):
    """Execute one specified case, preserving failures and bounded teardown."""
    from fractions import Fraction
    import mpmath as mp
    from flint import ctx
    from integration_bounds import hayes_models, tail_envelope, ball
    from integration_scaling import build_case
    import gubpi_encoding as enc
    if not args.dotnet or not args.gubpi_dll:
        raise ValueError('--dotnet and --gubpi-dll are required; dependencies are not downloaded')
    runtime, dll = args.dotnet.resolve(strict=True), args.gubpi_dll.resolve(strict=True)
    ctx.prec, mp.mp.dps = 128, 90
    model = hayes_models()[args.case] if args.case != 'balanced_8' else build_case(args.case)
    tail = tail_envelope(model, Fraction(800))
    proposal_constant = None
    if args.encoding == 'gaussian':
        proposal_constant = ball(1)
        for _, sd in enc.proposal_parameters(model):
            proposal_constant *= (2 * ball(1).pi()).sqrt() * ball(sd)
    settings = ('0.2', '0.1', '0.05') if args.encoding in ('pdf', 'arithmetic-polytopes') else (100, 1000, 10000)
    env = os.environ.copy()
    env.update(DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1',
               DOTNET_SYSTEM_GLOBALIZATION_INVARIANT='1', LANG='C', LC_ALL='C',
               DOTNET_ROOT=str(runtime.parent))
    for key in ('TMPDIR', 'XDG_DATA_HOME', 'XDG_CONFIG_HOME', 'XDG_CACHE_HOME',
                'DOTNET_CLI_HOME', 'NUGET_PACKAGES', 'NUGET_HTTP_CACHE_PATH'):
        folder = output / 'runtime_cache' / key.lower()
        folder.mkdir(parents=True)
        env[key] = str(folder)
    record = {'case': args.case, 'encoding': args.encoding, 'budget_seconds': 60,
              'attempts': [], 'scope': 'Specified third-party floating computation, not certified model equivalence',
              'dependency_files': {'runtime': file_identity(runtime), 'gubpi_dll': file_identity(dll)},
              'dependency_identity_limit': 'File fingerprints, not complete transitive build identity; runtime version and source commit are not independently probed'}
    elapsed = 0.0
    for setting in settings:
        dest = output / str(setting)
        dest.mkdir()
        if args.encoding == 'pdf':
            program, factors = enc.gaussian_program(model, setting)
            checks = enc.check_gaussian_encoding(factors)
        else:
            mode = 'polytopes' if args.encoding == 'arithmetic-polytopes' else 'boxes'
            proposal = 'gaussian' if args.encoding == 'gaussian' else 'uniform'
            program, factors = enc.arithmetic_program(model, mode, setting, proposal)
            checks = {'exact_control_points': enc.check_arithmetic_encoding(model, factors)}
        model_file = dest / 'model.spcf'
        model_file.write_text(program)
        save(dest / 'encoding.json', {'factors': factors, 'checks': checks})
        remaining = 60 - elapsed
        if remaining <= 0:
            break
        started, timed_out, exit_code = time.perf_counter(), False, None
        with (dest / 'stdout.log').open('w') as log:
            process = subprocess.Popen([str(runtime), str(dll), str(model_file)], cwd=dest,
                                       env=env, stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True)
            try:
                exit_code = process.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                timed_out = True
                os.killpg(process.pid, signal.SIGKILL)
                try:
                    exit_code = process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    pass
        seconds = time.perf_counter() - started
        elapsed += seconds
        row = {'setting': setting, 'seconds': seconds, 'cumulative_process_seconds': elapsed,
               'exit_code': exit_code, 'timeout': timed_out,
               'termination_confirmed': process.poll() is not None}
        within_budget = elapsed <= record['budget_seconds']
        row['within_time_budget'] = within_budget
        if not timed_out and exit_code == 0:
            row.update(enc.parse_output((dest / 'stdout.log').read_text(), len(model.leaves), tail, proposal_constant))
        if timed_out or not within_budget:
            row['status'] = 'TIME_BUDGET_EXCEEDED'
        elif exit_code != 0:
            row['status'] = 'PROCESS_ERROR'
        record['attempts'].append(row)
        save(output / 'gubpi_results.json', record)
        if timed_out or row['status'] not in ('REPORTED_BOUNDS_PARSED', 'INVALID_OR_NONFINITE_NORMALIZER') or row.get('reported_target_met'):
            break
    record['target_met'] = any(row.get('reported_target_met', False)
                               and row['within_time_budget']
                               and row['termination_confirmed']
                               and row['exit_code'] == 0
                               and not row['timeout']
                               for row in record['attempts'])
    record['completed_attempts_handled'] = True
    save(output / 'gubpi_results.json', record)
    return record


def save(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')


def verify_inventory(*, clean=False):
    """Check a complete payload inventory, not just whichever hashes are listed.

    Runtime directories are excluded only for a working checkout. A clean
    distribution audit rejects them too. Checksums detect corruption; they do
    not authenticate a project against a party able to rewrite all metadata.
    """
    manifest_path, sums_path = HERE / 'manifest.json', HERE / 'SHA256SUMS'
    if manifest_path.is_symlink() or sums_path.is_symlink():
        raise ValueError('Symlinked project metadata')
    manifest_raw = manifest_path.read_bytes()
    manifest = json.loads(manifest_raw)
    expected, sums = set(), []
    for item in manifest['files']:
        name = item['path']
        rel = PurePosixPath(name)
        if (not isinstance(name, str) or not name or rel.is_absolute()
                or '..' in rel.parts or rel.as_posix() != name or '\\' in name
                or '\n' in name or '\r' in name or name in expected
                or name in {'manifest.json', 'SHA256SUMS'}
                or rel.parts[0] in RUNTIME_DIRECTORIES):
            raise ValueError('Unsafe manifest member')
        path = HERE / name
        if (any(parent.is_symlink() for parent in (path, *path.parents) if parent != HERE)
                or not path.resolve().is_relative_to(HERE) or not path.is_file()):
            raise ValueError('Unsafe or missing project member: ' + name)
        raw = path.read_bytes()
        if (hashlib.sha256(raw).hexdigest() != item['sha256']
                or type(item['bytes']) is not int or len(raw) != item['bytes']):
            raise ValueError('Release hash/size mismatch: ' + name)
        expected.add(name)
        sums.append(f"{item['sha256']}  {name}\n")
    sums.append(hashlib.sha256(manifest_raw).hexdigest() + '  manifest.json\n')
    if sums_path.read_text() != ''.join(sums):
        raise ValueError('SHA256SUMS and manifest disagree')
    actual = set()
    def traversal_error(error):
        raise error
    for current, directories, files in os.walk(HERE, followlinks=False, onerror=traversal_error):
        root = Path(current)
        if root == HERE and not clean:
            directories[:] = [d for d in directories if d not in RUNTIME_DIRECTORIES]
        for name in directories + files:
            path = root / name
            if path.is_symlink():
                raise ValueError('Symlink in project payload: ' + str(path.relative_to(HERE)))
        actual.update((root / name).relative_to(HERE).as_posix() for name in files)
    if actual != expected | {'manifest.json', 'SHA256SUMS'}:
        raise ValueError('Release inventory has missing or unlisted files')
    if clean and any((HERE / name).exists() for name in RUNTIME_DIRECTORIES):
        raise ValueError('Runtime directory in clean project')
    return len(expected)


def run(args):
    if sys.flags.optimize:
        raise RuntimeError('Do not disable assertion-based numerical controls')
    output = args.output.resolve()
    if output.is_relative_to(HERE) and not output.is_relative_to(HERE / 'runs'):
        raise ValueError('Put local run outputs under runs/, or outside the project tree')
    if getattr(args, 'clean_tree', False) and output.is_relative_to(HERE):
        raise ValueError('A clean audit requires its output outside the project tree')
    output.mkdir(parents=True, exist_ok=False)
    report = {'started_at': datetime.now(timezone.utc).isoformat(), 'mode': args.mode,
              'status': 'RUNNING', 'limits': 'File checks, local controls and model calculations are distinct validation scopes'}
    try:
        report['files_checked'] = verify_inventory(clean=getattr(args, 'clean_tree', False))
        if args.mode in ('check', 'proofs', 'all', 'experiments'):
            mode = {'check': 'numerical', 'proofs': 'proofs', 'all': 'all', 'experiments': 'experiments'}[args.mode]
            command = [sys.executable, str(HERE / 'verify_integration.py'), '--mode', mode,
                       '--output', str(output / 'verification')]
            if args.mathlib:
                command += ['--mathlib', str(args.mathlib.resolve())]
            if mode in ('all', 'proofs') and not args.mathlib:
                raise ValueError('--mathlib is required for formal checks')
            with (output / 'verification.log').open('w') as log:
                cp = subprocess.run(command, cwd=output, stdout=log, stderr=subprocess.STDOUT)
            report['verification_exit_code'] = cp.returncode
            if cp.returncode:
                raise RuntimeError('Verification failed; inspect verification.log')
        if args.mode in ('artifacts', 'figures', 'all'):
            command = [sys.executable, str(HERE / 'make_integration_artifacts.py'),
                       '--output', str(output / 'artifacts')]
            if args.mode in ('figures', 'all'):
                command.append('--plots')
            with (output / 'artifacts.log').open('w') as log:
                cp = subprocess.run(command, cwd=output, stdout=log, stderr=subprocess.STDOUT)
            if cp.returncode:
                raise RuntimeError('Result export failed; inspect artifacts.log')
            report['artifacts'] = json.loads((output / 'artifacts/artifact_checks.json').read_text())
        if args.mode == 'analyze':
            from derive_query_contrasts import derive
            analysis, table = derive(HERE)
            save(output / 'derivation_checks.json', analysis)
            (output / 'query_contrasts.tex').write_text(table)
            report['analysis'] = analysis
        if args.mode == 'range':
            report['runs'] = []
            for method in ('scalar_factors', 'grid_recurrence'):
                for case in ('balanced_8', 'balanced_16', 'left_branching_16', 'stiff_balanced_16'):
                    dest = output / (method + '__' + case + '.json')
                    command = [sys.executable, str(HERE / 'integration_scaling.py'), '--worker', case,
                               '--method', method, '--output', str(dest)]
                    started, timed_out = time.perf_counter(), False
                    with dest.with_suffix('.log').open('w') as log:
                        try:
                            cp = subprocess.run(command, cwd=output, stdout=log, stderr=subprocess.STDOUT, timeout=60)
                            code = cp.returncode
                        except subprocess.TimeoutExpired:
                            code, timed_out = None, True
                    elapsed = time.perf_counter() - started
                    data = json.loads(dest.read_text()) if dest.exists() else {}
                    status = ('TIME_BUDGET_EXCEEDED' if timed_out else
                              'FEASIBILITY_TARGET_MET' if code == 0 and elapsed <= 60 and
                              data.get('status') == 'TOLERANCE_REACHED_PENDING_PROCESS_BUDGET_CHECK'
                              else 'FAILED_OR_UNRESOLVED')
                    report['runs'].append({'method': method, 'case': case, 'status': status,
                                           'seconds': elapsed, 'exit_code': code, 'result': dest.name})
                    save(output / 'run_checks.json', report)
            if any(x['status'] != 'FEASIBILITY_TARGET_MET' for x in report['runs']):
                raise RuntimeError('Some range-only cases did not meet the target; all attempt records are retained')
        if args.mode == 'gubpi-inputs':
            from fractions import Fraction
            import mpmath as mp
            from flint import ctx
            from integration_bounds import hayes_models
            from integration_scaling import build_case
            import gubpi_encoding as encoding
            ctx.prec, mp.mp.dps = 128, 90
            report['inputs'] = []
            for case in ('HH0H', 'HHLH', 'HLLH', 'balanced_8'):
                model = hayes_models()[case] if case != 'balanced_8' else build_case(case)
                for eps in ('0.2', '0.1', '0.05'):
                    program, factors = encoding.gaussian_program(model, eps)
                    checks = encoding.check_gaussian_encoding(factors)
                    dest = output / (case + '__pdf_' + eps + '.spcf')
                    dest.write_text(program)
                    report['inputs'].append({'file': dest.name, 'case': case, 'encoding': 'pdf', 'checks': checks})
                for proposal in ('uniform', 'gaussian'):
                    for splits in (100, 1000, 10000):
                        program, factors = encoding.arithmetic_program(model, 'boxes', splits, proposal)
                        checks = encoding.check_arithmetic_encoding(model, factors)
                        dest = output / (case + '__' + proposal + '_' + str(splits) + '.spcf')
                        dest.write_text(program)
                        report['inputs'].append({'file': dest.name, 'case': case, 'encoding': proposal, 'exact_control_points': checks})
            report['gubpi_limit'] = 'Input-generation checks only; see GUBPI.md. Gaussian factor tests are finite, not a uniform equivalence proof.'
        if args.mode == 'gubpi':
            report['gubpi'] = external_gubpi(args, output)
            if not report['gubpi']['target_met']:
                raise RuntimeError('The specified GuBPI case did not reach the .001 target; all attempts retained')
        report['status'] = 'PASSED_REQUESTED_SCOPE'
        code = 0
    except Exception as exc:
        report.update(status='FAILED_OR_UNRESOLVED', error=type(exc).__name__ + ': ' + str(exc))
        code = 1
    report['completed_at'] = datetime.now(timezone.utc).isoformat()
    save(output / 'run_checks.json', report)
    print(json.dumps({k: report.get(k) for k in ('mode', 'status', 'error')}), flush=True)
    return code


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['audit', 'analyze', 'artifacts', 'check', 'proofs', 'all', 'figures', 'experiments', 'range', 'gubpi-inputs', 'gubpi'], default='check')
    parser.add_argument('--clean-tree', action='store_true', help='Reject runtime/cache directories as well as unlisted payload files')
    parser.add_argument('--mathlib', type=Path)
    parser.add_argument('--dotnet', type=Path)
    parser.add_argument('--gubpi-dll', type=Path)
    parser.add_argument('--case', choices=['HH0H', 'HHLH', 'HLLH', 'balanced_8'], default='HH0H')
    parser.add_argument('--encoding', choices=['pdf', 'uniform', 'gaussian', 'arithmetic-polytopes'], default='pdf')
    parser.add_argument('--output', type=Path, required=True)
    raise SystemExit(run(parser.parse_args()))
