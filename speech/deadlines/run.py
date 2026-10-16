"""Compute and verify numerical decisions of a speech-duration planner."""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import argparse
import json
import math
import platform
import sys
import time

import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
from flint import ctx

HERE = Path(__file__).resolve().parent
APP = HERE.parent
sys.path.insert(0, str(APP.parent))
from integration_bounds import Node, Potential, fraction_string, tail_envelope
from integration_queries import event_integral
from integration_scaling import describe, save
from integration_point_queries import point_log_integral

ALPHA = F(19, 20)
FRAME = 10
CUBE = F(1600)
MESHES = [20., 10., 5., 2.5, 1.25]


def inputs():
    return json.loads((APP/'test_requests.json').read_text())


def tree(row):
    leaves = tuple(Node(Potential(((F(m), F(w), F(w)),)), name=f'x{i+1}')
                   for i, (m, w) in enumerate(zip(row['target_ms_exact'], row['unary_weight_exact'])))
    w = F(row['group_weight_exact'])
    return Node(Potential(((F(row['total_ms']), w, w),)), leaves, name='word')


def last_frame(row):
    return int(F(row['total_ms'])//FRAME)


def binary_search(cdf, last):
    if cdf(last*FRAME) < float(ALPHA):
        return None
    lo, hi = 0, last
    while hi-lo > 1:
        mid = (lo+hi)//2
        if cdf(mid*FRAME) >= float(ALPHA):
            hi = mid
        else:
            lo = mid
    return hi*FRAME


def quadratic_cdf(row):
    mu = np.array([float(F(x)) for x in row['target_ms_exact']])
    unary = np.array([2*float(F(w)) for w in row['unary_weight_exact']])
    shared = 2*float(F(row['group_weight_exact']))
    covariance = np.linalg.inv(np.diag(unary)+shared*np.ones((2, 2)))
    center = covariance@(mu*unary+shared*float(F(row['total_ms'])))
    scale = math.sqrt(covariance[0, 0])
    conditional_scale = math.sqrt(covariance[1, 1]-covariance[0, 1]**2/covariance[0, 0])
    def integrand(z):
        return math.exp(-.5*z*z)/math.sqrt(2*math.pi)*ndtr(
            (center[1]+covariance[0, 1]/scale*z)/conditional_scale)
    z0 = -center[0]/scale
    denominator, den_error = quad(integrand, z0, np.inf, epsabs=1e-12, epsrel=1e-10, limit=150)
    if not denominator > 0:
        raise ArithmeticError('Nonpositive quadratic normalizer')
    estimates = {}
    def cdf(deadline):
        if deadline <= 0:
            return 0.
        if deadline not in estimates:
            numerator, error = quad(integrand, z0, (deadline-center[0])/scale,
                                    epsabs=1e-12, epsrel=1e-10, limit=150)
            estimates[deadline] = {'probability': max(0., min(1., numerator/denominator)),
                                   'estimated_ratio_error': error/denominator+numerator*den_error/denominator**2}
        return estimates[deadline]['probability']
    return cdf, estimates


def ordinary_word(row):
    node = tree(row)
    last = last_frame(row)
    outputs = {}
    functions = {}
    for mesh in MESHES:
        start = time.perf_counter()
        diagnostics = {'negative_fft_coefficients': 0, 'clipped_negative_mass_sum': 0.}
        denominator = point_log_integral(node, mesh, float(CUBE), {}, diagnostics)
        estimates = {}
        def cdf(deadline, mesh=mesh, denominator=denominator, estimates=estimates, diagnostics=diagnostics):
            if deadline <= 0:
                return 0.
            if deadline not in estimates:
                numerator = point_log_integral(node, mesh, float(CUBE), {'x1': (F(0), F(deadline))}, diagnostics)
                estimates[deadline] = math.exp(numerator-denominator)
            return estimates[deadline]
        deadline = binary_search(cdf, last)
        checked = [last*FRAME] if deadline is None else [deadline-FRAME, deadline]
        endpoints = [cdf(x) for x in checked]
        outputs[str(mesh)] = {'deadline_ms': deadline, 'decision_cdf': dict(zip(checked, endpoints)),
                              'seconds': time.perf_counter()-start, 'all_cdf_queries': estimates,
                              'diagnostics': diagnostics}
        functions[mesh] = cdf
    stability = {'status': 'NO_STABLE_PAIR', 'deadline_ms': None}
    for old, new in zip(MESHES, MESHES[1:]):
        a, b = outputs[str(old)], outputs[str(new)]
        if a['deadline_ms'] != b['deadline_ms']:
            continue
        times = [last*FRAME] if b['deadline_ms'] is None else [b['deadline_ms']-FRAME, b['deadline_ms']]
        differences = {t: abs(functions[old](t)-functions[new](t)) for t in times}
        if max(differences.values()) <= .001:
            stability = {'status': 'STABLE_PAIR', 'deadline_ms': b['deadline_ms'],
                         'meshes_ms': [old, new], 'cdf_differences': differences}
            break
    start = time.perf_counter()
    cdf, estimates = quadratic_cdf(row)
    deadline = binary_search(cdf, last)
    checked = [last*FRAME] if deadline is None else [deadline-FRAME, deadline]
    for d in checked:
        cdf(d)
    quadratic = {'deadline_ms': deadline, 'all_cdf_queries': estimates, 'seconds': time.perf_counter()-start}
    return {'id': row['id'], 'word': row['word'], 'speaker': row['speaker'],
            'last_admissible_deadline_ms': last*FRAME, 'midpoint_fft': outputs,
            'stability': stability, 'quadratic': quadratic}


def ordinary():
    start = time.perf_counter()
    rows = inputs()
    results = [ordinary_word(r) for r in rows]
    save(HERE/'ordinary.json', {'requests': len(rows), 'wall_seconds': time.perf_counter()-start, 'rows': results})
    # This is descriptive only; no result is yet called an error.
    summary = {}
    for h in MESHES:
        summary[str(h)] = sum(r['midpoint_fft'][str(h)]['deadline_ms'] != r['quadratic']['deadline_ms'] for r in results)
    summary['stability'] = sum(r['stability']['status'] == 'STABLE_PAIR' and r['stability']['deadline_ms'] != r['quadratic']['deadline_ms'] for r in results)
    print(json.dumps({'ordinary_disagreements_with_quadratic_control': summary, 'seconds': time.perf_counter()-start}), flush=True)


def encode_interval(lo, hi):
    return {'lower': str(lo), 'upper': str(hi), 'lower_approx': float(lo), 'upper_approx': float(hi)}


def certify_word(item):
    row, proposal = item
    ctx.prec = 128
    node = tree(row)
    path = HERE/'certificates'/(row['id'].replace('/', '__')+'.json')
    if path.exists():
        cached = json.loads(path.read_text())
        if (cached.get('tree') == describe(node) and cached.get('status', '').startswith('CERTIFIED_')
                and cached.get('probability_requirement') == str(ALPHA)
                and cached.get('frame_ms') == FRAME and cached.get('cube_upper_ms') == str(CUBE)
                and cached.get('precision_bits') == 128):
            return {k: cached[k] for k in ['id', 'status', 'deadline_ms', 'seconds']}
    started = time.perf_counter()
    last = last_frame(row)
    frame = last if proposal is None else proposal//FRAME
    tail = tail_envelope(node, CUBE)
    record = {'id': row['id'], 'status': 'RUNNING',
              'tree': describe(node), 'probability_requirement': str(ALPHA),
              'frame_ms': FRAME, 'last_admissible_deadline_ms': last*FRAME,
              'initial_ordinary_proposal_ms': proposal, 'precision_bits': ctx.prec,
              'cube_upper_ms': str(CUBE), 'tail_upper_exact': fraction_string(tail),
              'attempts': [], 'deadline_ms': None}
    path.parent.mkdir(exist_ok=True)
    try:
        for j in range(11):
            mesh = F(20, 2**j)
            cache = {}
            zlo, zhi = event_integral(node, mesh, CUBE, {}, tail, cache)
            current = {'mesh_ms': str(mesh), 'normalizer': {'lower': fraction_string(zlo), 'upper': fraction_string(zhi)}, 'cdf': {}}
            record['attempts'].append(current)
            def bounds(k):
                if k == 0:
                    return F(0), F(0)
                key = str(k*FRAME)
                if key not in current['cdf']:
                    alo, ahi = event_integral(node, mesh, CUBE, {'x1': (F(0), F(k*FRAME))}, tail, cache)
                    if zlo > 0:
                        lo = max(F(0), F(fraction_string((alo/zhi).lower())))
                        hi = min(F(1), F(fraction_string((ahi/zlo).upper())))
                    else:
                        lo, hi = F(0), F(1)
                    assert lo <= hi
                    current['cdf'][key] = {**encode_interval(lo, hi), 'numerator': {'lower': fraction_string(alo), 'upper': fraction_string(ahi)}}
                b = current['cdf'][key]
                return F(b['lower']), F(b['upper'])
            for _ in range(last+1):
                lower, upper = bounds(frame)
                if upper < ALPHA:
                    if frame == last:
                        record.update(status='CERTIFIED_INFEASIBLE', deadline_ms=None,
                                      cdf_at_last=encode_interval(lower, upper))
                        break
                    frame += 1
                    continue
                before_lower, before_upper = bounds(frame-1)
                if before_lower >= ALPHA:
                    frame -= 1
                    continue
                if lower >= ALPHA and before_upper < ALPHA:
                    record.update(status='CERTIFIED_DEADLINE', deadline_ms=frame*FRAME,
                                  cdf_at_deadline=encode_interval(lower, upper),
                                  cdf_at_previous=encode_interval(before_lower, before_upper))
                break
            record['seconds'] = time.perf_counter()-started
            record['final_mesh_ms'] = str(mesh)
            if record['status'] == 'RUNNING' and record['seconds'] > 180:
                record['status'] = 'TIME_ALLOWANCE_EXCEEDED_AFTER_COMPUTATION'
            save(path, record)
            if record['status'] != 'RUNNING':
                break
        if record['status'] == 'RUNNING':
            record['status'] = 'MAXIMUM_MESH_UNRESOLVED'
    except Exception as error:
        record['status'] = 'ERROR'
        record['error'] = repr(error)
    record['seconds'] = time.perf_counter()-started
    save(path, record)
    return {k: record[k] for k in ['id', 'status', 'deadline_ms', 'seconds']}


def certify(workers):
    rows = inputs()
    ordinary_rows = {r['id']: r for r in json.loads((HERE/'ordinary.json').read_text())['rows']}
    requests = [(r, ordinary_rows[r['id']]['midpoint_fft']['10.0']['deadline_ms']) for r in rows]
    started = time.perf_counter()
    run = {'started_at': datetime.now(timezone.utc).isoformat(), 'workers': workers,
           'python': platform.python_version(), 'machine': platform.machine(),
           'requests': len(rows), 'cases': []}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(certify_word, item) for item in requests]
        for future in as_completed(futures):
            record = future.result()
            run['cases'].append(record)
            save(HERE/'certification_run.json', run)
            if len(run['cases']) % 25 == 0 or not record['status'].startswith('CERTIFIED_'):
                print(json.dumps({'completed': len(run['cases']), **record}), flush=True)
    run.update(wall_seconds=time.perf_counter()-started, status='ALL_REQUESTS_HANDLED', completed_at=datetime.now(timezone.utc).isoformat())
    save(HERE/'certification_run.json', run)
    print(json.dumps({'status': run['status'], 'cases': len(run['cases']), 'seconds': run['wall_seconds']}), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ordinary', action='store_true')
    parser.add_argument('--certify', action='store_true')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.ordinary:
        ordinary()
    if args.certify:
        certify(args.workers)


if __name__ == '__main__':
    main()
