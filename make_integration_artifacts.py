"""Generate tables and plot data from the local model-calculation records.

Table endpoints are rounded outward from exact rational strings. Optional PNG
plots show recorded measurements; no integration or parameter fitting is run.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
import os
from pathlib import Path
import sys
from typing import Any

from derive_query_contrasts import derive, decimal_bound

HERE = Path(__file__).resolve().parent
CASES = ['balanced_8', 'balanced_16', 'left_branching_16', 'stiff_balanced_16']
SHORT = dict(zip(CASES, ['B8', 'B16', 'L16', 'S16']))
QUERIES = ['HH0H', 'HHLH', 'HLLH', *CASES]


def build_tables(root: Path) -> tuple[dict[str, str], dict[str, Any], dict[str, Any]]:
    """Check recorded values and construct all table bodies and refinement data."""
    checks, contrast_table = derive(root)
    parameters = json.loads((root / 'published_timing_parameters.json').read_text())
    comparison = json.loads((root / 'results/moment_comparison.json').read_text())
    weights = parameters['weights_in_column_order_AU_to_BB']
    target_keys = ['mora', 'syllable', 'grid', 'hemimetron']
    target_labels = ['Mora', 'Syllable', 'Grid column', 'Hemimetron']
    rows = [r'\begin{tabular}{lrrr}\toprule',
            r'Target type & Target (ms) & Above & Below\\',
            r' & & \multicolumn{2}{c}{Weights in $10^{-4}\,\mathrm{ms}^{-2}$}\\\midrule']
    for i, (key, label) in enumerate(zip(target_keys, target_labels)):
        values = [parameters['targets'][key], *weights[2*i:2*i+2]]
        rows.append(label + ' & ' + ' & '.join(format(float(F(x)), '.6f') for x in values) + r' \\')
    rows.append(r'\bottomrule\end{tabular}')
    tables = {'parameter_table.tex': '\n'.join(rows) + '\n'}
    rows = [r'\begin{tabular}{llrrrr}\toprule',
            r'Case & Arrays & $h$ & Log width & Time & Peak\\\midrule']
    refinement: dict[str, Any] = {
        'description': 'Completed two- and three-moment normalizer meshes; exact rational coordinates',
        'target': '1/1000', 'cases': {}}
    for case in CASES:
        refinement['cases'][case] = {}
        for method, count in [('two', 2), ('three', 3)]:
            record = comparison['cases'][case][method]
            last = record['attempts'][-1]
            rows.append(f"{SHORT[case]} & {count} & {float(F(last['mesh_ms'])):g} & "
                        f"{decimal_bound(F(last['log_width_upper_exact_rational']), 7, upper=True)} & "
                        f"{record['process_seconds']:.2f} & "
                        f"{last['peak_resident_bytes_after_completed_mesh']/1e6:.1f} " + r'\\')
            refinement['cases'][case][method] = {
                'points': [{k: attempt[k] for k in ('mesh_ms', 'log_width_upper_exact_rational')}
                           for attempt in record['attempts']]}
    rows.append(r'\bottomrule\end{tabular}')
    tables['normalizer_table.tex'] = '\n'.join(rows) + '\n'
    rows = [r'\begin{tabular}{lrcrr}\toprule',
            r'Case & $h$ (ms) & Enclosure of $\Delta$ & Moments (s) & FFT (s)\\\midrule']
    for case in QUERIES:
        record = json.loads((root / f'results/queries/{case}.json').read_text())
        fft = json.loads((root / f'results/point_queries/{case}.json').read_text())
        last = record['attempts'][-1]
        interval = last['probabilities']['Delta']
        rows.append(f"{SHORT.get(case, case)} & {float(F(last['mesh_ms'])):g} & "
                    f"$[{decimal_bound(F(interval['lower']), 7)}, {decimal_bound(F(interval['upper']), 7, upper=True)}]$ & "
                    f"{record['process_elapsed_seconds']:.3f} & {fft['process_elapsed_seconds']:.3f} " + r'\\')
    rows.append(r'\bottomrule\end{tabular}')
    tables['query_table.tex'] = '\n'.join(rows) + '\n'
    record = json.loads((root / 'results/queries/HH0H.json').read_text())
    probabilities = record['attempts'][-1]['probabilities']
    rows = [r'\begin{tabular}{lc}\toprule Quantity & Enclosure\\\midrule']
    for label, name in zip([r'$P(A)$', r'$P(B)$', r'$P(A\cap B)$', r'$P(A\mid B)$', r'$\Delta$'],
                           ['P_A', 'P_B', 'P_AB', 'P_A_given_B', 'Delta']):
        value = probabilities[name]
        rows.append(label + ' & $[' + decimal_bound(F(value['lower']), 7)
                    + ', ' + decimal_bound(F(value['upper']), 7, upper=True) + r']$ \\')
    rows.append(r'\bottomrule\end{tabular}')
    tables['worked_query.tex'] = '\n'.join(rows) + '\n'
    tables['query_contrasts.tex'] = contrast_table
    return tables, refinement, checks


def plot_refinement(data: dict[str, Any], output: Path) -> dict[str, Any]:
    """Create one plot per model with Matplotlib's default font and colour cycle."""
    os.environ.setdefault('MPLCONFIGDIR', str(output / 'matplotlib_cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for case in CASES:
        fig, ax = plt.subplots(figsize=(6, 4))
        for method, label, marker in [('two', 'Two moments', 's'), ('three', 'Three moments', 'o')]:
            points = data['cases'][case][method]['points']
            ax.plot([float(F(p['mesh_ms'])) for p in points],
                    [float(F(p['log_width_upper_exact_rational'])) for p in points],
                    marker=marker, label=label)
        ax.axhline(float(F(data['target'])), linestyle=':', label='Target')
        ax.set_xscale('log', base=2)
        ax.set_yscale('log')
        ax.set_xlabel('Mesh h (ms)')
        ax.set_ylabel('Log-normalizer enclosure width')
        ax.set_title(SHORT[case])
        ax.legend()
        fig.tight_layout()
        fig.savefig(output / f'refinement_{case}.png', dpi=180)
        plt.close(fig)
    return {'matplotlib': matplotlib.__version__, 'plot_count': len(CASES),
            'scope': 'Completed recorded meshes; unfinished timed-out refinements are not plotted'}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=HERE)
    parser.add_argument('--output', type=Path, required=True, help='A new directory')
    parser.add_argument('--plots', action='store_true', help='Also create PNG refinement plots; requires Matplotlib')
    args = parser.parse_args()
    output = args.output.resolve()
    try:
        output.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        parser.error(f'Cannot create output directory: {exc}')
    try:
        tables, data, checks = build_tables(args.root.resolve(strict=True))
        for name, text in tables.items():
            (output / name).write_text(text)
        (output / 'refinement_data.json').write_text(json.dumps(data, indent=2) + '\n')
        report = {'status': 'PASSED', 'table_count': len(tables),
                  'completed_mesh_points': sum(len(m['points']) for c in data['cases'].values() for m in c.values()),
                  'record_checks': checks,
                  'scope': 'Exact-rational record analysis and descriptive output; no numerical integration or retiming'}
        if args.plots:
            report['plots'] = plot_refinement(data, output)
        (output / 'artifact_checks.json').write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps({k: report[k] for k in ['status', 'table_count', 'completed_mesh_points']}))
        return 0
    except Exception as exc:
        (output / 'artifact_checks.json').write_text(json.dumps({
            'status': 'FAILED', 'error': f'{type(exc).__name__}: {exc}'}, indent=2) + '\n')
        print(f'FAILED: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
