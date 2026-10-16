"""Verify exact deadline certificates and compare every ordinary decision."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import copy
import json
import statistics

from run import ALPHA, FRAME, HERE, APP, inputs, tree, describe, quadratic_cdf


def main():
    requests = inputs()
    ordinary = json.loads((HERE/'ordinary.json').read_text())
    controls = {r['id']: r for r in ordinary['rows']}
    run = json.loads((HERE/'certification_run.json').read_text())
    assert len(controls) == len(requests) == len(run['cases']) == 836
    compared = []
    methods = ['fft_20.0', 'fft_10.0', 'fft_5.0', 'fft_2.5', 'fft_1.25', 'stability', 'quadratic']
    summary = {m: {'different_decisions': 0, 'premature_deadline': 0, 'excessive_deadline': 0,
                   'false_infeasible': 0, 'false_feasible': 0, 'unresolved': 0, 'differences': []} for m in methods}
    certificates = []
    independent_queries = 0
    for row in requests:
        c = json.loads((HERE/'certificates'/(row['id'].replace('/', '__')+'.json')).read_text())
        certificates.append(c)
        assert c['id'] == row['id']
        assert c['tree'] == describe(tree(row))
        altered = copy.deepcopy(row)
        altered['syllables'][0]['duration_ms'] = '1'
        altered['syllables'][1]['duration_ms'] = str(F(row['total_ms'])-1)
        assert describe(tree(row)) == describe(tree(altered))
        assert c['status'] in ['CERTIFIED_DEADLINE', 'CERTIFIED_INFEASIBLE'], c['status']
        correct = c['deadline_ms']
        if correct is None:
            assert F(c['cdf_at_last']['upper']) < ALPHA
        else:
            assert correct % FRAME == 0 and 0 < correct <= F(row['total_ms'])
            assert F(c['cdf_at_deadline']['lower']) >= ALPHA
            assert F(c['cdf_at_previous']['upper']) < ALPHA
        cdf, _ = quadratic_cdf(row)
        for attempt in c['attempts']:
            for d, interval in attempt['cdf'].items():
                lo, hi = F(interval['lower']), F(interval['upper'])
                assert 0 <= lo <= hi <= 1
                assert lo <= F(cdf(int(d))) <= hi, (row['id'], d, interval, cdf(int(d)))
                independent_queries += 1
        r = controls[row['id']]
        item = {'id': row['id'], 'word': row['word'], 'certified_deadline_ms': correct, 'ordinary': {}}
        for method in methods:
            result = r['midpoint_fft'][method[4:]] if method.startswith('fft_') else r[method]
            predicted = result['deadline_ms']
            if method == 'stability' and result['status'] != 'STABLE_PAIR':
                summary[method]['unresolved'] += 1
                kind = 'UNRESOLVED'
            elif predicted == correct:
                kind = 'AGREES'
            else:
                kind = ('false_infeasible' if predicted is None else 'false_feasible' if correct is None
                        else 'premature_deadline' if predicted < correct else 'excessive_deadline')
                summary[method]['different_decisions'] += 1
                summary[method][kind] += 1
                summary[method]['differences'].append({'id': row['id'], 'word': row['word'],
                     'ordinary_deadline_ms': predicted, 'certified_deadline_ms': correct, 'kind': kind})
            item['ordinary'][method] = {'deadline_ms': predicted, 'comparison': kind}
        compared.append(item)
    first = summary['fft_10.0']['differences'][0] if summary['fft_10.0']['differences'] else None
    result = {'status': 'PASS', 'requests': len(requests), 'certified_feasible': sum(c['deadline_ms'] is not None for c in certificates),
              'certified_infeasible': sum(c['deadline_ms'] is None for c in certificates),
              'unresolved_certificates': 0, 'independent_cdf_checks_inside_enclosures': independent_queries,
              'methods': summary, 'final_mesh_counts': dict(Counter(c['final_mesh_ms'] for c in certificates)),
              'certificate_worker_seconds_median': statistics.median(c['seconds'] for c in certificates),
              'certificate_worker_seconds_maximum': max(c['seconds'] for c in certificates),
              'certificate_batch_wall_seconds': run['wall_seconds'],
              'ordinary_total_wall_seconds': ordinary['wall_seconds'],
              'ordinary_seconds_sum_by_mesh': {h: sum(r['midpoint_fft'][h]['seconds'] for r in controls.values()) for h in ['20.0','10.0','5.0','2.5','1.25']},
              'quadratic_seconds_sum': sum(r['quadratic']['seconds'] for r in controls.values()),
              'first_primary_disagreement': first,
              'first_primary_disagreement_ordinary': controls[first['id']] if first else None,
              'first_primary_disagreement_certificate': next(c for c in certificates if c['id'] == first['id']) if first else None}
    (HERE/'evaluation.json').write_text(json.dumps(result, indent=2)+'\n')
    (HERE/'decisions.json').write_text(json.dumps(compared, indent=2)+'\n')
    compact = {k: v for k, v in result.items() if k not in ['first_primary_disagreement_ordinary','first_primary_disagreement_certificate','methods']}
    compact['methods'] = {m: {k: v for k, v in r.items() if k != 'differences'} for m, r in summary.items()}
    print(json.dumps(compact, indent=2))


if __name__ == '__main__':
    main()
