"""Independent checks of the two illustrated speech-planning decisions."""
from fractions import Fraction as F
from pathlib import Path
import json
import mpmath as mp
import numpy as np

from run import HERE, inputs, quadratic_cdf


def real(value):
    value = F(value)
    return mp.mpf(value.numerator)/value.denominator


def high_precision_cdf(row, deadline, digits):
    mp.mp.dps = digits
    a, b = [2*F(w) for w in row['unary_weight_exact']]
    c = 2*F(row['group_weight_exact'])
    determinant = a*b+c*(a+b)
    v00, v01, v11 = (b+c)/determinant, -c/determinant, (a+c)/determinant
    mu = list(map(F, row['target_ms_exact']))
    h0, h1 = a*mu[0]+c*F(row['total_ms']), b*mu[1]+c*F(row['total_ms'])
    m0, m1 = real(v00*h0+v01*h1), real(v01*h0+v11*h1)
    scale = mp.sqrt(real(v00))
    conditional_scale = mp.sqrt(real(v11-v01*v01/v00))
    z0 = -m0/scale
    def f(z):
        z2 = (m1+real(v01)/scale*z)/conditional_scale
        return mp.exp(-z*z/2)/mp.sqrt(2*mp.pi)*mp.erfc(-z2/mp.sqrt(2))/2
    def knots(a, b):
        return [a]+[mp.mpf(x) for x in [-8,-4,0,4,8] if a<x<b]+[b]
    denominator = mp.quad(f, knots(z0, mp.inf))
    return mp.quad(f, knots(z0, (real(deadline)-m0)/scale))/denominator


def dense_midpoint(row, deadline, mesh):
    mu = np.array([float(F(v)) for v in row['target_ms_exact']])
    w = np.array([float(F(v)) for v in row['unary_weight_exact']])
    shared = float(F(row['group_weight_exact']))
    total = float(F(row['total_ms']))
    left = np.arange(round(1600/mesh))*mesh
    x = left+mesh/2
    density = np.exp(-w[0]*(x[:,None]-mu[0])**2-w[1]*(x[None,:]-mu[1])**2
                     -shared*(x[:,None]+x[None,:]-total)**2)
    fraction = np.clip((deadline-left)/mesh, 0., 1.)
    return float((density*fraction[:,None]).sum()/density.sum())


def main():
    rows = {r['id']: r for r in inputs()}
    ordinary = {r['id']: r for r in json.loads((HERE/'ordinary.json').read_text())['rows']}
    result = {'status': 'PASS', 'cases': []}
    for identifier, deadlines in [('slt/arctic_a0120/8',[340,350]),('bdl/arctic_b0340/2',[240])]:
        row = rows[identifier]
        certificate = json.loads((HERE/'certificates'/(identifier.replace('/','__')+'.json')).read_text())
        case = {'id': identifier, 'word': row['word'], 'probabilities': {}}
        cdf, _ = quadratic_cdf(row)
        for deadline in deadlines:
            p50 = high_precision_cdf(row, deadline, 50)
            p80 = high_precision_cdf(row, deadline, 80)
            assert abs(p50-p80) < mp.mpf('1e-45')
            assert abs(mp.mpf(cdf(deadline))-p80) < mp.mpf('1e-12')
            bounds = certificate['attempts'][-1]['cdf'][str(deadline)]
            assert real(bounds['lower']) <= p80 <= real(bounds['upper'])
            case['probabilities'][str(deadline)] = mp.nstr(p80,70)
        comparisons = {}
        for h in [10.,5.]:
            d = deadlines[0]
            direct = dense_midpoint(row,d,h)
            fft = ordinary[identifier]['midpoint_fft'][str(h)]['all_cdf_queries'][str(d)]
            assert abs(direct-fft)<1e-13
            comparisons[str(h)] = {'dense_midpoint': direct, 'fft': fft, 'absolute_difference': abs(direct-fft)}
        case['independent_dense_midpoint_checks'] = comparisons
        result['cases'].append(case)
    (HERE/'independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
