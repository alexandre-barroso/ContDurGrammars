"""Independent numerical controls of moment enclosures; not a formal proof."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import time

import mpmath as mp
import numpy as np
from flint import arb, ctx
from integration_bounds import Node, Potential, fraction_string
from integration_moments import (leaf_moments, factor_derivative_bounds,
                                 enclose_moments, partition_enclosure_moments)


def mm(x):
    x=F(x); return mp.mpf(x.numerator)/x.denominator


def bounds(x):
    return mm(F(fraction_string(x.lower()))), mm(F(fraction_string(x.upper())))


def own_cost(terms,x):
    return sum((mm(a if x>=mm(t) else b)*(x-mm(t))**2 for t,a,b in terms),mp.mpf(0))


def own_value(terms,x):
    return mp.exp(-own_cost(terms,x))


def run_leaf_controls():
    cases = [
        ((),F(-3),F(5)),
        (((F(1,3),F(1,2),F(2,3)),),F(-2),F(3)),
        (((F(-4,5),F(0),F(3)),(F(6,5),F(2),F(0))),F(-2),F(3)),
        (((F(-1),F(3),F(1)),(F(1,2),F(1,8),F(2)),(F(7,5),F(0),F(7,9))),F(-1),F(2)),
        (((F(1000001),F(2,7),F(3,7)),),F(1000000),F(1000002)),
        (((F(1,2),F(7),F(11)),),F(499,1000),F(501,1000))]
    count=0
    for terms,lo,hi in cases:
        pot=Potential(terms); center=(lo+hi)/2
        actual=leaf_moments(pot,lo,hi,center)
        cuts=[mm(lo),*(mm(t) for t in sorted({t for t,_,_ in terms}) if lo<t<hi),mm(hi)]
        for j,enclosure in enumerate(actual):
            reference=mp.quad(lambda x:(x-mm(center))**j*own_value(terms,x),cuts)
            low,high=bounds(enclosure)
            # At exact zero, quadrature can have an 80-digit cancellation
            # residual. Record a tiny independent quadrature tolerance only.
            slack=mp.mpf('1e-65')*max(1,abs(reference))
            assert low-slack<=reference<=high+slack,(terms,j,low,reference,high)
            count+=1
    return {'cases':len(cases),'moment_checks':count,'quadrature_digits':mp.mp.dps,
            'reference_tolerance':'1e-65*max(1,abs(reference)); not an integration proof'}


def run_derivative_controls():
    terms=((F(-7,10),F(3,7),F(4,5)),(F(1,3),F(0),F(2)),
           (F(8,5),F(5,2),F(1,9)))
    pot=Potential(terms); count=0
    for lo,hi in [(F(-2),F(3)),(F(-7,10),F(1,3)),(F(1,3),F(8,5)),
                  (F(1599,1000),F(1601,1000))]:
        _,_,k1,k2=factor_derivative_bounds(pot,lo,hi,(lo+hi)/2)
        cap1,cap2=bounds(k1)[1],bounds(k2)[1]
        points={lo+(hi-lo)*F(i,100) for i in range(101)}
        points.update(t for t,_,_ in terms if lo<=t<=hi)
        for q in points:
            x=mm(q); value=own_value(terms,x)
            for right_branch in (True,False):
                weights=[a if (q>=t if right_branch else q>t) else b for t,a,b in terms]
                d1=sum((2*mm(w)*(x-mm(t)) for (t,_,_),w in zip(terms,weights)),mp.mpf(0))
                d2=2*sum((mm(w) for w in weights),mp.mpf(0))
                # At the outer endpoint, the branch outside the cell does not
                # contribute a Lipschitz bound inside the interval.
                if (q==lo and not right_branch) or (q==hi and right_branch):
                    continue
                assert abs(-d1*value)<=cap1+mp.mpf('1e-65')
                assert abs((d1*d1-d2)*value)<=cap2+mp.mpf('1e-65')
                count+=2
    trap=Potential(((F(2),F(1),F(1)),))
    _,_,_,cap=factor_derivative_bounds(trap,F(1),F(3),F(2))
    assert bounds(cap)[1]>=2  # |g''(2)|=2 exceeds the endpoint value 2/e.
    return {'derivative_checks':count,'interior_curvature_counterexample':'checked',
            'scope':'finite samples and knot one-sided limits'}


def nested_model():
    p=Potential(((F(6,5),F(1,7),F(1,5)),))
    q=Potential(((F(5,2),F(1,9),F(2,7)),))
    a=Node(Potential(((F(17,5),F(1,31),F(1,19)),)),(Node(p),Node(q)))
    b=Node(Potential(((F(22,5),F(1,37),F(1,23)),)),(Node(q),Node(p)))
    return Node(Potential(((F(29,4),F(1,91),F(1,57)),)),(a,b))


def quadrature_cells(order):
    """Independent tensor quadrature; does not use messages or production costs."""
    points,weights=np.polynomial.legendre.leggauss(order)
    points=(points+1)/2;weights=weights/2
    grid=np.array(list(itertools.product(range(order),repeat=4)),dtype=int)
    offsets=points[grid];w=np.prod(weights[grid],axis=1)
    totals=np.zeros((3,13))
    def penalty(x,t,above,below):
        return np.where(x>=t,above,below)*(x-t)**2
    for cells in itertools.product(range(4),repeat=4):
        x=offsets+np.array(cells)
        cost=(penalty(x[:,0],1.2,1/7,1/5)+penalty(x[:,1],2.5,1/9,2/7)
              +penalty(x[:,2],2.5,1/9,2/7)+penalty(x[:,3],1.2,1/7,1/5)
              +penalty(x[:,0]+x[:,1],3.4,1/31,1/19)
              +penalty(x[:,2]+x[:,3],4.4,1/37,1/23)
              +penalty(x.sum(axis=1),7.25,1/91,1/57))
        k=sum(cells);d=x.sum(axis=1)-(k+2)
        mass=w*np.exp(-cost)
        for j in range(3): totals[j,k]+=np.dot(mass,d**j)
    return totals


def run_nested_controls():
    model=nested_model(); message=enclose_moments(model,F(1),F(0),F(4),{})
    q12,q20=quadrature_cells(12),quadrature_cells(20)
    max_delta=float(np.max(np.abs(q12-q20)))
    assert max_delta<1e-4,max_delta
    checks=0
    for j,row in enumerate(message.moments):
        for k,x in enumerate(row):
            lo,hi=map(float,bounds(x))
            assert lo<=q12[j,k]<=hi and lo<=q20[j,k]<=hi,(j,k,lo,q12[j,k],q20[j,k],hi)
            checks+=1
    coarse=partition_enclosure_moments(model,F(1),F(4))
    fine=partition_enclosure_moments(model,F(1,2),F(4))
    # Taylor intervals need not be nested under refinement; verify overlap and
    # actual improvement here, without making a universal nesting claim.
    assert max(float(coarse[0]),float(fine[0]))<=min(float(coarse[1]),float(fine[1]))
    assert float(fine[3])<float(coarse[3])
    return {'cartesian_leaf_boxes':256,'checked_aggregate_moments':checks,
            'quadrature_orders':[12,20],'maximum_reference_difference':max_delta,
            'coarse_log_width':float(coarse[3]),'fine_log_width':float(fine[3]),
            'scope':'independent numerical quadrature checks; not certified oracle values'}


def main(output):
    before=time.perf_counter();ctx.prec=128;mp.mp.dps=80
    result={'started_at':datetime.now(timezone.utc).isoformat(),
            'status':'RUNNING','arithmetic_precision_bits':ctx.prec,
            'source_sha256':hashlib.sha256(Path(__file__).with_name('integration_moments.py').read_bytes()).hexdigest()}
    result['leaf_controls']=run_leaf_controls()
    result['derivative_controls']=run_derivative_controls()
    result['nested_controls']=run_nested_controls()
    result.update(status='ALL_PLANNED_NUMERICAL_CONTROLS_PASSED',seconds=time.perf_counter()-before)
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
