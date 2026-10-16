"""Numerical controls for restricted-cell queries, independently integrated."""
from __future__ import annotations

from datetime import datetime,timezone
from fractions import Fraction as F
import json
from pathlib import Path
import time

import mpmath as mp
from scipy.integrate import quad
from flint import arb,ctx
from integration_bounds import ROOT,Node,Potential,hayes_models
from integration_moments import leaf_moments
from integration_queries import event_message,event_integral,event_definitions,probability_results,encode_interval


def m(q):
    q=F(q);return mp.mpf(q.numerator)/q.denominator


def endpoints(a):
    return m(F(str(a.lower().fmpq()))),m(F(str(a.upper().fmpq())))


def direct_density(potential,x):
    return mp.exp(-sum(m(a if x>=m(t) else b)*(x-m(t))**2
                       for t,a,b in potential.terms))


def main(output):
    started=time.perf_counter();ctx.prec=128;mp.mp.dps=85;checks=[]
    potential=Potential(((F(3,5),F(1,3),F(7,4)),(F(7,10),F(2),F(1,5))))
    for left,right,center in [(F(3,10),F(1,2),F(1,2)),(F(11,20),F(9,10),F(1,2)),
                              (F(0),F(1,7),F(1,2)),(F(3,5),F(7,10),F(1,2))]:
        values=leaf_moments(potential,left,right,center)
        cuts=[m(left),*(m(t) for t in potential.knots if left<t<right),m(right)]
        for j,value in enumerate(values):
            truth=mp.quad(lambda x:(x-m(center))**j*direct_density(potential,x),cuts)
            lower,upper=endpoints(value)
            assert lower<=truth<=upper,(left,right,j,truth,lower,upper)
            checks.append({'kind':'partial_cell_moment','order':j,'left':str(left),'right':str(right),
                           'center':str(center),'reference':mp.nstr(truth,70),'inside':True})
    # Repeated equal potentials exercise event-aware caching in a nested tree.
    leaves=tuple(Node(potential,name=f'x{i}') for i in range(4))
    node=Node(Potential(),(Node(Potential(),leaves[:2]),Node(Potential(),leaves[2:])))
    events=[{}, {'x0':(F(3,7),None)}, {'x1':(F(1,3),F(5,6))},
            {'x0':(F(3,7),None),'x1':(F(1,3),F(5,6))}]
    cache={}
    for restriction in events:
        message=event_message(node,F(1,4),F(0),F(1),restriction,cache)
        value=sum(message.moments[0],arb(0));truth=mp.mpf(1)
        for leaf in leaves:
            a,b=restriction.get(leaf.name,(F(0),F(1)));b=F(1) if b is None else b
            cuts=[m(a),*(m(t) for t in potential.knots if a<t<b),m(b)]
            truth*=mp.quad(lambda x:direct_density(potential,x),cuts)
        lower,upper=endpoints(value);assert lower<=truth<=upper,(restriction,truth,lower,upper)
        checks.append({'kind':'nested_independent_event','restricted_leaves':list(restriction),
                       'reference':mp.nstr(truth,70),'inside':True})
    full=event_message(node,F(1,4),F(0),F(1),{},{}).moments[0]
    all_event=event_message(node,F(1,4),F(0),F(1),{'x0':(F(-1),F(2))},{}).moments[0]
    # Arb == asks for certain real equality, so two identical nonzero-radius
    # balls compare false. Compare their exact stored endpoints instead.
    assert all(a.lower().fmpq()==b.lower().fmpq() and a.upper().fmpq()==b.upper().fmpq()
               for a,b in zip(full,all_event))
    empty=event_message(node,F(1,4),F(0),F(1),{'x0':(F(2),None)},{}).moments[0]
    assert all(v==0 for v in empty)
    for bad in [{'unknown':(F(0),None)},{'x0':(F(1),F(0))}]:
        try:event_integral(node,F(1,4),F(1),bad)
        except ValueError:pass
        else:raise AssertionError('invalid event accepted')
    checks.append({'kind':'empty_full_unknown_reversed_event_controls','passed':True})

    # Independent adaptive quadrature for the actual two-variable case. The
    # returned QUADPACK error estimates are heuristics, not enclosing bounds.
    import math
    grammar=hayes_models()['HH0H'];definitions=event_definitions(grammar)
    def scalar_cost(p,x):
        return sum(float(a if x>=float(t) else b)*(x-float(t))**2 for t,a,b in p.terms)
    estimates={};reported={}
    for label,restriction in definitions.items():
        al=float(restriction.get(grammar.leaves[0].name,(F(0),None))[0]);bl=float(restriction.get(grammar.leaves[1].name,(F(0),None))[0])
        def outer(x):
            cuts=sorted({bl,800.0,*[float(t) for t in grammar.leaves[1].potential.knots if bl<float(t)<800],
                         *[float(t)-x for t in grammar.potential.knots if bl<float(t)-x<800]})
            return sum(quad(lambda y:math.exp(-scalar_cost(grammar.leaves[0].potential,x)
                                             -scalar_cost(grammar.leaves[1].potential,y)
                                             -scalar_cost(grammar.potential,x+y)),a,b,
                            epsabs=1e-10,epsrel=1e-10,limit=150)[0] for a,b in zip(cuts,cuts[1:]))
        cuts=sorted({al,800.0,*[float(t) for t in grammar.leaves[0].potential.knots if al<float(t)<800],
                     *[float(t)-bl for t in grammar.potential.knots if al<float(t)-bl<800]})
        values=[quad(outer,a,b,epsabs=1e-8,epsrel=1e-10,limit=150) for a,b in zip(cuts,cuts[1:])]
        estimates[label]=sum(v[0] for v in values);reported[label]=sum(v[1] for v in values)
    integrals={label:event_integral(grammar,F(5),F(800),event) for label,event in definitions.items()}
    for label,(lo,hi) in integrals.items():
        assert float(lo)<=estimates[label]<=float(hi),(label,estimates[label],lo,hi)
    truth_probs={'P_A':estimates['A']/estimates['Z'],'P_B':estimates['B']/estimates['Z'],
                 'P_AB':estimates['AB']/estimates['Z'],'P_A_given_B':estimates['AB']/estimates['B']}
    truth_probs['Delta']=truth_probs['P_A_given_B']-truth_probs['P_A']
    bounds=probability_results(integrals)
    for name,value in bounds.items():assert float(value.lower())<=truth_probs[name]<=float(value.upper())
    result={'completed_at':datetime.now(timezone.utc).isoformat(), 'status':'ALL_CONTROLS_PASSED', 'seconds':time.perf_counter()-started, 'controls':checks, 'HH0H_independent_QUADPACK':{'cube_integrals':estimates,'outer_reported_error_estimates':reported,
                                       'probabilities':truth_probs,'enclosures':{k:encode_interval(v) for k,v in bounds.items()}}, 'limitations':'Finite numerical controls, not independent certificates or a whole-program proof'}
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'seconds':result['seconds'],'partial_moments':12,'independent_event_integrals':4}))


if __name__=='__main__':
    import argparse
    import sys
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if sys.flags.optimize:
        parser.error('Do not disable assertion-based controls with -O')
    if args.output.exists():
        parser.error('Output already exists; choose a new file')
    main(args.output)

