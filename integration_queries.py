"""Event-integral enclosures for continuous duration models.

Product-interval events retain the leaf-grid geometry. Four full-domain
integrals yield marginal, joint, conditional and difference probability bounds.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

from flint import arb, ctx
from integration_bounds import ROOT, ball, fraction_string, hayes_models, tail_envelope
from integration_moments import (MomentMessage, combine, interval, leaf_moments,
                                 multiply_factor)
from integration_scaling import CASES, build_case, describe, save

QUERY_CASES = ['HH0H', 'HHLH', 'HLLH', *(c[0] for c in CASES)]


def query_model(case):
    return hayes_models()[case] if case in ('HH0H', 'HHLH', 'HLLH') else build_case(case)


def validate_restrictions(node, restrictions):
    names = [leaf.name for leaf in node.leaves]
    if len(names) != len(set(names)):
        raise ValueError('query leaf names must be unique')
    if set(restrictions) - set(names):
        raise ValueError('unknown restricted leaf')
    for lower, upper in restrictions.values():
        if upper is not None and upper < lower:
            raise ValueError('reversed event interval')


def event_message(node, h, lo, hi, restrictions, cache):
    h, lo, hi = F(h), F(lo), F(hi)
    if not node.children:
        restriction = restrictions.get(node.name, (lo, hi))
        key = (node.potential, h, lo, hi, restriction)
        if key not in cache:
            cells = (hi-lo)/h
            if h <= 0 or cells.denominator != 1 or cells <= 0:
                raise ValueError('positive aligned mesh required')
            lower, upper = restriction
            lower, upper = F(lower), hi if upper is None else F(upper)
            values = []
            for k in range(int(cells)):
                left = max(lo+k*h, lower)
                right = min(lo+(k+1)*h, upper)
                center = lo+(F(k)+F(1,2))*h
                values.append(leaf_moments(node.potential,left,right,center)
                              if right > left else [arb(0),arb(0),arb(0)])
            cache[key] = MomentMessage([list(x) for x in zip(*values)],1,lo)
        return cache[key]
    children = [event_message(c,h,lo,hi,restrictions,cache) for c in node.children]
    while len(children) > 1:
        children = [combine(children[i],children[i+1]) if i+1<len(children) else children[i]
                    for i in range(0,len(children),2)]
    message = children[0]
    radius = message.leaves*h/2
    result = [[],[],[]]
    for k,moments in enumerate(zip(*message.moments)):
        center = message.origin+k*h+radius
        values = multiply_factor(moments,node.potential,center,radius)
        for row,value in zip(result,values):
            row.append(value)
    return MomentMessage(result,message.leaves,message.origin)


def event_integral(node,h,hi,restrictions,tail=None,cache=None):
    h, hi = F(h), F(hi)
    validate_restrictions(node,restrictions)
    message = event_message(node,h,F(0),hi,restrictions,{} if cache is None else cache)
    mass = sum(message.moments[0],arb(0))
    tail = tail_envelope(node,hi) if tail is None else tail
    return mass.lower().max(0),(mass.upper()+tail).upper()


def event_definitions(node):
    a,b = node.leaves[:2]
    return {'Z':{},'A':{a.name:(a.potential.terms[0][0],None)},
            'B':{b.name:(b.potential.terms[0][0],None)},
            'AB':{a.name:(a.potential.terms[0][0],None),
                  b.name:(b.potential.terms[0][0],None)}}


def probability_results(integrals):
    def ratio(numerator,denominator):
        l,u=integrals[numerator];dl,du=integrals[denominator]
        if not dl > 0:
            raise ArithmeticError('conditioning lower mass must be positive')
        return interval((l/du).lower(),(u/dl).upper()).intersection(interval(arb(0),arb(1)))
    pa=ratio('A','Z');pb=ratio('B','Z');pab=ratio('AB','Z');conditional=ratio('AB','B')
    return {'P_A':pa,'P_B':pb,'P_AB':pab,'P_A_given_B':conditional,
            'Delta':conditional-pa}


def encode_interval(value):
    lo,hi=value.lower(),value.upper()
    return {'lower':fraction_string(lo),'upper':fraction_string(hi),
            'lower_approx':float(lo),'upper_approx':float(hi),
            'width_upper':fraction_string((hi-lo).upper())}


def worker(case,path):
    ctx.prec=128
    started=time.perf_counter();node=query_model(case);events=event_definitions(node)
    tree=describe(node);tail=tail_envelope(node,F(800))
    result={'case':case, 'status':'RUNNING', 'precision_bits':ctx.prec, 'attempts':[], 'started_at':datetime.now(timezone.utc).isoformat(), 'tree':tree, 'events':{label:{k:[str(v[0]),None if v[1] is None else str(v[1])]
                             for k,v in restriction.items()}
                      for label,restriction in events.items()}, 'tail_upper_exact_rational':fraction_string(tail), 'rounding_and_scope':'Arb endpoints under audited mathematical assumptions; not full implementation verification'}
    path.parent.mkdir(parents=True,exist_ok=True);save(path,result)
    for j in range(15):
        h=F(20,2**j);result['in_progress_mesh_ms']=str(h);save(path,result)
        before=time.perf_counter();cache={};integrals={}
        for label,restriction in events.items():
            result['in_progress_event']=label;save(path,result)
            integrals[label]=event_integral(node,h,F(800),restriction,tail,cache)
        probs=probability_results(integrals)
        widths={k:(v.upper()-v.lower()).upper() for k,v in probs.items()}
        reached=all(w<=ball(F(1,500)) for w in widths.values())
        delta=probs['Delta']
        sign='NEGATIVE' if delta.upper()<0 else 'POSITIVE' if delta.lower()>0 else 'UNRESOLVED'
        row={'mesh_ms':str(h),'seconds':time.perf_counter()-before,
             'cumulative_worker_seconds':time.perf_counter()-started,
             'integrals':{k:{'lower':fraction_string(v[0]),'upper':fraction_string(v[1])}
                          for k,v in integrals.items()},
             'probabilities':{k:encode_interval(v) for k,v in probs.items()},
             'all_widths_at_most_0_002':reached,'delta_sign':sign,
             'peak_resident_bytes_after_completed_mesh':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        result['attempts'].append(row)
        del result['in_progress_mesh_ms'];del result['in_progress_event']
        if reached:result['status']='QUERY_TOLERANCE_REACHED_PENDING_BUDGET_CHECK'
        save(path,result)
        print(json.dumps({'case':case,'mesh_ms':str(h),'max_width':max(float(w) for w in widths.values()),
                          'delta_sign':sign,'seconds':row['cumulative_worker_seconds']}),flush=True)
        if reached:return
    result['status']='MAXIMUM_MESH_WITHOUT_TARGET';save(path,result)


def main(output):
    output.mkdir(parents=True,exist_ok=False)
    sources=['integration_queries.py','integration_moments.py','integration_bounds.py','integration_scaling.py']
    result={'status':'RUNNING', 'started_at':datetime.now(timezone.utc).isoformat(), 'cases':[], 'process_budget_seconds_per_case':180, 'interpretation':'Conditional predictions of fixed models; not empirical observations'}
    for case in QUERY_CASES:
        path=(output/(case+'.json')).resolve();before=time.perf_counter();timeout=False
        try:
            cp=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker',case,'--output',str(path)],
                              timeout=180,check=False)
            code=cp.returncode
        except subprocess.TimeoutExpired:
            timeout=True;code=None
        elapsed=time.perf_counter()-before
        record=json.loads(path.read_text()) if path.is_file() else {'case':case,'attempts':[]}
        record.update(process_elapsed_seconds=elapsed,process_exit_code=code)
        if timeout:record['status']='TIME_BUDGET_EXCEEDED'
        elif code:record['status']='PROCESS_ERROR'
        elif record['status']=='QUERY_TOLERANCE_REACHED_PENDING_BUDGET_CHECK':
            record['status']='QUERY_TOLERANCE_MET' if elapsed<=180 else 'TIME_BUDGET_EXCEEDED'
        save(path,record)
        result['cases'].append({'case':case,'status':record['status'],'seconds':elapsed,'result':path.name})
        save(output/'results.json',result)
    result.update(status='ALL_PLANNED_CASES_HANDLED',completed_at=datetime.now(timezone.utc).isoformat())
    save(output/'results.json',result)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--worker',choices=QUERY_CASES);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();worker(a.worker,a.output) if a.worker else main(a.output)
