"""Run synthetic nested-tree normalizer calculations with hard case budgets.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import json
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time

import flint
from flint import ctx

from integration_bounds import (ROOT, Node, Potential, ball, fraction_string,
                                hayes_models, partition_enclosure)

CASES = [('balanced_8',8,False,1), ('balanced_16',16,False,1),
         ('left_branching_16',16,True,1), ('stiff_balanced_16',16,False,4)]


def build_case(name):
    _,n,left_branching,stiffness = next(c for c in CASES if c[0]==name)
    source=hayes_models()['HLLH']
    heavy,light=source.leaves[:2]
    target,above,below=source.potential.terms[0]
    nodes=[(Node(heavy.potential if i%2==0 else light.potential,name=f'x{i+1}'),
            2 if i%2==0 else 1) for i in range(n)]
    def join(a,b):
        m=a[1]+b[1]
        pot=Potential(((target*m/4,above*F(4,m)*stiffness,below*F(4,m)*stiffness),))
        return Node(pot,(a[0],b[0]),name=f'M{m}'),m
    if left_branching:
        root=nodes[0]
        for other in nodes[1:]:
            root=join(root,other)
    else:
        while len(nodes)>1:
            nodes=[join(nodes[i],nodes[i+1]) for i in range(0,len(nodes),2)]
        root=nodes[0]
    return root[0]


def describe(node):
    return {'name':node.name,'terms':[[str(x) for x in t] for t in node.potential.terms],
            'children':[describe(c) for c in node.children]}


def save(path,record):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(record,indent=2)+'\n')
    temporary.replace(path)


def worker(case,path,method):
    if method=='grid_recurrence':
        from integration_grid_factors import partition_enclosure_grid
        enclosure=partition_enclosure_grid
    else:
        enclosure=partition_enclosure
    ctx.prec=128
    node=build_case(case)
    description=describe(node)
    started=time.perf_counter()
    r={'case':case, 'status':'RUNNING', 'method':method, 'dimension':len(node.leaves), 'tree':description, 'attempts':[], 'started_at':datetime.now(timezone.utc).isoformat(), 'python_flint':flint.__version__, 'arithmetic_precision_bits':ctx.prec}
    path.parent.mkdir(parents=True,exist_ok=True)
    save(path,r)
    for j in range(15):
        h=F(20,2**j)
        r['in_progress_mesh_ms']=str(h)
        save(path,r)
        before=time.perf_counter()
        lo,up,tail,width,ncoeff=enclosure(node,h,F(800))
        reached=bool(width<=ball(F(1,1000)))
        row={'mesh_ms':str(h),'seconds':time.perf_counter()-before,
             'cumulative_worker_seconds':time.perf_counter()-started,
             'Z_lower_exact_rational':fraction_string(lo),'Z_upper_exact_rational':fraction_string(up),
             'tail_upper_exact_rational':fraction_string(tail),
             'log_width_upper_exact_rational':fraction_string(width),
             'log_width_approx':float(width),'target_met_by_implemented_bounds':reached,
             'coefficients_at_root':ncoeff,'peak_resident_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        r['attempts'].append(row)
        del r['in_progress_mesh_ms']
        if reached:
            r['status']='TOLERANCE_REACHED_PENDING_PROCESS_BUDGET_CHECK'
        save(path,r)
        print(json.dumps({'case':case,'mesh_ms':str(h),'log_width':float(width),
                          'cumulative_seconds':row['cumulative_worker_seconds']}),flush=True)
        if reached:
            return
    r['status']='MAXIMUM_MESH_WITHOUT_TARGET'
    save(path,r)


def main(output,method):
    output.mkdir(parents=True,exist_ok=False)
    results={'started_at':datetime.now(timezone.utc).isoformat(), 'method':method, 'machine':platform.machine(), 'platform':platform.platform(), 'seconds_budget_per_case':60, 'status':'RUNNING', 'cases':[], 'interpretation':'Measured normalizer accuracy and elapsed time for synthetic models'}
    for case,*_ in CASES:
        path=(output/(case+'.json')).resolve()
        before=time.perf_counter()
        timed_out=False
        try:
            proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker',case,
                                 '--output',str(path),'--method',method],timeout=60,check=False)
            exit_code=proc.returncode
        except subprocess.TimeoutExpired:
            timed_out=True
            exit_code=None
        elapsed=time.perf_counter()-before
        r=json.loads(path.read_text()) if path.is_file() else {'case':case,'attempts':[]}
        r['process_elapsed_seconds']=elapsed
        r['process_exit_code']=exit_code
        if timed_out:
            r['status']='TIME_BUDGET_EXCEEDED'
        elif exit_code:
            r['status']='PROCESS_ERROR'
        elif r.get('status')=='TOLERANCE_REACHED_PENDING_PROCESS_BUDGET_CHECK':
            r['status']='FEASIBILITY_TARGET_MET' if elapsed<=60 else 'TIME_BUDGET_EXCEEDED'
        save(path,r)
        results['cases'].append({'case':case,'status':r['status'],
                                'process_elapsed_seconds':elapsed,'result':path.name})
        save(output/'results.json',results)
    results['status']='COMPLETED_WITH_ALL_PLANNED_CASES'
    results['completed_at']=datetime.now(timezone.utc).isoformat()
    save(output/'results.json',results)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--worker',choices=[c[0] for c in CASES])
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--method',choices=['scalar_factors','grid_recurrence'],default='scalar_factors')
    args=parser.parse_args()
    if args.worker:
        worker(args.worker,args.output,args.method)
    else:
        main(args.output,args.method)
