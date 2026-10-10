"""Run the fixed moment-enclosure feasibility experiment with hard timeouts."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time

import flint
from flint import ctx
from integration_bounds import ROOT,ball,fraction_string
from integration_moments import partition_enclosure_moments
from integration_scaling import CASES,build_case,describe,save


def worker(case,path):
    ctx.prec=128
    node=build_case(case);tree=describe(node);started=time.perf_counter()
    result={'case':case,'status':'RUNNING','dimension':len(node.leaves),
            'tree':tree,'tree_sha256':hashlib.sha256(json.dumps(tree,sort_keys=True).encode()).hexdigest(),
            'attempts':[],'started_at':datetime.now(timezone.utc).isoformat(),
            'method':'centered_moment_balls','precision_bits':ctx.prec,'python_flint':flint.__version__,
            'source_sha256':hashlib.sha256(Path(__file__).with_name('integration_moments.py').read_bytes()).hexdigest(),
            'memory_note':'macOS ru_maxrss bytes after completed meshes; timeout record is only a lower bound on actual process peak'}
    path.parent.mkdir(parents=True,exist_ok=True);save(path,result)
    for j in range(15):
        h=F(20,2**j);result['in_progress_mesh_ms']=str(h);save(path,result)
        before=time.perf_counter()
        lo,up,tail,width,count=partition_enclosure_moments(node,h,F(800))
        reached=bool(width.is_finite() and width<=ball(F(1,1000)))
        row={'mesh_ms':str(h),'seconds':time.perf_counter()-before,
             'cumulative_worker_seconds':time.perf_counter()-started,
             'Z_lower_exact_rational':fraction_string(lo),'Z_upper_exact_rational':fraction_string(up),
             'tail_upper_exact_rational':fraction_string(tail),
             'log_width_upper_exact_rational':fraction_string(width) if width.is_finite() else None,
             'log_width_approx':float(width) if width.is_finite() else None,
             'target_met_by_implemented_bounds':reached,'coefficients_at_root':count,
             'peak_resident_bytes_after_completed_mesh':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        result['attempts'].append(row);del result['in_progress_mesh_ms']
        if reached:result['status']='TOLERANCE_REACHED_PENDING_PROCESS_BUDGET_CHECK'
        save(path,result)
        print(json.dumps({'case':case,'mesh_ms':str(h),'log_width':row['log_width_approx'],
                          'seconds':row['cumulative_worker_seconds']}),flush=True)
        if reached:return
    result['status']='MAXIMUM_MESH_WITHOUT_TARGET';save(path,result)


def main(output):
    output.mkdir(parents=True,exist_ok=False)
    result={'started_at':datetime.now(timezone.utc).isoformat(),
            'status':'RUNNING','cases':[],'process_budget_seconds_per_case':60,
            'driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'platform':platform.platform(),'machine':platform.machine(),
            'interpretation':'Measured normalizer accuracy and elapsed time for the specified model and implementation'}
    for case,*_ in CASES:
        path=(output/(case+'.json')).resolve();before=time.perf_counter();timeout=False
        try:
            cp=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker',case,
                               '--output',str(path)],timeout=60,check=False)
            code=cp.returncode
        except subprocess.TimeoutExpired:
            timeout=True;code=None
        elapsed=time.perf_counter()-before
        record=json.loads(path.read_text()) if path.is_file() else {'case':case,'attempts':[]}
        record.update(process_elapsed_seconds=elapsed,process_exit_code=code)
        if timeout:record['status']='TIME_BUDGET_EXCEEDED'
        elif code:record['status']='PROCESS_ERROR'
        elif record['status']=='TOLERANCE_REACHED_PENDING_PROCESS_BUDGET_CHECK':
            record['status']='FEASIBILITY_TARGET_MET' if elapsed<=60 else 'TIME_BUDGET_EXCEEDED'
        save(path,record)
        result['cases'].append({'case':case,'status':record['status'],
                                'seconds':elapsed,'result':path.name});save(output/'results.json',result)
    result.update(status='ALL_PLANNED_CASES_HANDLED',completed_at=datetime.now(timezone.utc).isoformat())
    save(output/'results.json',result)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--worker',choices=[c[0] for c in CASES]);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    worker(a.worker,a.output) if a.worker else main(a.output)
