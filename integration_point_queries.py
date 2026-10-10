"""Ordinary midpoint/FFT comparison. It provides no certified error bound."""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
from scipy.signal import fftconvolve
from integration_bounds import ROOT
from integration_queries import QUERY_CASES, query_model, event_definitions, validate_restrictions
from integration_scaling import describe,save


def cost(potential,x):
    result=np.zeros_like(x,dtype=float)
    for target,above,below in potential.terms:
        d=x-float(target)
        result+=np.where(d>=0,float(above),float(below))*d*d
    return result


def point_log_integral(node,h,hi,restrictions,diagnostics):
    validate_restrictions(node,restrictions)
    count=int(hi/h)
    if h<=0 or hi/h!=count:
        raise ValueError('positive aligned mesh required')
    def normalize(masses,scale):
        total=float(masses.sum())
        if not math.isfinite(total) or total<=0:
            raise ArithmeticError('invalid point-estimate mass')
        return masses/total,scale+math.log(total)
    def visit(current):
        if not current.children:
            left=np.arange(count,dtype=float)*h
            center=left+h/2
            a,b=restrictions.get(current.name,(F(0),F(hi)))
            b=hi if b is None else float(b)
            lengths=np.maximum(0,np.minimum(left+h,b)-np.maximum(left,float(a)))
            masses=lengths*np.exp(-cost(current.potential,center))
            masses,scale=normalize(masses,0.0)
            return masses,scale,1
        children=[visit(c) for c in current.children]
        while len(children)>1:
            combined=[]
            for i in range(0,len(children),2):
                if i+1==len(children):combined.append(children[i]);continue
                a,sa,na=children[i];b,sb,nb=children[i+1]
                masses=fftconvolve(a,b)
                negative=masses[masses<0]
                diagnostics['negative_fft_coefficients']+=int(negative.size)
                diagnostics['clipped_negative_mass_sum']+=float(-negative.sum())
                masses=np.maximum(masses,0)
                masses,scale=normalize(masses,sa+sb)
                combined.append((masses,scale,na+nb))
            children=combined
        masses,scale,n=children[0]
        centers=(np.arange(masses.size,dtype=float)+n/2)*h
        masses=masses*np.exp(-cost(current.potential,centers))
        masses,scale=normalize(masses,scale)
        return masses,scale,n
    _,scale,_=visit(node)
    return scale


def worker(case,path):
    node=query_model(case);events=event_definitions(node);tree=describe(node)
    result={'case':case,'status':'RUNNING','attempts':[],
            'tree_sha256':hashlib.sha256(json.dumps(tree,sort_keys=True).encode()).hexdigest(),
            'method':'ordinary midpoint convolution with FFT; partial-event cell fractions',
            'guarantee':'NONE; cube-only estimates; clipping of negative roundoff coefficients recorded'}
    started=time.perf_counter();path.parent.mkdir(parents=True,exist_ok=True);save(path,result)
    previous=None
    for j in range(11):
        h=20/(2**j);before=time.perf_counter();logs={};diagnostics={'negative_fft_coefficients':0,'clipped_negative_mass_sum':0.0}
        for label,event in events.items():logs[label]=point_log_integral(node,h,800.0,event,diagnostics)
        probs={'P_A':math.exp(logs['A']-logs['Z']),'P_B':math.exp(logs['B']-logs['Z']),
               'P_AB':math.exp(logs['AB']-logs['Z']),'P_A_given_B':math.exp(logs['AB']-logs['B'])}
        probs['Delta']=probs['P_A_given_B']-probs['P_A']
        row={'mesh_ms':str(F(20,2**j)),'log_integrals':logs,'probabilities':probs,
             'successive_probability_changes':None if previous is None else {k:probs[k]-previous[k] for k in probs},
             'seconds':time.perf_counter()-before,'cumulative_worker_seconds':time.perf_counter()-started,
             'peak_resident_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,**diagnostics}
        result['attempts'].append(row);previous=probs;save(path,result)
    result['status']='ALL_ELEVEN_MESHES_FINISHED';save(path,result)


def main(output):
    output.mkdir(parents=True,exist_ok=False)
    result={'status':'RUNNING','started_at':datetime.now(timezone.utc).isoformat(),'cases':[],
            'process_budget_seconds_per_case':180,
            'driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'guarantee':'NONE; ordinary estimates are not a validated-quadrature competitor'}
    for case in QUERY_CASES:
        path=(output/(case+'.json')).resolve();before=time.perf_counter();timeout=False
        try:
            cp=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker',case,'--output',str(path)],timeout=180)
            code=cp.returncode
        except subprocess.TimeoutExpired:timeout=True;code=None
        record=json.loads(path.read_text()) if path.is_file() else {'case':case,'attempts':[]}
        elapsed=time.perf_counter()-before;record.update(process_elapsed_seconds=elapsed,process_exit_code=code)
        if timeout:record['status']='TIME_BUDGET_EXCEEDED'
        elif code:record['status']='PROCESS_ERROR'
        save(path,record);result['cases'].append({'case':case,'status':record['status'],'seconds':elapsed,'result':path.name})
        save(output/'results.json',result)
    result.update(status='ALL_PLANNED_CASES_HANDLED',completed_at=datetime.now(timezone.utc).isoformat())
    save(output/'results.json',result)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--worker',choices=QUERY_CASES);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();worker(a.worker,a.output) if a.worker else main(a.output)
