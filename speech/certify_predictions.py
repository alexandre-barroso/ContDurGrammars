"""Compute enclosed midpoint probabilities with the duration-grammar integrator.

Run after prepare_data.py and fit_models.py. Each word gets an audit JSON with
all mesh attempts, exact enclosing endpoints, and the supplied model inputs.
"""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import argparse
import json
import platform
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from flint import ctx
import flint
from integration_bounds import Node, Potential, ball, fraction_string, tail_envelope
from integration_queries import event_integral
from integration_scaling import describe, save


def result_path(identifier):
    return HERE / 'certificates' / (identifier.replace('/', '__') + '.json')


def model_tree(row):
    leaves = tuple(Node(Potential(((F(m), F(w), F(w)),)), name=f'x{i+1}')
                   for i, (m,w) in enumerate(zip(row['target_ms_exact'], row['unary_weight_exact'])))
    return Node(Potential(((F(row['total_ms']), F(row['group_weight_exact']), F(row['group_weight_exact'])),)), leaves, name='word')


def certify(row):
    ctx.prec = 128
    path = result_path(row['id'])
    root = model_tree(row)
    if path.exists():
        prior = json.loads(path.read_text())
        if (prior.get('status') == 'TOLERANCE_MET' and prior.get('tree') == describe(root)
                and prior.get('event') == {'x1_lower_ms':str(F(row['total_ms'])/2),'x1_upper_ms':None}
                and prior.get('precision_bits') == 128 and prior.get('cube_upper_ms') == '1600'):
            return {k: prior[k] for k in ['id','status','seconds','attempt_count']}
    started = time.perf_counter()
    cube = F(1600)
    tail = tail_envelope(root,cube)
    record = {'id':row['id'], 'status':'RUNNING', 'precision_bits':ctx.prec,
              'cube_upper_ms':str(cube),
              'tree':describe(root), 'event':{'x1_lower_ms':str(F(row['total_ms'])/2),'x1_upper_ms':None},
              'tail_upper_exact':fraction_string(tail), 'attempts':[]}
    path.parent.mkdir(exist_ok=True)
    save(path,record)
    try:
        for j in range(11):
            mesh = F(20,2**j)
            before = time.perf_counter()
            cache = {}
            zlo,zhi = event_integral(root,mesh,cube,{},tail,cache)
            alo,ahi = event_integral(root,mesh,cube,{'x1':(F(row['total_ms'])/2,None)},tail,cache)
            attempt = {'mesh_ms':str(mesh), 'seconds':time.perf_counter()-before,
                       'Z':{'lower':fraction_string(zlo),'upper':fraction_string(zhi)},
                       'A':{'lower':fraction_string(alo),'upper':fraction_string(ahi)}}
            if zlo > 0:
                raw_lower = F(fraction_string((alo/zhi).lower()))
                raw_upper = F(fraction_string((ahi/zlo).upper()))
                # Intersect serialized exact endpoints with the probability
                # range. Reboxing in an Arb ball could round outside [0,1].
                lower,upper=max(F(0),raw_lower),min(F(1),raw_upper)
                if lower>upper:raise ArithmeticError('Empty probability enclosure')
                attempt['probability'] = {'lower':str(lower),'upper':str(upper),
                                          'lower_approx':float(lower),'upper_approx':float(upper),
                                          'width_upper':str(upper-lower)}
                reached = upper-lower <= F(1,500)
            else:
                reached = False
                attempt['reason'] = 'NONPOSITIVE_NORMALIZER_LOWER_BOUND'
            record['attempts'].append(attempt)
            if reached:
                record['status'] = 'TOLERANCE_MET'
                record['probability'] = attempt['probability']
            elif time.perf_counter()-started > 180:
                record['status'] = 'TIME_BUDGET_EXCEEDED_AFTER_MESH'
            record['seconds'] = time.perf_counter()-started
            record['attempt_count'] = len(record['attempts'])
            save(path,record)
            if record['status'] != 'RUNNING':
                break
        if record['status'] == 'RUNNING':
            record['status'] = 'MAXIMUM_MESH_WITHOUT_TARGET'
    except Exception as error:
        record['status'] = 'ERROR'
        record['error'] = repr(error)
    record['seconds'] = time.perf_counter()-started
    record['attempt_count'] = len(record['attempts'])
    save(path,record)
    return {k:record[k] for k in ['id','status','seconds','attempt_count']}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--limit',type=int)
    args=parser.parse_args()
    rows=json.loads((HERE/'test_requests.json').read_text())
    if args.limit is not None:rows=rows[:args.limit]
    start=time.perf_counter()
    run={'started_at':datetime.now(timezone.utc).isoformat(), 'workers':args.workers,
         'python':platform.python_version(),'python_flint':flint.__version__,
         'machine':platform.machine(),'platform':platform.system(),
         'requests':len(rows),'cases':[]}
    runpath=HERE/('certification_run.json' if args.limit is None else 'partial_run.json')
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(certify,row) for row in rows]
        for completed in as_completed(futures):
            result=completed.result();run['cases'].append(result)
            save(runpath,run)
            count=len(run['cases'])
            if count%20==0 or result['status']!='TOLERANCE_MET' or args.limit:
                print(json.dumps({'completed':count,'of':len(rows),**result}),flush=True)
    run['wall_seconds']=time.perf_counter()-start
    run['completed_at']=datetime.now(timezone.utc).isoformat()
    run['status']='ALL_REQUESTS_HANDLED'
    save(runpath,run)
    print(json.dumps({'status':run['status'],'seconds':run['wall_seconds'],'cases':len(run['cases'])}),flush=True)


if __name__=='__main__':main()
