"""Score every held-out request, retain all comparisons and uncertainty sources."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent


def summarize(rows):
    out={'words':len(rows),'prompts':len({r['utterance'] for r in rows}),
         'event_frequency':float(np.mean([r['event'] for r in rows]))}
    for name in ['grammar','independent','proportion']:
        p=np.array([r[name] for r in rows]);y=np.array([r['event'] for r in rows])
        out[name]={'brier':float(np.mean((p-y)**2)),
                   'log_loss':float(np.mean(-y*np.log(np.clip(p,1e-15,1-1e-15))-(1-y)*np.log(np.clip(1-p,1e-15,1-1e-15))))}
    lows=[];highs=[]
    for r in rows:
        l,u=F(r['lower_exact']),F(r['upper_exact'])
        if r['event']:l,u=1-u,1-l
        lows.append(l*l);highs.append(u*u)
    out['grammar']['brier_enclosure']=[str(sum(lows)/len(rows)),str(sum(highs)/len(rows))]
    out['grammar']['brier_enclosure_approx']=list(map(lambda s:float(F(s)),out['grammar']['brier_enclosure']))
    out['max_probability_width']=max(r['width'] for r in rows)
    out['mesh_counts_ms']=dict(Counter(r['mesh_ms'] for r in rows))
    out['seconds_median']=float(np.median([r['seconds'] for r in rows]))
    out['seconds_max']=max(r['seconds'] for r in rows)
    return out


def main():
    requests=json.loads((HERE/'test_requests.json').read_text());rows=[];failures=[];mismatches=[]
    for req in requests:
        cert=json.loads((HERE/'certificates'/(req['id'].replace('/','__')+'.json')).read_text())
        if cert['status']!='TOLERANCE_MET':
            failures.append({'id':req['id'],'status':cert['status']});continue
        interval=cert['probability'];lo,hi=map(F,[interval['lower'],interval['upper']])
        assert hi-lo<=F(1,500),req['id']
        control=req['quadratic_probability_control']
        if not float(lo)-1e-10<=control<=float(hi)+1e-10:mismatches.append(req['id'])
        rows.append({'id':req['id'],'speaker':req['speaker'],'utterance':req['utterance'],
                     'word':req['word'],'event':int(F(req['syllables'][0]['duration_ms'])>F(req['total_ms'])/2),
                     'grammar':float((lo+hi)/2),'independent':req['independent_probability'],
                     'proportion':req['proportion_probability'],
                     'quadratic_control':control,'lower_exact':str(lo),'upper_exact':str(hi),
                     'width':float(hi-lo),'mesh_ms':cert['attempts'][-1]['mesh_ms'],'seconds':cert['seconds']})
    result={'requested_words':len(requests),'successful_words':len(rows),'failures':failures,
            'quadratic_control_outside_enclosure':mismatches}
    assert not failures and not mismatches,(failures,mismatches)
    for voice in ['slt','bdl']:
        result[voice]=summarize([r for r in rows if r['speaker']==voice])
    result['pooled']=summarize(rows)
    # The test split has identical word counts and prompt membership in both
    # voices. Compute equal-voice means explicitly for every bootstrap sample.
    prompts=sorted({r['utterance'] for r in rows})
    byvoice={s:{u:[r for r in rows if r['speaker']==s and r['utterance']==u] for u in prompts} for s in ['slt','bdl']}
    rng=np.random.default_rng(20261016);bs={b:[] for b in ['independent','proportion']}
    # Cluster sums/counts avoid materializing resampled row lists.
    count=np.array([[len(byvoice[s][u]) for u in prompts] for s in byvoice])
    diffs={b:np.array([[sum((r['grammar']-r['event'])**2-(r[b]-r['event'])**2 for r in byvoice[s][u]) for u in prompts] for s in byvoice]) for b in bs}
    for _ in range(2000):
        idx=rng.integers(0,len(prompts),len(prompts))
        for b in bs:bs[b].append(float(np.mean(np.sum(diffs[b][:,idx],axis=1)/np.sum(count[:,idx],axis=1))))
    result['brier_differences']={b:{'grammar_minus_baseline':result['pooled']['grammar']['brier']-result['pooled'][b]['brier'],
                                  'paired_prompt_bootstrap_95_percentile_interval':np.quantile(bs[b],[.025,.975]).tolist()}
                                  for b in bs}
    result['bootstrap']={'replicates':2000,'seed':20261016,'unit':'shared prompt, all words and both voices retained','unique_test_prompts':len(prompts)}
    result['equal_voice_mean_brier']={b:float(np.mean([result[s][b]['brier'] for s in ['slt','bdl']])) for b in ['grammar','independent','proportion']}
    result['first_test_example']={'request':requests[0],'prediction':rows[0]}
    (HERE/'test_predictions.json').write_text(json.dumps(rows,indent=2)+'\n')
    (HERE/'evaluation.json').write_text(json.dumps(result,indent=2)+'\n')
    display={k:v for k,v in result.items() if k!='first_test_example'}
    for k in ['slt','bdl','pooled']:display[k]['grammar'].pop('brier_enclosure')
    print(json.dumps(display,indent=2))


if __name__=='__main__':main()
