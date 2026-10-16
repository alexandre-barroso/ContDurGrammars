"""Check split isolation, input independence, model identities and probability enclosures."""
from pathlib import Path
from fractions import Fraction as F
import copy
import json
from fit_models import model_inputs
from certify_predictions import model_tree, describe

HERE=Path(__file__).resolve().parent


def main():
    raw=json.loads((HERE/'prepared_words.json').read_text())
    req=json.loads((HERE/'test_requests.json').read_text())
    fit=json.loads((HERE/'fitted_models.json').read_text())
    assert len({r['id'] for r in raw})==len(raw)
    prompt_split={}
    assignments=json.loads((HERE/'split_assignments.json').read_text())
    for row in raw:
        u=row['utterance'];expected=assignments[u]
        assert row['split']==expected
        assert prompt_split.setdefault(u,expected)==expected
        assert len(row['syllables'])==2
        assert F(row['total_ms'])==sum(F(s['duration_ms']) for s in row['syllables'])
        assert all(F(s['duration_ms'])>0 for s in row['syllables'])
    assert {r['id'] for r in req}=={r['id'] for r in raw if r['split']=='test'}
    assert len(req)==836
    best=int(min(fit['validation_equal_voice_mean'],key=fit['validation_equal_voice_mean'].get))
    assert best==fit['selected_group_tolerance_ms']==40
    for row in req:
        fm=fit['fits'][row['speaker']]
        training=[r for r in raw if r['split']=='train' and r['speaker']==row['speaker']]
        assert fm['training_words']==len(training)
        assert set(fm['phones'])=={p for r in training for s in r['syllables'] for p in s['phones']}
        perturbed=copy.deepcopy(row)
        perturbed['syllables'][0]['duration_ms']='1'
        perturbed['syllables'][1]['duration_ms']=str(F(row['total_ms'])-1)
        assert model_inputs(row,fm)==model_inputs(perturbed,fm),'Internal timing leaked into features'
        means,sds=model_inputs(row,fm)
        assert list(map(str,map(F,means)))==row['target_ms_exact']
        assert [str(F(1/(2*s*s))) for s in sds]==row['unary_weight_exact']
        certificate=json.loads((HERE/'certificates'/(row['id'].replace('/','__')+'.json')).read_text())
        assert certificate['tree']==describe(model_tree(row))
        assert certificate['event']=={'x1_lower_ms':str(F(row['total_ms'])/2),'x1_upper_ms':None}
        assert certificate['status']=='TOLERANCE_MET'
        lo,hi=map(F,[certificate['probability']['lower'],certificate['probability']['upper']])
        assert 0<=lo<=hi<=1 and hi-lo<=F(1,500)
        assert float(lo)<=row['quadratic_probability_control']<=float(hi)
    result={'status':'PASS','data_rows':len(raw),'test_requests':len(req),
            'unique_test_prompts':len({r['utterance'] for r in req}),
            'split_isolation':True,'internal_duration_not_used_as_feature':True,
            'all_exact_probability_widths_at_most_0_002':True,
            'all_quadratic_controls_inside_enclosures':True,
            'certificate_models_match_requests':True,
            'reference_midpoint_ties_scored_zero':sum(F(r['syllables'][0]['duration_ms'])==F(r['total_ms'])/2 for r in req)}
    (HERE/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
