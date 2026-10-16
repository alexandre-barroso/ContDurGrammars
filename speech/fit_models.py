"""Fit the prespecified small timing models; select tolerance on validation only."""
from pathlib import Path
from fractions import Fraction as F
import json,math
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr

HERE=Path(__file__).resolve().parent

def feature(syllable,position,phones):
    stress=syllable['stress']
    return np.array([1,stress,position,stress*position]+[syllable['phones'].count(p) for p in phones],dtype=float)

def unary_probability(mu,sigma,threshold):
    return float(ndtr((mu-threshold)/sigma)/ndtr(mu/sigma))

def proportion_probability(r,mean_residual,sigma):
    mu=r+mean_residual
    return float((ndtr((mu-.5)/sigma)-ndtr((mu-1)/sigma))/(ndtr((1-mu)/sigma)-ndtr(-mu/sigma)))

def joint_probability(mu,sd,total,tau):
    # Independent check for the quadratic special case only. The reported
    # certificate is produced separately by the general integrator.
    precision=np.diag(1/np.square(sd))+np.ones((2,2))/tau**2
    covariance=np.linalg.inv(precision)
    center=covariance@(np.array(mu)/np.square(sd)+total/tau**2)
    s=math.sqrt(covariance[0,0])
    conditional_sd=math.sqrt(covariance[1,1]-covariance[0,1]**2/covariance[0,0])
    def integrand(z):
        return math.exp(-.5*z*z)/math.sqrt(2*math.pi)*ndtr((center[1]+covariance[0,1]/s*z)/conditional_sd)
    def mass(lo):
        return quad(integrand,(lo-center[0])/s,np.inf,epsabs=1e-12,epsrel=1e-10,limit=150)[0]
    den=mass(0);num=mass(total/2)
    if not den>0:raise ArithmeticError('empty normalizing control')
    return max(0.,min(1.,num/den))

def model_inputs(row,fit):
    means=[max(5.,float(feature(s,i,fit['phones'])@fit['coefficients'])) for i,s in enumerate(row['syllables'])]
    sds=[fit['residual_sd_ms'][f'{s["stress"]}/{i}'] for i,s in enumerate(row['syllables'])]
    return means,sds

def main():
    rows=json.loads((HERE/'prepared_words.json').read_text());fits={};validation={}
    for speaker in ['slt','bdl']:
        train=[r for r in rows if r['speaker']==speaker and r['split']=='train']
        phones=sorted({p for r in train for s in r['syllables'] for p in s['phones']})
        items=[(r,i,s) for r in train for i,s in enumerate(r['syllables'])]
        X=np.array([feature(s,i,phones) for _,i,s in items]);Y=np.array([float(F(s['duration_ms'])) for _,i,s in items])
        penalty=10*np.eye(X.shape[1]);penalty[0,0]=0
        beta=np.linalg.solve(X.T@X+penalty,X.T@Y)
        residual=Y-np.maximum(5,X@beta)
        variances={}
        for stress in sorted({s['stress'] for _,_,s in items}):
            for pos in [0,1]:
                e=np.array([residual[j] for j,(_,i,s) in enumerate(items) if s['stress']==stress and i==pos])
                variances[f'{stress}/{pos}']=max(10.,float(np.sqrt(np.mean(e**2)))) if len(e) else max(10.,float(np.sqrt(np.mean(residual**2))))
        fit={'phones':phones,'coefficients':beta.tolist(),'residual_sd_ms':variances,'ridge_penalty':10,'training_words':len(train)}
        residual_ratios=[]
        for row in train:
            mu,sd=model_inputs(row,fit)
            residual_ratios.append(float(F(row['syllables'][0]['duration_ms'])/F(row['total_ms']))-mu[0]/sum(mu))
        fit['ratio_residual_mean']=float(np.mean(residual_ratios));fit['ratio_residual_sd']=max(.01,float(np.std(residual_ratios)))
        fits[speaker]=fit;validation[speaker]={}
        val=[r for r in rows if r['speaker']==speaker and r['split']=='validation']
        for tau in [5,10,20,40]:
            loss=[]
            for row in val:
                mu,sd=model_inputs(row,fit);total=float(F(row['total_ms']))
                y=float(F(row['syllables'][0]['duration_ms'])>F(row['total_ms'])/2)
                p=joint_probability(mu,sd,total,tau);loss.append((p-y)**2)
            validation[speaker][str(tau)]=float(np.mean(loss))
    pooled={str(t):float(np.mean([validation[s][str(t)] for s in fits])) for t in [5,10,20,40]}
    tau=int(min(pooled,key=pooled.get))
    result={'fits':fits,'validation_brier':validation,'validation_equal_voice_mean':pooled,'selected_group_tolerance_ms':tau}
    (HERE/'fitted_models.json').write_text(json.dumps(result,indent=2)+'\n')
    requests=[]
    for row in rows:
        if row['split']!='test':continue
        fit=fits[row['speaker']];mu,sd=model_inputs(row,fit);total=float(F(row['total_ms']))
        unknown=set(p for s in row['syllables'] for p in s['phones'])-set(fit['phones'])
        assert not unknown,(row['id'],unknown)
        requests.append({**row,'target_ms_exact':[str(F(m)) for m in mu],
                         'unary_weight_exact':[str(F(1/(2*s*s))) for s in sd],
                         'group_weight_exact':str(F(1,2*tau*tau)),
                         'independent_probability':unary_probability(mu[0],sd[0],total/2),
                         'proportion_probability':proportion_probability(mu[0]/sum(mu),fit['ratio_residual_mean'],fit['ratio_residual_sd']),
                         'quadratic_probability_control':joint_probability(mu,sd,total,tau)})
    (HERE/'test_requests.json').write_text(json.dumps(requests,indent=2)+'\n')
    print('Validation only',json.dumps(validation,indent=2),flush=True)
    print('Selected tolerance',tau,'ms; requests',len(requests),flush=True)

if __name__=='__main__':main()
