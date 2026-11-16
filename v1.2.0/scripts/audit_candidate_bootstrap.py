"""Patient-stratified Hedges-g interval and rank-stability audit for WNT7A/VSTM2L."""
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind

def hedges(a,b):
    na,nb=a.shape[-1],b.shape[-1]
    v=((na-1)*np.var(a,axis=-1,ddof=1)+(nb-1)*np.var(b,axis=-1,ddof=1))/(na+nb-2)
    return (1-3/(4*(na+nb-2)-1))*(np.mean(a,axis=-1)-np.mean(b,axis=-1))/np.sqrt(v)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--zscore',required=True,help='GSE153125 SecAct z-score TSV/TSV.GZ')
    ap.add_argument('--manifest',required=True,help='GSE153125 patient manifest TSV')
    ap.add_argument('--out',default='GSE153125_candidate_stability.json')
    ap.add_argument('--effect-bootstrap',type=int,default=20000)
    ap.add_argument('--rank-bootstrap',type=int,default=5000)
    a=ap.parse_args()
    m=pd.read_csv(a.manifest,sep='\t'); z=pd.read_csv(a.zscore,sep='\t',index_col=0); z.index=z.index.astype(str)
    cr=m.loc[m.response.eq('CR'),'sample'].tolist(); pr=m.loc[m.response.eq('PR'),'sample'].tolist()
    x=z[cr].to_numpy(float); y=z[pr].to_numpy(float); features=z.index.tolist()
    out={}
    for protein,seed in [('WNT7A',20481115),('VSTM2L',20481116)]:
        i=features.index(protein); rng=np.random.default_rng(seed)
        ix=rng.integers(0,x.shape[1],size=(a.effect_bootstrap,x.shape[1])); iy=rng.integers(0,y.shape[1],size=(a.effect_bootstrap,y.shape[1]))
        gb=hedges(x[i,ix],y[i,iy])
        ranks=[]; loo_sign=[]
        for lo in range(0,a.rank_bootstrap,100):
            n=min(100,a.rank_bootstrap-lo)
            bx=x[:,rng.integers(0,x.shape[1],size=(n,x.shape[1]))].transpose(1,0,2)
            by=y[:,rng.integers(0,y.shape[1],size=(n,y.shape[1]))].transpose(1,0,2)
            bg=hedges(bx,by); candidate=bg[:,i]
            ranks.extend((1+(np.abs(bg)>np.abs(candidate[:,None])).sum(axis=1)).tolist())
        for j in range(x.shape[1]): loo_sign.append(float(np.sign(hedges(np.delete(x,j,axis=1),y)[i])))
        for j in range(y.shape[1]): loo_sign.append(float(np.sign(hedges(x,np.delete(y,j,axis=1))[i])))
        out[protein]={'hedges_g':float(hedges(x[i],y[i])),'bootstrap_95ci_percentile':[float(v) for v in np.quantile(gb,[.025,.975])], 'leave_one_patient_out_same_direction':int(sum(v==np.sign(hedges(x[i],y[i])) for v in loo_sign)), 'bootstrap_abs_g_rank_median':float(np.median(ranks)), 'bootstrap_top67_fraction':float(np.mean(np.asarray(ranks)<=67)), 'effect_bootstrap_replicates':a.effect_bootstrap, 'rank_bootstrap_replicates':a.rank_bootstrap}
    Path(a.out).write_text(json.dumps(out,indent=2),encoding='utf-8'); print(json.dumps(out,indent=2))
if __name__=='__main__':main()
