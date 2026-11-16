from pathlib import Path
import argparse, json, zipfile
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu

BASE=Path(r'C:\Users\Lenovo\Documents\Codex\2026-11-15\new-chat')
OUT=BASE/'outputs'/'Supplementary_Audit_Files'
OUT.mkdir(parents=True,exist_ok=True)

def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); a=p[o]*len(p)/np.arange(1,len(p)+1); a=np.minimum.accumulate(a[::-1])[::-1]; q=np.empty_like(a); q[o]=np.minimum(a,1.0); return q

def screen(zpath, mpath):
    manifest=pd.read_csv(mpath,sep='\t')
    z=pd.read_csv(zpath,sep='\t',index_col=0); z.index=z.index.astype(str)
    cr=manifest.loc[manifest.response.eq('CR'),'sample'].tolist(); pr=manifest.loc[manifest.response.eq('PR'),'sample'].tolist()
    x=z[cr].to_numpy(float); y=z[pr].to_numpy(float); n1,n2=x.shape[1],y.shape[1]
    pooled=((n1-1)*x.var(axis=1,ddof=1)+(n2-1)*y.var(axis=1,ddof=1))/(n1+n2-2)
    d=(x.mean(axis=1)-y.mean(axis=1))/np.sqrt(pooled)
    j=1-3/(4*(n1+n2-2)-1); g=j*d
    wp=np.array([ttest_ind(x[i],y[i],equal_var=False).pvalue for i in range(len(z))])
    mp=np.array([mannwhitneyu(x[i],y[i],alternative='two-sided',method='exact').pvalue for i in range(len(z))])
    res=pd.DataFrame({'protein':z.index,'n_CR':n1,'mean_CR':x.mean(axis=1),'sd_CR':x.std(axis=1,ddof=1),'n_PR':n2,'mean_PR':y.mean(axis=1),'sd_PR':y.std(axis=1,ddof=1),'cohens_d_CR_minus_PR':d,'hedges_g_CR_minus_PR':g,'welch_p':wp,'welch_BH_q':bh(wp),'mann_whitney_exact_p':mp,'mann_whitney_BH_q':bh(mp)})
    return res

def main():
    ap=argparse.ArgumentParser(description='Reproduce the 1,170-program response screen from the public GSE153125 activity matrix and patient manifest.')
    ap.add_argument('--zscore',required=True,help='Tab-separated GSE153125 SecAct z-score matrix (proteins in rows, samples in columns).')
    ap.add_argument('--manifest',required=True,help='Patient manifest TSV with sample and response columns.')
    ap.add_argument('--out',default='GSE153125_all_program_results.tsv')
    a=ap.parse_args(); res=screen(a.zscore,a.manifest); res.to_csv(a.out,sep='\t',index=False,float_format='%.10g')
    s={'n_programs':int(len(res)),'welch_BH_q_lt_0.05':int((res.welch_BH_q<.05).sum()),'exact_MW_BH_q_lt_0.05':int((res.mann_whitney_BH_q<.05).sum()),'welch_only_BH_q_lt_0.05':int(((res.welch_BH_q<.05)&(res.mann_whitney_BH_q>=.05)).sum())}
    print(json.dumps(s,indent=2))
if __name__=='__main__':main()
