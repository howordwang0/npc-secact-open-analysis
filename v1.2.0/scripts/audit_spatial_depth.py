"""Section-level depth/detected-gene adjusted spatial sensitivity audit."""
import argparse,json,math
from pathlib import Path
import numpy as np,pandas as pd

def rank(x): return pd.Series(np.asarray(x,float)).rank(method='average').to_numpy()
def corr(a,b): return float(np.corrcoef(np.asarray(a,float),np.asarray(b,float))[0,1])
def sign_p(pos,n):
    k=min(pos,n-pos); return min(1.0,2.0*sum(math.comb(n,j) for j in range(k+1))/(2**n))
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--spatial-dir',required=True,help='Directory with GSE206245_spot_scores.tsv.gz and ST*_selected_group_activity.tsv.gz')
    ap.add_argument('--out-dir',default='.'); a=ap.parse_args(); root=Path(a.spatial_dir); out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    spots=pd.read_csv(root/'GSE206245_spot_scores.tsv.gz',sep='\t'); rows=[]
    for f in sorted(root.glob('ST*_selected_group_activity.tsv.gz')):
        sample=f.name.split('_',1)[0]; meta=spots.loc[(spots['sample']==sample)&spots['passes_fixed_spot_qc']].set_index('barcode')
        act=pd.read_csv(f,sep='\t',index_col=0); act.index=act.index.astype(str).str.upper(); shared=act.columns.intersection(meta.index)
        cov=np.column_stack([np.log1p(meta.loc[shared,'nCount'].to_numpy(float)),meta.loc[shared,'nFeature'].to_numpy(float)])
        design=np.column_stack([np.ones(len(shared)),np.column_stack([rank(cov[:,j]) for j in range(2)])]); marker=meta.loc[shared,'myeloid_context'].to_numpy(float)
        for prot in ['WNT7A','VSTM2L']:
            activity=act.loc[prot,shared].to_numpy(float)
            rx=rank(activity)-design@np.linalg.lstsq(design,rank(activity),rcond=None)[0]
            ry=rank(marker)-design@np.linalg.lstsq(design,rank(marker),rcond=None)[0]
            rows.append({'sample':sample,'protein':prot,'n_spots':len(shared),'raw_rho':corr(rank(activity),rank(marker)),'depth_adjusted_partial_rho':corr(rx,ry)})
    df=pd.DataFrame(rows); df.to_csv(out/'candidate_spatial_depth_adjustment.tsv',sep='\t',index=False)
    summary={}
    for prot,g in df.groupby('protein'):
        r=g.depth_adjusted_partial_rho.to_numpy(float); loo=[float(np.median(np.delete(r,i))) for i in range(len(r))]
        summary[prot]={'n_sections':len(r),'adjusted_positive_sections':int((r>0).sum()),'adjusted_median_rho':float(np.median(r)),'adjusted_exact_sign_two_sided_p':sign_p(int((r>0).sum()),len(r)),'leave_one_section_out_median_min':min(loo),'leave_one_section_out_median_max':max(loo)}
    (out/'candidate_spatial_depth_adjustment.json').write_text(json.dumps(summary,indent=2),encoding='utf-8'); print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
