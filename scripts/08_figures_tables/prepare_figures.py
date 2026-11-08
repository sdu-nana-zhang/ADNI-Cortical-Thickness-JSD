import sys,os,shutil
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from project import REPO,stage

out=stage('figures')/'data'
out.mkdir(exist_ok=True)
sources={
    'adni_scores.csv':stage('score')/'results/scored_MRI.csv',
    'adni_trajectory_effects.csv':stage('core')/'results/trajectory_results.csv',
    'adni_trajectory_included.csv':stage('core')/'audit/trajectory_included.csv',
    'adni_cognition_effects.csv':stage('core')/'results/cognition_results.csv',
    'adni_cognition_included.csv':stage('core')/'audit/cognition_included.csv',
    'adni_edge_effects.csv':stage('edges')/'results/all_edge_effects.csv',
    'roi_dictionary.csv':REPO/'config/ROI_dictionary.csv',
}
for name,path in sources.items():shutil.copy2(path,out/name)
edges=pd.read_csv(out/'adni_edge_effects.csv',float_precision='round_trip')
primary=edges[edges.MODEL.eq('L0')]
sig=primary[primary.CONTRAST.eq('AD-CN')&primary.Q_BH1891.lt(.05)].copy()
for group,col in [('CN','CN_SLOPE'),('AD','AD_SLOPE')]:
    values=primary[primary.CONTRAST.eq(group)][['EDGE_ID','ESTIMATE']].rename(columns={'ESTIMATE':col})
    sig=sig.merge(values,on='EDGE_ID',validate='one_to_one')
sig['QUADRANT']=np.select([
    sig.CN_SLOPE.lt(0)&sig.AD_SLOPE.lt(0),sig.CN_SLOPE.gt(0)&sig.AD_SLOPE.gt(0),
    sig.CN_SLOPE.lt(0)&sig.AD_SLOPE.gt(0),sig.CN_SLOPE.gt(0)&sig.AD_SLOPE.lt(0)],
    ['Both decreasing','Both increasing','CN decreasing; AD increasing','CN increasing; AD decreasing'],default='Zero slope')
sig.to_csv(out/'significant277_edges.csv',index=False,float_format='%.17g')
sig[sig.HOMOTOPIC].to_csv(out/'homotopic4_edges.csv',index=False,float_format='%.17g')
roi=pd.read_csv(out/'roi_dictionary.csv')
degree=pd.concat([sig.ROI1,sig.ROI2]).value_counts()
roi['DEGREE']=roi.ROI.map(degree).fillna(0).astype(int)
roi.to_csv(out/'regional_involvement.csv',index=False)
summary=stage('edges')/'results'
sig.to_csv(summary/'primary_discoveries.csv',index=False,float_format='%.17g')
sig[sig.HOMOTOPIC].to_csv(summary/'homotopic_discoveries.csv',index=False,float_format='%.17g')
roi.to_csv(summary/'regional_involvement.csv',index=False)
sig.QUADRANT.value_counts().rename_axis('PATTERN').reset_index(name='N_CONNECTIONS').to_csv(summary/'pattern_counts.csv',index=False)
trajectory=pd.read_csv(out/'adni_trajectory_effects.csv')
trajectory['Q_BH3']=np.nan
for model,g in trajectory[trajectory.TYPE.eq('SLOPE_DIFFERENCE')].groupby('MODEL'):
    p=g.P.to_numpy();ix=np.argsort(p)
    q=np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    original=np.empty_like(q);original[ix]=np.minimum(q,1)
    trajectory.loc[g.index,'Q_BH3']=original
trajectory.to_csv(stage('core')/'results/trajectory_with_q.csv',index=False,float_format='%.17g')
anatomy=os.environ.get('ADNI_ANATOMY_DIR')
if anatomy:
    target=out/'anatomy';target.mkdir(exist_ok=True)
    for hemi in ['lh','rh']:
        for sub,ext in [('surf','pial'),('label','aparc.DKTatlas.annot')]:
            shutil.copy2(Path(anatomy)/sub/f'{hemi}.{ext}',target/f'{hemi}.{ext}')
print('Prepared figure inputs from saved model outputs')
