from pathlib import Path
import numpy as np
import shutil
import pandas as pd
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from project import REPO,stage
root=stage('score')
edge=stage('edges')
(edge/"inputs").mkdir(parents=True,exist_ok=True);(edge/"models").mkdir(exist_ok=True);(edge/"work").mkdir(exist_ok=True);(edge/"results").mkdir(exist_ok=True);(edge/"audit").mkdir(exist_ok=True);(edge/"scripts").mkdir(exist_ok=True)
p=root/"results"/"scored_MRI.csv"
cols=['PTID','MRI_IMAGE_ID','MRI_DATE','DX_STD','AGE','MALE','EDUCATION','CORTEX_CM3','DKT62_AREA_WEIGHTED_THICKNESS_MM']+[f'JSD_{i:04d}' for i in range(1,1892)]
d=pd.read_csv(p,usecols=cols,encoding='utf-8-sig',low_memory=False)
d=d.sort_values(['PTID','MRI_DATE','MRI_IMAGE_ID'],kind='mergesort').reset_index(drop=True)
first=d.drop_duplicates('PTID',keep='first')[['PTID','MRI_DATE','DX_STD','AGE','MALE','EDUCATION']].rename(columns={'MRI_DATE':'BASELINE_MRI_DATE','DX_STD':'BASELINE_DX','AGE':'BASELINE_AGE','MALE':'BASELINE_MALE','EDUCATION':'BASELINE_EDUCATION'})
d=d.merge(first,on='PTID',how='left',validate='many_to_one',sort=False)
ids=first.loc[first.BASELINE_DX.isin(['CN','AD']) & first[['BASELINE_AGE','BASELINE_MALE','BASELINE_EDUCATION']].notna().all(axis=1),'PTID']
q=d[d.PTID.isin(set(ids))].copy().sort_values(['PTID','MRI_DATE','MRI_IMAGE_ID'],kind='mergesort').reset_index(drop=True)
q['TIME']=(pd.to_datetime(q.MRI_DATE)-pd.to_datetime(q.BASELINE_MRI_DATE)).dt.days/365.25
q['AGE_C']=q.BASELINE_AGE-75;q['MALE']=q.BASELINE_MALE;q['EDUC_C']=q.BASELINE_EDUCATION-16;q['CORTEX100']=q.CORTEX_CM3/100;q['THICKNESS']=q.DKT62_AREA_WEIGHTED_THICKNESS_MM
q['PTID']=q.PTID.astype(str);q['BASELINE_DX']=q.BASELINE_DX.astype(str)
q['CORTEX_BETWEEN']=q.groupby('PTID').CORTEX100.transform('mean');q['CORTEX_WITHIN']=q.CORTEX100-q.CORTEX_BETWEEN
q['THICKNESS_BETWEEN']=q.groupby('PTID').THICKNESS.transform('mean');q['THICKNESS_WITHIN']=q.THICKNESS-q.THICKNESS_BETWEEN
assert not q.duplicated(['PTID','MRI_IMAGE_ID']).any() and set(q.BASELINE_DX)=={'CN','AD'}
keep=['PTID','MRI_IMAGE_ID','MRI_DATE','BASELINE_DX','TIME','AGE_C','MALE','EDUC_C','CORTEX100','THICKNESS','CORTEX_BETWEEN','CORTEX_WITHIN','THICKNESS_BETWEEN','THICKNESS_WITHIN']
q[keep].to_csv(edge/'inputs'/'covariates.csv',index=False,encoding='utf-8-sig',float_format='%.17g')
Y=q[[f'JSD_{i:04d}' for i in range(1,1892)]].to_numpy(dtype=np.float64)
(edge/'inputs'/'JSD_columnmajor_float64.bin').write_bytes(Y.tobytes(order='F'))

shutil.copy2(REPO/'config/edge_dictionary.csv',edge/'inputs/edge_dictionary.csv')
shutil.copy2(REPO/'config/ROI_dictionary.csv',edge/'inputs/ROI_dictionary.csv')
(edge/'protocol').mkdir(exist_ok=True)
cfg=json.loads((REPO/'config/edges.json').read_text(encoding='utf-8'))
cfg.update(n_patients=int(q.PTID.nunique()),n_MRI=len(q))
def rv(x):
    if isinstance(x,str):return json.dumps(x)
    if isinstance(x,list):return 'list('+','.join(rv(v) for v in x)+')'
    if isinstance(x,dict):return 'list('+','.join(json.dumps(k)+'='+rv(v) for k,v in x.items())+')'
    return str(x)
(edge/'protocol/frozen_config.R').write_text(rv(cfg),encoding='utf-8')
(edge/'protocol/frozen_config.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
print({'n_patients':int(q.PTID.nunique()),'n_mri':len(q)})
