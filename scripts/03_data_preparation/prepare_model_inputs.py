import sys
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from project import REPO,stage
ROOT=stage('score')
for name in ['inputs','protocol','models','results','audit','work']:(ROOT/name).mkdir(exist_ok=True)
master=Path(sys.argv[1]) if len(sys.argv)>1 else stage('data')/'clinical_jsd.csv'
def save(df,p):
    p.parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(p,index=False,encoding='utf-8-sig',float_format='%.17g',date_format='%Y-%m-%d')
dic=pd.read_csv(REPO/'config/feature_dictionary.csv')
source_cols=['PTID','MRI_IMAGE_ID','MRI_DATE','DX.bl','AGE_AT_MRI_EST','PTGENDER','PTEDUCAT','CORTEX_CM3','DKT62_AREA_WEIGHTED_THICKNESS_MM','ADAS13','MMSE','CDRSB','MOCA','EXAMDATE','MRI_YEARS_FROM_BASELINE','MRI_COUNT_VALID','MULTI_VISIT_ELIGIBLE','RID','VISCODE','SITE','COLPROT','MRI_Visit','DX']
mean_cols=dic.loc[dic.BLOCK.eq('ROI_MEAN'),'SOURCE_COLUMN'].tolist()
sd_cols=dic.loc[dic.BLOCK.eq('ROI_SD'),'SOURCE_COLUMN'].tolist()
jsd_cols=[f'JSD_{i:04d}' for i in range(1,1892)]
use=[]
for c in source_cols+mean_cols+sd_cols+jsd_cols:
 if c not in use: use.append(c)
print('reading columns',len(use),flush=True)
df=pd.read_csv(master,usecols=use,float_precision='round_trip').copy()
print('read rows',len(df),'participants',df.PTID.nunique(),flush=True)
df['MRI_DATE']=pd.to_datetime(df['MRI_DATE'])
df=df.sort_values(['PTID','MRI_DATE','MRI_IMAGE_ID'],kind='mergesort').reset_index(drop=True)
df['DX_STD']=df['DX.bl'].replace({'LMCI':'MCI','EMCI':'MCI'})
df['MALE']=df['PTGENDER'].eq('Male').astype(int)
df['AGE']=pd.to_numeric(df['AGE_AT_MRI_EST'],errors='coerce')
feature_names=dic.FEATURE.tolist()
feature_sources=[]
for _,r in dic.iterrows():
 if r.BLOCK=='JSD': feature_sources.append(f'JSD_{int(r.FEATURE.split("_")[-1]):04d}')
 else: feature_sources.append(r.SOURCE_COLUMN)
missing=[c for c in feature_sources if c not in df.columns]
assert not missing, missing[:10]
assert df[feature_sources].notna().all().all(), 'missing feature values'
assert np.isfinite(df[feature_sources].to_numpy(float)).all(), 'nonfinite feature values'
assert (df[jsd_cols].to_numpy(float)>=0).all() and (df[jsd_cols].to_numpy(float)<=1).all()
first=df.drop_duplicates('PTID',keep='first').copy()
first['CLASS']=first['DX_STD']
first=first[first.CLASS.isin(['CN','AD'])].copy()
first['Y']=first.CLASS.eq('AD').astype(int)
need=['AGE','MALE','PTEDUCAT','CORTEX_CM3','DKT62_AREA_WEIGHTED_THICKNESS_MM']+feature_sources
first=first.loc[first[need].notna().all(axis=1)].copy()
first=first.sort_values('PTID').reset_index(drop=True)
model=pd.DataFrame({'PTID':first.PTID.astype(str),'MRI_IMAGE_ID':first.MRI_IMAGE_ID.astype(str),'MRI_DATE':first.MRI_DATE,'CLASS':first.CLASS,'Y':first.Y.astype(int),'AGE':first.AGE.astype(float),'MALE':first.MALE.astype(int),'EDUCATION':first.PTEDUCAT.astype(float),'CORTEX_CM3':first.CORTEX_CM3.astype(float),'DKT62_AREA_WEIGHTED_THICKNESS_MM':first.DKT62_AREA_WEIGHTED_THICKNESS_MM.astype(float)})
model=pd.concat([model,pd.DataFrame({f:first[s].to_numpy(float) for f,s in zip(feature_names,feature_sources)})],axis=1)
save(model,ROOT/'inputs/classifier_input.csv')
allout=pd.DataFrame({'PTID':df.PTID.astype(str),'MRI_IMAGE_ID':df.MRI_IMAGE_ID.astype(str),'MRI_DATE':df.MRI_DATE,'DX_STD':df.DX_STD,'AGE':df.AGE,'MALE':df.MALE,'EDUCATION':df.PTEDUCAT.astype(float),'CORTEX_CM3':df.CORTEX_CM3.astype(float),'DKT62_AREA_WEIGHTED_THICKNESS_MM':df.DKT62_AREA_WEIGHTED_THICKNESS_MM.astype(float),'ADAS13':df.ADAS13,'MMSE':df.MMSE,'CDRSB':df.CDRSB,'MOCA':df.MOCA,'EXAMDATE':df.EXAMDATE,'MRI_YEARS_FROM_BASELINE':df.MRI_YEARS_FROM_BASELINE,'MRI_COUNT_VALID':df.MRI_COUNT_VALID,'MULTI_VISIT_ELIGIBLE':df.MULTI_VISIT_ELIGIBLE,'RID':df.RID,'VISCODE':df.VISCODE,'SITE':df.SITE,'COLPROT':df.COLPROT,'MRI_Visit':df.MRI_Visit,'DX':df.DX,'DX.bl':df['DX.bl']})
allout=pd.concat([allout,pd.DataFrame({f:df[s].to_numpy(float) for f,s in zip(feature_names,feature_sources)})],axis=1)
save(allout,ROOT/'inputs/all_scans.csv')
save(dic,ROOT/'inputs/feature_dictionary.csv')
cfg=json.loads((REPO/'config/classification.json').read_text(encoding='utf-8'))
cfg.update({'root':str(ROOT),'n':len(model),'class_counts':{'CN':int((model.Y==0).sum()),'AD':int((model.Y==1).sum())},'frozen_utc':datetime.now(timezone.utc).isoformat(),'point_estimate':f'average five held-out probabilities per patient, then calculate metrics on {len(model)} unique patients'})
(ROOT/'protocol/frozen_protocol.json').write_text(json.dumps(cfg,ensure_ascii=False,indent=2),encoding='utf-8')
def rvalue(x):
 if isinstance(x,bool): return 'TRUE' if x else 'FALSE'
 if isinstance(x,str): return json.dumps(x,ensure_ascii=False)
 if isinstance(x,(int,float)): return str(x)
 if isinstance(x,list): return 'c('+','.join(rvalue(v) for v in x)+')'
 if isinstance(x,dict): return 'list('+','.join(json.dumps(str(k))+'='+rvalue(v) for k,v in x.items())+')'
 raise TypeError(type(x))
(ROOT/'protocol/frozen_protocol.R').write_text(rvalue(cfg),encoding='utf-8')
def folds(y,k,seed):
 rng=np.random.default_rng(seed); v=np.empty(len(y),int)
 for c in [0,1]:
  ix=np.flatnonzero(y==c); rng.shuffle(ix); v[ix]=np.arange(len(ix))%k
 return v+1
outer=[]; inner=[]; y=model.Y.to_numpy(int)
for rep,seed in enumerate(cfg['repeat_seeds'],1):
 of=folds(y,5,seed)
 for i in range(len(model)): outer.append({'PTID':model.PTID[i],'ROW1':i+1,'REPEAT':rep,'OUTER_FOLD':int(of[i]),'Y':int(y[i])})
 for fo in range(1,6):
  tr=np.flatnonzero(of!=fo); iv=folds(y[tr],4,seed+100+fo-1)
  for loc,ix in enumerate(tr): inner.append({'PTID':model.PTID[ix],'ROW1':int(ix+1),'REPEAT':rep,'OUTER_FOLD':fo,'INNER_FOLD':int(iv[loc])})
save(pd.DataFrame(outer),ROOT/'protocol/outer_patient_folds.csv'); save(pd.DataFrame(inner),ROOT/'protocol/inner_patient_folds.csv')
summary={'rows':len(df),'participants':int(df.PTID.nunique()),'first_cn_ad_complete':len(model),'class_counts':model.CLASS.value_counts().to_dict(),'feature_count':len(feature_names),'jsd_count':len(jsd_cols)}
(ROOT/'audit/input_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False),flush=True)
