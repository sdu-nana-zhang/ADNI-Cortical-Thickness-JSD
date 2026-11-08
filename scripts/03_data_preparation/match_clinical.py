import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from project import DATA,stage
import pandas as pd

features = pd.read_csv(stage('jsd')/'JSD_features.csv',float_precision='round_trip').drop(columns=['ROW','DX'],errors='ignore')
clinical = pd.read_csv(DATA/'clinical_mri.csv',float_precision='round_trip',low_memory=False)
keys = ['PTID','MRI_IMAGE_ID']
if features.duplicated(keys).any() or clinical.duplicated(keys).any():
    raise ValueError('MRI keys must be unique in both input tables')
counts = features.groupby('PTID').size()
features = features[features.PTID.isin(counts[counts.ge(2)].index)].copy()
features['MRI_COUNT_VALID'] = features.PTID.map(counts)
features['MULTI_VISIT_ELIGIBLE'] = True
features = features.copy()
clinical_columns = [c for c in clinical if c not in features or c in keys]
merged = features.merge(clinical[clinical_columns],on=keys,how='left',validate='one_to_one',indicator=True,sort=False)
if not merged['_merge'].eq('both').all():
    raise ValueError('Every retained MRI requires an exact clinical-table key match')
merged.drop(columns='_merge').to_csv(stage('data')/'clinical_jsd.csv',index=False,encoding='utf-8-sig',float_format='%.17g')
print('Matched',merged.PTID.nunique(),'participants and',len(merged),'MRI records')
