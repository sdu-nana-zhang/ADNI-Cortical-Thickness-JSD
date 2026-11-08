import sys,json
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from project import REPO,stage

expected=json.loads((REPO/'config/expected_results.json').read_text())
checks=[]
def check(name,ok):
    checks.append({'CHECK':name,'PASS':bool(ok)})
    if not ok:raise AssertionError(name)
s=pd.read_csv(stage('score')/'results/scored_MRI.csv')
check('Scored cohort',s.PTID.nunique()==expected['n_participants'] and len(s)==expected['n_mri'])
b=s[s.JSD_BASELINE_ROW]
check('Baseline groups',b.DX_STD.value_counts().to_dict()==expected['baseline_groups'])
check('Final score scale',abs(b.JSD_SCORE_Z.mean())<1e-10 and abs(b.JSD_SCORE_Z.std(ddof=1)-1)<1e-10)
cv=pd.read_csv(stage('score')/'results/classifier_performance.csv').set_index('MODEL').loc['M3']
for key in ['AUC','CI_LOW','CI_HIGH']:check('Full classifier '+key,abs(cv[key]-expected['classifier'][key])<1e-12)
t=pd.read_csv(stage('core')/'results/trajectory_results.csv')
for g,v in expected['trajectory']['slopes'].items():
    r=t[t.MODEL.eq('L0')&t.CONTRAST.eq(g)].iloc[0]
    check('Trajectory '+g,abs(r.ESTIMATE-v)<1e-8 and r.N_PATIENTS==expected['trajectory']['n'] and r.N_MRI==expected['trajectory']['n_mri'])
c=pd.read_csv(stage('core')/'results/cognition_results.csv')
for outcome,z in expected['cognition'].items():
    for term,key in [('SCORE_BETWEEN','between'),('SCORE_WITHIN','within')]:
        r=c[c.OUTCOME.eq(outcome)&c.TERM.eq(term)].iloc[0]
        check(outcome+' '+term,abs(r.BETA-z[key])<1e-8 and r.N_PATIENTS==z['n'] and r.N_MRI==z['n_mri'])
e=pd.read_csv(stage('edges')/'results/all_edge_effects.csv')
for (m,contrast),d in e.groupby(['MODEL','CONTRAST']):check('Full edge family '+m+contrast,len(d)==expected['connections']['tested'] and d.EDGE_ID.is_unique)
primary=e[e.MODEL.eq('L0')&e.CONTRAST.eq('AD-CN')]
sig=primary[primary.Q_BH1891.lt(.05)]
check('Discovery count',len(sig)==expected['connections']['significant'])
check('Homotopic targets',sorted(sig[sig.HOMOTOPIC].EDGE_ID.tolist())==expected['connections']['homotopic'])
degree=pd.concat([sig.ROI1,sig.ROI2]).value_counts()
check('Endpoint count',degree.sum()==expected['connections']['endpoint_sum'])
pd.DataFrame(checks).to_csv(stage('validation')/'study_results.csv',index=False)
print('PASS:',len(checks),'study-output checks')
