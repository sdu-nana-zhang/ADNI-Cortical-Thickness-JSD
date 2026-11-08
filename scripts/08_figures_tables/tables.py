import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from project import stage

scans=pd.read_csv(stage('score')/'results/scored_MRI.csv',float_precision='round_trip')
base=scans[scans.JSD_BASELINE_ROW].copy()
n=len(base)
rows=[['Participants, n',str(n)],['MRI examinations, n',f'{len(scans):,}']]
for label,col in [('Age, mean (SD), y','AGE'),('Male sex, n (%)','MALE'),('Education, mean (SD), y','EDUCATION')]:
    x=base[col].dropna()
    value=f'{int(x.sum())} ({100*x.mean():.1f}%)' if col=='MALE' else f'{x.mean():.1f} ({x.std(ddof=1):.1f})'
    rows.append([label,value])
for group,label in [('CN','CN'),('MCI','MCI'),('AD','Clinical AD')]:
    count=int(base.DX_STD.eq(group).sum());rows.append([label+', n (%)',f'{count} ({100*count/n:.1f}%)'])
x=base.MMSE.dropna();rows.append(['Baseline MMSE, mean (SD)',f'{x.mean():.1f} ({x.std(ddof=1):.1f})'])
t1=pd.DataFrame(rows,columns=['Characteristic','ADNI'])
trajectory=pd.read_csv(stage('core')/'results/trajectory_with_q.csv',float_precision='round_trip')
cognition=pd.read_csv(stage('core')/'results/cognition_results.csv',float_precision='round_trip')
records=[]
for group in ['CN','MCI','AD','AD-CN']:
    r=trajectory[trajectory.MODEL.eq('L0')&trajectory.CONTRAST.eq(group)].iloc[0].to_dict()
    r.update(ANALYSIS_TYPE='Annual score slope' if group!='AD-CN' else 'Slope difference',OUTCOME_OR_CONTRAST={'AD':'Clinical AD','AD-CN':'Clinical AD - CN'}.get(group,group),UNIT='ADNI baseline score SD/year')
    records.append(r)
for outcome in ['ADAS13','MMSE']:
    for term,association in [('SCORE_BETWEEN','Between-person'),('SCORE_WITHIN','Within-person')]:
        r=cognition[cognition.OUTCOME.eq(outcome)&cognition.TERM.eq(term)].iloc[0].to_dict()
        r.update(ANALYSIS_TYPE=association,OUTCOME_OR_CONTRAST=outcome,ESTIMATE=r['BETA'],UNIT=outcome+' points/ADNI score SD')
        records.append(r)
t2=pd.DataFrame(records)[['ANALYSIS_TYPE','OUTCOME_OR_CONTRAST','ESTIMATE','CR1_SE','CI_LOW','CI_HIGH','P','Q_BH3','DF','N_PATIENTS','N_MRI','UNIT']]
out=stage('tables')
t1.to_csv(out/'Table1_ADNI.csv',index=False,encoding='utf-8-sig')
t2.to_csv(out/'Table2_ADNI.csv',index=False,encoding='utf-8-sig',float_format='%.17g')
display=t2[['ANALYSIS_TYPE','OUTCOME_OR_CONTRAST']].copy()
display['Estimate (95% CI)']=[f'{r.ESTIMATE:.3f} ({r.CI_LOW:.3f}, {r.CI_HIGH:.3f})' for r in t2.itertuples()]
display['P']=[f'{x:.3g}' for x in t2.P]
page='''<!doctype html><meta charset="utf-8"><title>ADNI main tables</title>
<style>body{font:10pt Arial,sans-serif;max-width:1000px;margin:40px auto;color:#263746}table{border-collapse:collapse;width:100%;margin:16px 0 36px}th,td{text-align:left;padding:7px 10px}thead{border-top:1.2px solid;border-bottom:1px solid}tbody{border-bottom:1px solid}h2{font-size:12pt}@media print{body{margin:0}}</style>'''
page+='<h2>Table 1. ADNI baseline characteristics</h2>'+t1.to_html(index=False,border=0)
page+='<h2>Table 2. ADNI primary results</h2>'+display.to_html(index=False,border=0)
(out/'Main_tables_ADNI.html').write_text(page,encoding='utf-8')
print('Exported ADNI main tables')
