import sys,csv,random
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from project import stage
root=stage('score')
with (root/'results/OOF_predictions_long.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
def auc(y,s):
    order=sorted(range(len(y)),key=lambda i:s[i]);ranks=[0.]*len(y);i=0
    while i<len(y):
        j=i
        while j+1<len(y) and s[order[j+1]]==s[order[i]]:j+=1
        rank=(i+j+2)/2
        for k in range(i,j+1):ranks[order[k]]=rank
        i=j+1
    n1=sum(y);n0=len(y)-n1
    return (sum(ranks[i] for i in range(len(y)) if y[i]==1)-n1*(n1+1)/2)/(n1*n0)
P={};ids=None;labels={}
for mo in [f'M{i}' for i in range(6)]:
    v={};labs={}
    for r in rows:
        if r['MODEL']==mo:v.setdefault(r['PTID'],[]).append(float(r['P_AD']));labs[r['PTID']]=int(r['Y'])
    if ids is None:ids=list(v);labels=labs
    assert list(v)==ids
    P[mo]=[sum(v[k])/len(v[k]) for k in ids]
y=[labels[k] for k in ids];a0=[i for i,x in enumerate(y) if x==0];a1=[i for i,x in enumerate(y) if x==1]
rng=random.Random(2026102401);boot={mo:[] for mo in P}
for _ in range(2000):
    ix=[rng.choice(a0) for _ in a0]+[rng.choice(a1) for _ in a1];yy=[y[i] for i in ix]
    for mo in P:boot[mo].append(auc(yy,[P[mo][i] for i in ix]))
summary=[]
for mo,s in P.items():
    b=sorted(boot[mo]);summary.append(dict(MODEL=mo,AUC=auc(y,s),CI_LOW=b[50],CI_HIGH=b[1950],N_PATIENTS=len(y),N_CN=len(a0),N_AD=len(a1),N_BOOTSTRAP=2000,SEED=2026102401))
pd.DataFrame(summary).to_csv(root/'results/classifier_performance.csv',index=False,float_format='%.17g',encoding='utf-8-sig')
repeat_summary=[]
for (model,repeat),d in pd.DataFrame(rows).groupby(['MODEL','REPEAT'],sort=False):
    repeat_summary.append(dict(MODEL=model,REPEAT=int(repeat),AUC=auc(d.Y.astype(int).tolist(),d.P_AD.astype(float).tolist())))
pd.DataFrame(repeat_summary).to_csv(root/'results/repeat_AUC.csv',index=False,float_format='%.17g')
