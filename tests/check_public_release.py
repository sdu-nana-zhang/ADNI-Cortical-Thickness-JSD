from pathlib import Path
import re,json

root=Path(__file__).resolve().parents[1]
excluded={'data','analysis','local','.git','.venv','__pycache__'}
allowed={'.py','.R','.json','.csv','.md','.yml','.txt'}
findings=[];count=0
for p in root.rglob('*'):
    if not p.is_file() or any(x in excluded for x in p.relative_to(root).parts):continue
    count+=1
    if p.name!='.gitignore' and p.suffix not in allowed:findings.append([str(p.relative_to(root)),'Unexpected file type']);continue
    text=p.read_text(encoding='utf-8-sig')
    for label,pattern in [
        ('Personal absolute path',r'[A-Za-z]:[\\/](?!/)|'+'/'+'Users/|'+'/'+'home/'),
        ('Participant identifier',r'\b\d{3}_S_\d{4}\b|\bNACC\d{4,}\b|\bOAS2_\d{3,}\b'),
        ('Access credential',r'gh[pousr]_[A-Za-z0-9]{24,}|sk-[A-Za-z0-9_-]{24,}|AKIA[A-Z0-9]{16}')]:
        if re.search(pattern,text):findings.append([str(p.relative_to(root)),label])
    if p.suffix=='.csv' and p.parent.name not in {'config','schemas','docs'}:findings.append([str(p.relative_to(root)),'CSV outside dictionary/schema/report locations'])
if findings:raise SystemExit(json.dumps(findings,indent=2))
print('PASS:',count,'public source/configuration files; no detected participant identifiers, personal paths, credentials or binary data')
