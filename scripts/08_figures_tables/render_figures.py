import sys,subprocess
from pathlib import Path
root=Path(__file__).resolve().parent
for n in [3,1,2]:subprocess.run([sys.executable,str(root/f'figure{n}.py')],check=True)
