from pathlib import Path
import os

REPO = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get('ADNI_DATA_DIR', REPO / 'data')).resolve()
RUN = Path(os.environ.get('ADNI_WORK_DIR', REPO / 'analysis')).resolve()

def stage(name):
    p = RUN / name
    p.mkdir(parents=True, exist_ok=True)
    return p
