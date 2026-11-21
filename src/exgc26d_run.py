"""EX-GC26D driver for the PI's PC: four boosted runs in parallel, then the analysis. Resumable.
Needs the three input arrays in local_raw/EX-GC26D (hashes in results/EX-GC26D/inputs.json; built in the cloud from
the closed 26C end state and the 26E1 rest state) and, on Windows, the MSVC-built kernels.

  python src/exgc26d_run.py
"""
import hashlib, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; SRC = ROOT / 'src'; OUT = ROOT / 'results/EX-GC26D'; LOGS = ROOT / 'local_raw/EX-GC26D/logs'
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')


def check():
    man = json.loads((OUT / 'inputs.json').read_text())
    for name in ('c26_p.npy', 'c26_Fcrop.npy', 'rest_B1_h0.2_L7.npy'):
        got = hashlib.sha256((ROOT / 'local_raw/EX-GC26D' / name).read_bytes()).hexdigest()
        assert got == man[name], (name, got, man[name])
    assert man['rest_B1_h0.2_L7.npy'] == json.loads((ROOT / 'results/EX-GC26E1/static_B1.json').read_text())['sha256']


def job(u, GN):
    tag = f'u{u:g}_GN{GN:g}'
    if (OUT / f'{tag}.json').exists():
        print('skip', tag, flush=True); return 0
    LOGS.mkdir(parents=True, exist_ok=True); t = time.time()
    with open(LOGS / f'{tag}.log', 'w') as log:
        rc = subprocess.run([sys.executable, str(SRC / 'exgc26d.py'), 'run', '--u', str(u), '--GN', str(GN)], cwd=SRC, env=ENV, stdout=log, stderr=subprocess.STDOUT).returncode
    print(tag, 'rc', rc, f'{time.time() - t:.0f} s wall', flush=True); return rc


if __name__ == '__main__':
    check(); print('inputs verified', flush=True)
    with ThreadPoolExecutor(4) as ex:
        rc = list(ex.map(lambda a: job(*a), [(0.3, 0.01), (0.1, 0.01), (0.3, 0), (0.1, 0)]))
    print('return codes', rc, flush=True)
    subprocess.run([sys.executable, str(SRC / 'exgc26d.py'), 'analyse'], cwd=SRC, env=ENV, check=True)
