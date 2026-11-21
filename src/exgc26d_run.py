"""EX-GC26D driver for the PI's PC: four boosted runs in parallel, then the analysis. Resumable.
Needs local_raw/EX-GC26E1/rest_B1_h0.2_L7.npy and local_raw/EX-GC26D/B1_GN0.01_end.npz (the cloud 26C end state,
hash-checked against results/EX-GC26C/B1_GN0.01.json), and the MSVC-built kernels on Windows.

  python src/exgc26d_run.py
"""
import hashlib, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; SRC = ROOT / 'src'; OUT = ROOT / 'results/EX-GC26D'; LOGS = ROOT / 'local_raw/EX-GC26D/logs'
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')


def check():
    want = json.loads((ROOT / 'results/EX-GC26C/B1_GN0.01.json').read_text())['end_sha256']
    got = hashlib.sha256((ROOT / 'local_raw/EX-GC26D/B1_GN0.01_end.npz').read_bytes()).hexdigest()
    assert got == want, ('26C end state hash mismatch', got, want)
    want = json.loads((ROOT / 'results/EX-GC26E1/static_B1.json').read_text())['sha256']
    got = hashlib.sha256((ROOT / 'local_raw/EX-GC26E1/rest_B1_h0.2_L7.npy').read_bytes()).hexdigest()
    assert got == want, ('26E1 rest state hash mismatch', got, want)


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
