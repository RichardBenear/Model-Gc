"""EX-GC26F driver for the PI's PC: four two-knot runs in parallel, then the analysis. Resumable.
Reuses the 26D input arrays c26_p.npy and c26_Fcrop.npy in local_raw/EX-GC26D (hashes in results/EX-GC26D/inputs.json,
from the closed 26C end state) and, on Windows, the MSVC-built kernels (already built for 26D).

  python src/exgc26f_run.py
"""
import hashlib, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; SRC = ROOT / 'src'; OUT = ROOT / 'results/EX-GC26F'; LOGS = ROOT / 'local_raw/EX-GC26F/logs'; IN = ROOT / 'results/EX-GC26D'
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')


def check():
    man = json.loads((IN / 'inputs.json').read_text())
    for name in ('c26_p.npy', 'c26_Fcrop.npy'):
        got = hashlib.sha256((ROOT / 'local_raw/EX-GC26D' / name).read_bytes()).hexdigest()
        assert got == man[name], (name, got, man[name])


def job(d, GN):
    tag = f'd{d:g}_GN{GN:g}'
    if (OUT / f'{tag}.json').exists():
        print('skip', tag, flush=True); return 0
    LOGS.mkdir(parents=True, exist_ok=True); t = time.time()
    with open(LOGS / f'{tag}.log', 'w') as log:
        rc = subprocess.run([sys.executable, str(SRC / 'exgc26f.py'), 'run', '--d', str(d), '--GN', str(GN)], cwd=SRC, env=ENV, stdout=log, stderr=subprocess.STDOUT).returncode
    print(tag, 'rc', rc, f'{time.time() - t:.0f} s wall', flush=True); return rc


if __name__ == '__main__':
    check(); print('inputs verified', flush=True)
    with ThreadPoolExecutor(4) as ex:
        rc = list(ex.map(lambda a: job(*a), [(9.2, 0.01), (9.2, 0), (6.0, 0.01), (6.0, 0)]))
    print('return codes', rc, flush=True)
    subprocess.run([sys.executable, str(SRC / 'exgc26f.py'), 'analyse'], cwd=SRC, env=ENV, check=True)
