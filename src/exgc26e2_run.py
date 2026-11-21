"""EX-GC26E2 driver for the PI's PC: validate, relax B = 1 and B = 2 at h = 0.15, then the six fall runs in parallel,
then the analysis. Resumable: finished records are skipped. Workers from MODELG_WORKERS (default 6, at most 8).

  python src/exgc26e2_run.py
"""
import os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; SRC = ROOT / 'src'; OUT = ROOT / 'results/EX-GC26E2'; LOGS = ROOT / 'local_raw/EX-GC26E2/logs'
WORKERS = min(8, int(os.environ.get('MODELG_WORKERS', '6')))
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')


def job(args, tag):
    rec = OUT / f'{tag}.json'
    if rec.exists():
        print('skip', tag, flush=True); return 0
    LOGS.mkdir(parents=True, exist_ok=True); t = time.time()
    with open(LOGS / f'{tag}.log', 'w') as log:
        rc = subprocess.run([sys.executable, str(SRC / 'exgc26e1.py')] + args + ['--h', '0.15', '--cpu-limit', '0'], cwd=SRC, env=ENV, stdout=log, stderr=subprocess.STDOUT).returncode
    print(tag, 'rc', rc, f'{time.time() - t:.0f} s wall', flush=True)
    return rc


def main():
    print('workers', WORKERS, flush=True)
    assert job(['validate'], 'validation') == 0
    with ThreadPoolExecutor(2) as ex:
        rc = list(ex.map(lambda a: job(*a), [(['static', '--B', '1', '--maxiter', '300'], 'static_B1'), (['static', '--B', '2', '--maxiter', '600'], 'static_B2')]))
    assert rc == [0, 0], rc
    runs = [(['run', '--B', B, '--mode', m], f'B{B}_{m}') for B, m in (('1', 'energy'), ('2', 'energy'), ('2', 'null'), ('2', 'control'), ('1', 'control'), ('1', 'null'))]
    with ThreadPoolExecutor(WORKERS) as ex:
        rc = list(ex.map(lambda a: job(*a), runs))
    print('run return codes', rc, flush=True)
    subprocess.run([sys.executable, str(SRC / 'exgc26e1_analysis.py'), 'E2'], cwd=SRC, env=ENV, check=True)


if __name__ == '__main__':
    main()
