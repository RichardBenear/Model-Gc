"""Re-evaluate GC35B's original h=.2 rest states with the GC35D mode cut."""
import argparse
import json
from pathlib import Path
import time

from exgc35d import Knot
from exgc35b_refinement import Knot as OriginalKnot

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("degree", type=int, choices=(1,2))
    b=parser.parse_args().degree
    start=time.process_time()
    knot=Knot(.2,b,6)
    p,_,_,rest=OriginalKnot.rest(knot)
    measured=knot.measure(p)
    result=dict(degree=b,rest=rest,measure=measured,cpu_seconds=time.process_time()-start)
    path=Path(__file__).resolve().parents[1]/"results"/"EX-GC35D"
    path.mkdir(parents=True,exist_ok=True)
    (path/f"baseline_b{b}.json").write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)
