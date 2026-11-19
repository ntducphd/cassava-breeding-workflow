"""Run the configured programme-design experiment."""
from pathlib import Path
import argparse
import os
import shutil
import sys

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
from threadpoolctl import threadpool_limits
from uavbreed.programme import run

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--config',type=Path,default=ROOT/'configs/pilot.json')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--rscript',type=Path)
    parser.add_argument('--r-lib',type=Path,default=ROOT/'analysis/env/Rlib')
    args=parser.parse_args()
    rscript=args.rscript or shutil.which('Rscript')
    if not rscript and os.name=='nt':
        matches=sorted(Path('C:/Program Files/R').glob('R-*/bin/Rscript.exe'))
        if matches:rscript=matches[-1]
    if not rscript:raise SystemExit('Rscript is required; provide --rscript')
    with threadpool_limits(limits=1):
        out=run(ROOT,args.config.resolve(),Path(rscript),args.r_lib.resolve(),args.output)
    print(out)

if __name__=='__main__':main()
