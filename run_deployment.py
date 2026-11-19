"""Run the selected-genotype management case from frozen programme outputs."""
from pathlib import Path
import argparse
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from uavbreed.deployment import run_deployment


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=ROOT / "analysis/results/validated_pilot")
    parser.add_argument("--output", type=Path, default=ROOT / "analysis/results/genotype_deployment")
    args = parser.parse_args()
    print(run_deployment(ROOT, args.source, args.output))


if __name__ == "__main__":
    main()
