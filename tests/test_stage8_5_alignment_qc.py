from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"

sys.path.insert(0, str(SRC_DIR))

from primerforge.pipeline.stage8_5_alignment_qc import stage8_5_alignment_qc


def main():

    exp = "example"
    run_tag = "run2"

    results_dir = Path("results") / exp

    stage8_5_alignment_qc(
        exp=exp,
        run_tag=run_tag,
        results_exp_dir=results_dir,
    )


if __name__ == "__main__":
    main()