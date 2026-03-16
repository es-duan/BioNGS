from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"

sys.path.insert(0, str(SRC_DIR))

from primerforge.pipeline.stage8_alignment import stage8_alignment


def main():

    exp = "example"
    run_tag = "run2"

    results_dir = Path("results") / exp

    reference_index = Path(
        "input_data/reference_docs/rpoB_bowtie_index/rpoBsequence_index"
    )

    stage8_alignment(
        exp=exp,
        run_tag=run_tag,
        results_exp_dir=results_dir,
        reference_index=reference_index,
    )


if __name__ == "__main__":
    main()