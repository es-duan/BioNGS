# test_stage7.py
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_DIR))

from primerforge.pipeline.stage7_family_strategy import stage7_family_strategy


def main():
    exp = "example"
    run_tag = "run2"

    results_exp_dir = Path("results") / exp
    metrics_dir = results_exp_dir / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)

    payload = stage7_family_strategy(
        exp=exp,
        run_tag=run_tag,
        results_exp_dir=results_exp_dir,
        metrics_dir=metrics_dir,
        trimmed_version=None,   # 自动选最新 v*
    )

    print("\nStage 7 finished.")
    print("Returned keys:", payload.keys())
    print("Populations:", list(payload["per_population"].keys()))


if __name__ == "__main__":
    main()