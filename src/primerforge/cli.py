from __future__ import annotations

import sys


def main() -> None:
    """
    PrimerForge CLI wrapper.

    Usage:
      primerforge run [args...]
        -> forwards args to primerforge.run_pipeline.main()

    Rationale:
      Your existing run_pipeline.py already has argparse + parameters.
      We keep it as the single source of truth for args.
    """
    if len(sys.argv) <= 1 or sys.argv[1] in ("-h", "--help"):
        print(
            "primerforge: PrimerForge pipeline CLI\n\n"
            "Commands:\n"
            "  run        Run the full pipeline (forwards args to run_pipeline)\n\n"
            "Examples:\n"
            "  primerforge run --exp example --run_tag run1 --r1 ... --r2 ... \\\n"
            "    --multiplex_csv ... --umi_primers_csv ...\n"
        )
        raise SystemExit(0)

    cmd = sys.argv[1]

    if cmd == "run":
        # Strip subcommand and forward the rest to the existing parser in run_pipeline.py
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        from primerforge.run_pipeline import main as run_main
        run_main()
        return

    print(f"Unknown command: {cmd}\nTry: primerforge --help")
    raise SystemExit(2)


if __name__ == "__main__":
    main()
