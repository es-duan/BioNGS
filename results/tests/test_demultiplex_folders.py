"""Smoke tests for the BioNGS workflow scripts."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from biongs import step1_demultiplex_folders as demultiplex_folders


def test_demultiplex_folders_smoke(tmp_path, monkeypatch):
    """Create a minimal experiment and verify demultiplex folders are generated."""
    experiment_name = "smoke_experiment"
    input_dir = tmp_path / "input_data" / experiment_name
    input_dir.mkdir(parents=True)

    multiplexing_csv = input_dir / "test_multiplexing_info.csv"
    with multiplexing_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["GW_name", "Population", "R1_index", "R2_index"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "GW_name": "GW_A",
                "Population": "1",
                "R1_index": "ACGT",
                "R2_index": "TGCA",
            }
        )
        writer.writerow(
            {
                "GW_name": "GW_B",
                "Population": "2",
                "R1_index": "GGTT",
                "R2_index": "CCAA",
            }
        )

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        demultiplex_folders,
        "setup_terminal_logging",
        lambda experiment_name, script_name: None,
    )
    monkeypatch.setattr(sys, "argv", ["demultiplex-folders", experiment_name])

    demultiplex_folders.main()

    demultiplex_dir = tmp_path / "results" / experiment_name / "demultiplexing"
    assert demultiplex_dir.is_dir()

    for population in ("P1", "P2"):
        population_dir = demultiplex_dir / population
        assert population_dir.is_dir()
        assert (population_dir / f"{population}_R1.fastq").is_file()
        assert (population_dir / f"{population}_R2.fastq").is_file()