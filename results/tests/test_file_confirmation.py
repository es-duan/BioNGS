"""Tests for Step 0 file confirmation logic."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from biongs import step0_file_confirmation as file_confirmation


def test_file_confirmation_fastq_only(tmp_path, monkeypatch):
    """Only fastq uploads with fastqc installed should allow Step 0a."""
    experiment_name = "exp_a"
    experiment_dir = tmp_path / "input_data" / experiment_name
    fastq_dir = experiment_dir / f"{experiment_name}_fastq"
    fastq_dir.mkdir(parents=True)

    (fastq_dir / "P22R1_R1_001.fastq").write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    (fastq_dir / "P22R1_R2_001.fastq").write_text("@r1\nTGCA\n+\n!!!!\n", encoding="utf-8")

    monkeypatch.setattr(file_confirmation.shutil, "which", lambda name: "/usr/bin/fastqc" if name == "fastqc" else None)

    summary, pairs = file_confirmation.collect_presence(experiment_name, tmp_path)
    report = file_confirmation.build_report(experiment_name, summary, pairs, [])

    assert summary.experiment_dir_exists is True
    assert summary.has_fastq_pairs is True
    assert summary.multiplexing_csv_exists is False
    assert summary.umi_csv_exists is False
    assert summary.fastqc_installed is True
    assert len(pairs) == 1
    assert "Step 0a can run" in report


def test_file_confirmation_all_required_inputs_and_typo_hint(tmp_path, monkeypatch):
    """All expected inputs with fastqc installed should allow full pipeline and include typo suggestions."""
    experiment_name = "exp_b"
    experiment_dir = tmp_path / "input_data" / experiment_name
    fastq_dir = experiment_dir / f"{experiment_name}_fastq"
    fastq_dir.mkdir(parents=True)

    (fastq_dir / "sample_R1.fastq").write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    (fastq_dir / "sample_R2.fastq").write_text("@r1\nTGCA\n+\n!!!!\n", encoding="utf-8")
    (experiment_dir / f"{experiment_name}_multiplexing_info.csv").write_text(
        "GW_name,Population,R1_index,R2_index\nA,1,AAAA,TTTT\n",
        encoding="utf-8",
    )
    (experiment_dir / f"{experiment_name}_UMI_primers.csv").write_text(
        "name,primer\nR1,AAAANNNNNNNNNN\n",
        encoding="utf-8",
    )
    (experiment_dir / f"{experiment_name}_multiplex_info.csv").write_text(
        "placeholder\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(file_confirmation.shutil, "which", lambda name: "/usr/bin/fastqc" if name == "fastqc" else None)

    summary, pairs = file_confirmation.collect_presence(experiment_name, tmp_path)
    suggestions = file_confirmation.find_possible_typos(experiment_name, tmp_path)
    report = file_confirmation.build_report(experiment_name, summary, pairs, suggestions)

    assert summary.experiment_dir_exists is True
    assert summary.has_fastq_pairs is True
    assert summary.multiplexing_csv_exists is True
    assert summary.umi_csv_exists is True
    assert len(pairs) == 1
    assert "Step 6 can run" in report
    assert any(observed == f"{experiment_name}_multiplex_info.csv" for observed, _ in suggestions)


def test_file_confirmation_fastqc_not_installed(tmp_path, monkeypatch):
    """When fastqc is absent the report must warn Step 0a cannot run."""
    experiment_name = "exp_c"
    experiment_dir = tmp_path / "input_data" / experiment_name
    fastq_dir = experiment_dir / f"{experiment_name}_fastq"
    fastq_dir.mkdir(parents=True)

    (fastq_dir / "P22R1_R1_001.fastq").write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    (fastq_dir / "P22R1_R2_001.fastq").write_text("@r1\nTGCA\n+\n!!!!\n", encoding="utf-8")

    monkeypatch.setattr(file_confirmation.shutil, "which", lambda name: None)

    summary, pairs = file_confirmation.collect_presence(experiment_name, tmp_path)
    report = file_confirmation.build_report(experiment_name, summary, pairs, [])

    assert summary.fastqc_installed is False
    assert "NOT INSTALLED" in report
    assert "Step 0a cannot run" in report
