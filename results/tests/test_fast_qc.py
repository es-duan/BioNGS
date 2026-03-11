"""Tests for Step 0a FastQC runner."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from biongs import step0a_fast_qc as fast_qc


def test_find_input_fastqs(tmp_path):
    experiment = "exp_qc"
    fastq_dir = tmp_path / "input_data" / experiment / f"{experiment}_fastq"
    fastq_dir.mkdir(parents=True)

    (fastq_dir / "a_R1_001.fastq").write_text("@r\nA\n+\n!\n", encoding="utf-8")
    (fastq_dir / "a_R2_001.fastq.gz").write_text("placeholder", encoding="utf-8")
    (fastq_dir / "notes.txt").write_text("ignore me", encoding="utf-8")

    found = fast_qc.find_input_fastqs(experiment, tmp_path)
    assert [path.name for path in found] == ["a_R1_001.fastq", "a_R2_001.fastq.gz"]


def test_cleanup_non_html_fastqc_outputs(tmp_path):
    outdir = tmp_path / "results" / "exp_qc" / "qc"
    outdir.mkdir(parents=True)

    (outdir / "sample_fastqc.html").write_text("html", encoding="utf-8")
    (outdir / "sample_fastqc.zip").write_text("zip", encoding="utf-8")
    extracted = outdir / "sample_fastqc"
    extracted.mkdir()
    (outdir / "other.tmp").write_text("tmp", encoding="utf-8")

    fast_qc.cleanup_non_html_fastqc_outputs(outdir, {"sample_fastqc.html"})

    assert (outdir / "sample_fastqc.html").exists()
    assert not (outdir / "sample_fastqc.zip").exists()
    assert not extracted.exists()
    assert not (outdir / "other.tmp").exists()


def test_run_script_0_5_with_mocked_fastqc(tmp_path, monkeypatch):
    experiment = "exp_qc"
    fastq_dir = tmp_path / "input_data" / experiment / f"{experiment}_fastq"
    fastq_dir.mkdir(parents=True)

    r1 = fastq_dir / "P22R1_R1_001.fastq"
    r2 = fastq_dir / "P22R1_R2_001.fastq"
    r1.write_text("@r\nA\n+\n!\n", encoding="utf-8")
    r2.write_text("@r\nT\n+\n!\n", encoding="utf-8")

    def fake_run_fastqc(fastq_files, outdir):
        outdir.mkdir(parents=True, exist_ok=True)
        for path in fastq_files:
            html_name = f"{fast_qc._strip_fastq_ext(path.name)}_fastqc.html"
            (outdir / html_name).write_text("<html></html>", encoding="utf-8")
            (outdir / f"{fast_qc._strip_fastq_ext(path.name)}_fastqc.zip").write_text(
                "zip",
                encoding="utf-8",
            )

    monkeypatch.setattr(fast_qc, "run_fastqc", fake_run_fastqc)

    html_paths = fast_qc.run_script_0_5(experiment, tmp_path)
    html_names = [path.name for path in html_paths]

    assert html_names == ["P22R1_R1_001_fastqc.html", "P22R1_R2_001_fastqc.html"]
    assert not (tmp_path / "results" / experiment / "qc" / "P22R1_R1_001_fastqc.zip").exists()
    assert not (tmp_path / "results" / experiment / "qc" / "P22R1_R2_001_fastqc.zip").exists()
