# src/utils/fastq_io.py
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO, Tuple


@dataclass
class FastqRecord:
    header: str
    seq: str
    plus: str
    qual: str


class FastqFormatError(RuntimeError):
    def __init__(self, message: str, read_id: str | None = None, record_no: int | None = None):
        super().__init__(message)
        self.read_id = read_id
        self.record_no = record_no


def open_text_maybe_gz(path: Path, mode: str) -> TextIO:
    # mode: "rt" or "wt"
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", errors="replace")  # type: ignore
    return path.open(mode, encoding="utf-8", errors="replace")


def normalize_read_id(header_line: str) -> str:
    """
    从 '@...' header 中提取 read id，用于配对一致性检查。
    - 去掉开头 '@'
    - 取空格前的 token
    - 去掉结尾的 /1 或 /2（如果有）
    """
    h = header_line.strip()
    if h.startswith("@"):
        h = h[1:]
    token = h.split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def read_fastq(path: Path) -> Iterator[FastqRecord]:
    """
    逐条读取 FASTQ。遇到格式问题抛 FastqFormatError。
    record_no 从 1 开始计数。
    """
    record_no = 0
    with open_text_maybe_gz(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline()
            p = f.readline()
            q = f.readline()
            record_no += 1

            if not (s and p and q):
                rid = normalize_read_id(h) if h else None
                raise FastqFormatError(
                    "Truncated FASTQ record (missing lines).",
                    read_id=rid,
                    record_no=record_no,
                )

            h = h.rstrip("\n")
            s = s.rstrip("\n")
            p = p.rstrip("\n")
            q = q.rstrip("\n")

            rid = normalize_read_id(h)

            if not h.startswith("@"):
                raise FastqFormatError("Invalid FASTQ header (does not start with '@').", rid, record_no)
            if not p.startswith("+"):
                raise FastqFormatError("Invalid FASTQ plus line (does not start with '+').", rid, record_no)
            if len(s) != len(q):
                raise FastqFormatError(
                    f"Sequence/quality length mismatch: len(seq)={len(s)} len(qual)={len(q)}",
                    rid,
                    record_no,
                )

            yield FastqRecord(h, s, p, q)


def write_fastq_record(out: TextIO, rec: FastqRecord) -> None:
    out.write(rec.header + "\n")
    out.write(rec.seq + "\n")
    out.write(rec.plus + "\n")
    out.write(rec.qual + "\n")