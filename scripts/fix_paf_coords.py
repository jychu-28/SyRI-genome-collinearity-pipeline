#!/usr/bin/env python3
"""Restore chunk-level PAF coordinates to chromosome-level coordinates."""

from __future__ import annotations

import argparse
import os
import re

DEFAULT_CHUNK_SIZE = 100_000_000
CHUNK_RE = re.compile(r"(.+?)_chunk(\d+)$")


def load_fai(fai_path: str) -> dict[str, int]:
    """Load FASTA index and return sequence lengths."""
    lengths: dict[str, int] = {}
    with open(fai_path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    return lengths


def build_fai_dict(fai_paths: list[str]) -> dict[str, int]:
    """Merge sequence lengths from multiple .fai files."""
    lengths: dict[str, int] = {}
    for fai_path in fai_paths:
        if not os.path.exists(fai_path):
            raise FileNotFoundError(fai_path)
        lengths.update(load_fai(fai_path))
    return lengths


def fix_name_and_coords(
    name: str,
    length: str,
    start: str,
    end: str,
    fai_dict: dict[str, int],
    chunk_size: int,
) -> tuple[str, str, str, str]:
    """Remove chunk suffix and add chunk offset to coordinates."""
    match = CHUNK_RE.match(name)
    if not match:
        return name, length, start, end

    chrom_name = match.group(1)
    chunk_index = int(match.group(2))
    offset = chunk_index * chunk_size
    chrom_length = str(fai_dict.get(chrom_name, int(length)))

    return (
        chrom_name,
        chrom_length,
        str(int(start) + offset),
        str(int(end) + offset),
    )


def fix_coordinates(
    input_paf: str,
    output_paf: str,
    fai_dict: dict[str, int],
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> int:
    """Rewrite PAF records with chromosome-level coordinates."""
    lines_out = 0
    with open(input_paf, encoding="utf-8") as fin, open(output_paf, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            fields[0], fields[1], fields[2], fields[3] = fix_name_and_coords(
                fields[0], fields[1], fields[2], fields[3], fai_dict, chunk_size
            )
            fields[5], fields[6], fields[7], fields[8] = fix_name_and_coords(
                fields[5], fields[6], fields[7], fields[8], fai_dict, chunk_size
            )

            fout.write("\t".join(fields) + "\n")
            lines_out += 1

    return lines_out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Restore chunk-level PAF coordinates to chromosome-level coordinates."
    )
    parser.add_argument("--input", required=True, help="Input PAF file.")
    parser.add_argument("--output", required=True, help="Output corrected PAF file.")
    parser.add_argument(
        "--fai",
        action="append",
        required=True,
        help="FASTA index (.fai). Provide multiple times for reference and query genomes.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help=f"Chunk size used during splitting. Default: {DEFAULT_CHUNK_SIZE}.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fai_dict = build_fai_dict(args.fai)
    if not fai_dict:
        raise SystemExit("No chromosome lengths were loaded from .fai files.")

    n_lines = fix_coordinates(args.input, args.output, fai_dict, args.chunk_size)
    size_mb = os.path.getsize(args.output) / (1024 * 1024)
    print(f"{args.output}: {n_lines:,} lines, {size_mb:.1f} MB")


if __name__ == "__main__":
    main()

