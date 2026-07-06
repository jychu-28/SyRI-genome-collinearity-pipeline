#!/usr/bin/env python3
"""Convert minimap2 cs:Z: tags in PAF records to cg:Z: CIGAR tags."""

from __future__ import annotations

import argparse
import os
import re
import sys

CS_TOKEN_RE = re.compile(r"(:\d+|\*[a-z][a-z]|[+\-][A-Za-z]+|~[A-Za-z]{2}\d+[A-Za-z]{2}|=[A-Za-z]+)")


def cs_to_cigar(cs: str) -> str:
    """Convert minimap2 cs tag content to an extended CIGAR string."""
    cigar: list[str] = []

    for token in CS_TOKEN_RE.findall(cs):
        op = token[0]
        value = token[1:]

        if op == ":":
            cigar.append(f"{int(value)}=")
        elif op == "=":
            cigar.append(f"{len(value)}=")
        elif op == "*":
            cigar.append("1X")
        elif op == "+":
            cigar.append(f"{len(value)}I")
        elif op == "-":
            cigar.append(f"{len(value)}D")
        elif op == "~":
            match = re.match(r"[A-Za-z]{2}(\d+)[A-Za-z]{2}$", value)
            if match:
                cigar.append(f"{int(match.group(1))}N")

    return "".join(cigar) if cigar else "*"


def process_paf(input_paf: str, output_paf: str, max_lines: int | None = None) -> tuple[int, int, int]:
    """Add cg:Z: tags to PAF records that contain cs:Z: tags."""
    lines_in = 0
    lines_out = 0
    skipped = 0

    with open(input_paf, encoding="utf-8") as fin, open(output_paf, "w", encoding="utf-8") as fout:
        for line in fin:
            lines_in += 1
            if max_lines is not None and lines_in > max_lines:
                break

            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13:
                skipped += 1
                continue

            cs_index = next((i for i in range(12, len(fields)) if fields[i].startswith("cs:Z:")), None)
            if cs_index is None:
                skipped += 1
                continue

            cs = fields[cs_index][5:]
            cigar = cs_to_cigar(cs)
            fields.insert(cs_index + 1, f"cg:Z:{cigar}")

            fout.write("\t".join(fields) + "\n")
            lines_out += 1

            if lines_out % 100_000 == 0:
                print(f"  processed {lines_in:,} lines, wrote {lines_out:,}", file=sys.stderr)

    return lines_in, lines_out, skipped


def run_tests() -> None:
    tests = [
        (":20*ag:11", "20=1X11="),
        (":10+actg:5", "10=4I5="),
        (":10-atcg:5", "10=4D5="),
        (":20*ag*tc:10", "20=1X1X10="),
        ("=ACGT*ag=TT", "4=1X2="),
        ("~gt100ag", "100N"),
    ]
    for cs, expected in tests:
        observed = cs_to_cigar(cs)
        if observed != expected:
            raise AssertionError(f"{cs}: observed {observed}, expected {expected}")
    print("All tests passed.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert PAF cs:Z: tags to cg:Z: CIGAR tags.")
    parser.add_argument("--input", help="Input PAF file.")
    parser.add_argument("--output", help="Output PAF file with cg:Z: tags.")
    parser.add_argument("--max-lines", type=int, default=None, help="Only process the first N lines.")
    parser.add_argument("--test", action="store_true", help="Run built-in tests.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.test:
        run_tests()
        return

    if not args.input or not args.output:
        raise SystemExit("--input and --output are required unless --test is used.")

    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    if os.path.exists(args.input):
        print(f"Input size: {os.path.getsize(args.input) / (1024 ** 3):.2f} GB")

    lines_in, lines_out, skipped = process_paf(args.input, args.output, args.max_lines)
    output_size = os.path.getsize(args.output) / (1024 ** 3)

    print(f"Done: {lines_in:,} lines read, {lines_out:,} written, {skipped:,} skipped")
    print(f"Output size: {output_size:.2f} GB")


if __name__ == "__main__":
    main()

