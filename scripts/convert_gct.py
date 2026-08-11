#!/usr/bin/env python3
"""Convert a validated GCT v1.2 count matrix to the RNA TSV contract."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import sys
import tempfile
from pathlib import Path


class GCTValidationError(ValueError):
    """Raised when a GCT v1.2 file does not meet the input contract."""


def split_tab(line: str, line_number: int) -> list[str]:
    if "\t" not in line:
        raise GCTValidationError(f"GCT line {line_number} must be tab-delimited")
    return line.rstrip("\r\n").split("\t")


def parse_dimensions(line: str) -> tuple[int, int]:
    fields = split_tab(line, 2)
    if len(fields) != 2:
        raise GCTValidationError(
            "GCT dimensions line must contain exactly two integers"
        )
    try:
        n_genes, n_samples = (int(value) for value in fields)
    except ValueError as error:
        raise GCTValidationError(
            "GCT dimensions line must contain two integers"
        ) from error
    if n_genes < 1 or n_samples < 1:
        raise GCTValidationError("GCT dimensions must both be positive")
    return n_genes, n_samples


def parse_count(value: str, row_number: int, sample_id: str) -> float:
    try:
        count = float(value)
    except ValueError as error:
        raise GCTValidationError(
            f"GCT row {row_number}, sample {sample_id!r} has a nonnumeric count"
        ) from error
    if not math.isfinite(count):
        raise GCTValidationError(
            f"GCT row {row_number}, sample {sample_id!r} must have a finite count"
        )
    if count < 0:
        raise GCTValidationError(
            f"GCT row {row_number}, sample {sample_id!r} must have a non-negative count"
        )
    return count


def convert_gct(input_path: Path, output_path: Path) -> None:
    """Stream a validated GCT to an atomically replaced counts TSV."""
    temporary_path: Path | None = None
    try:
        with input_path.open("rb") as probe:
            is_gzip = probe.read(2) == b"\x1f\x8b"
        input_stream = (
            gzip.open(input_path, "rt", encoding="utf-8", newline="")
            if is_gzip
            else input_path.open("r", encoding="utf-8", newline="")
        )
        with input_stream as input_handle:
            version_line = input_handle.readline()
            dimensions_line = input_handle.readline()
            header_line = input_handle.readline()
            if not version_line or not dimensions_line or not header_line:
                raise GCTValidationError(
                    "GCT must contain a version, dimensions, and header line"
                )
            if version_line.rstrip("\r\n").lstrip("\ufeff") != "#1.2":
                raise GCTValidationError(
                    "GCT must start with the standard #1.2 header"
                )

            n_genes, n_samples = parse_dimensions(dimensions_line)
            header = split_tab(header_line, 3)
            if header[:2] != ["Name", "Description"]:
                raise GCTValidationError(
                    "GCT header must begin with Name and Description"
                )
            if len(header) != n_samples + 2:
                raise GCTValidationError(
                    "GCT header sample count does not match the declared dimensions"
                )
            sample_ids = header[2:]
            if any(not sample_id for sample_id in sample_ids):
                raise GCTValidationError("GCT sample IDs must not be empty")
            if len(set(sample_ids)) != len(sample_ids):
                raise GCTValidationError("GCT sample IDs must be unique")

            output_path.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                newline="",
                dir=output_path.parent,
                prefix=f".{output_path.name}.",
                suffix=".tmp",
                delete=False,
            ) as temporary:
                temporary_path = Path(temporary.name)
                writer = csv.writer(temporary, delimiter="\t", lineterminator="\n")
                writer.writerow(["gene_id", *sample_ids])
                expected_fields = n_samples + 2
                seen_gene_ids: set[str] = set()
                row_count = 0
                for row_number, line in enumerate(input_handle, start=4):
                    row_count += 1
                    if row_count > n_genes:
                        raise GCTValidationError(
                            "GCT gene count does not match the declared dimensions"
                        )
                    fields = split_tab(line, row_number)
                    if len(fields) != expected_fields:
                        raise GCTValidationError(
                            f"GCT row {row_number} has {len(fields)} fields; "
                            f"expected {expected_fields}"
                        )
                    gene_id, description = fields[:2]
                    if not gene_id or not description:
                        raise GCTValidationError(
                            f"GCT row {row_number} must contain nonempty Name and "
                            "Description values"
                        )
                    if gene_id in seen_gene_ids:
                        raise GCTValidationError(
                            f"GCT contains duplicate gene ID {gene_id!r}"
                        )
                    seen_gene_ids.add(gene_id)
                    counts = [
                        parse_count(value, row_number, sample_id)
                        for sample_id, value in zip(sample_ids, fields[2:])
                    ]
                    writer.writerow([gene_id, *counts])
                if row_count != n_genes:
                    raise GCTValidationError(
                        "GCT gene count does not match the declared dimensions"
                    )

            os.replace(temporary_path, output_path)
            temporary_path = None
    except OSError as error:
        raise GCTValidationError(
            f"could not convert GCT input {input_path} to {output_path}: {error}"
        ) from error
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a validated GCT v1.2 count matrix to gene_id TSV format."
    )
    parser.add_argument("--input", required=True, type=Path, help="input GCT v1.2 file")
    parser.add_argument("--output", required=True, type=Path, help="output counts TSV")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        convert_gct(args.input, args.output)
    except GCTValidationError as error:
        print(f"convert_gct.py: error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
