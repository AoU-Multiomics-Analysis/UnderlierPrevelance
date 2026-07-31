#!/usr/bin/env python3
"""Aggregate ClinVar P/LP allele frequencies into per-gene q-squared estimates.

``carrier_freq = 1 - (1 - q)^2`` is the frequency of individuals with at
least one pathogenic allele, including affected homozygotes.
"""

from __future__ import annotations

import argparse
import csv
import glob
import gzip
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Iterator


OUTPUT_COLUMNS = [
    "gene",
    "n_plp_alleles_clinvar",
    "n_plp_alleles_observed",
    "q",
    "incidence",
    "carrier_freq",
]


def open_vcf(path: Path):
    """Open either a plain-text or gzip/bgzip-compressed VCF."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_info(value: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    if value in {"", "."}:
        return fields
    for field in value.split(";"):
        key, separator, field_value = field.partition("=")
        if separator:
            fields[key] = field_value
    return fields


def is_pathogenic_or_likely_pathogenic(info: dict[str, str]) -> bool:
    clinical_significance = info.get("CLNSIG", "")
    normalized = clinical_significance.lower().replace("_", " ").replace("%2f", "/")
    if "conflicting" in normalized:
        return False
    return bool(re.search(r"(^|[|,/])(?:likely )?pathogenic($|[|,/])", normalized))


def genes_from_info(info: dict[str, str]) -> set[str]:
    """Return ClinVar gene symbols from the direct or conventional GENEINFO tag."""
    gene_value = info.get("GENE") or info.get("GENEINFO") or ""
    genes: set[str] = set()
    for gene in re.split(r"[|,]", gene_value):
        symbol = gene.split(":", 1)[0].strip()
        if symbol and symbol != ".":
            genes.add(symbol)
    return genes


def is_small_variant(ref: str, alt: str) -> bool:
    return (
        ref not in {"", "."}
        and alt not in {"", "."}
        and not alt.startswith("<")
        and not any(marker in alt for marker in ("[", "]", "*"))
        and len(ref) <= 50
        and len(alt) <= 50
    )


def parse_clinvar_alleles(path: Path, gene_whitelist: set[str] | None) -> dict[str, set[tuple[str, str, str, str]]]:
    alleles_by_gene: dict[str, set[tuple[str, str, str, str]]] = defaultdict(set)
    with open_vcf(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"malformed VCF record in {path}: {line.rstrip()}")
            chrom, pos, _, ref, alts, _, _, info_text = fields[:8]
            info = parse_info(info_text)
            if not is_pathogenic_or_likely_pathogenic(info):
                continue
            genes = genes_from_info(info)
            if gene_whitelist is not None:
                genes.intersection_update(gene_whitelist)
            for alt in alts.split(","):
                if not is_small_variant(ref, alt):
                    continue
                allele = (chrom, pos, ref, alt)
                for gene in genes:
                    alleles_by_gene[gene].add(allele)
    return alleles_by_gene


def numeric_af(info: dict[str, str], alt_index: int) -> float | None:
    values = info.get("AF", "").split(",")
    if alt_index >= len(values) or values[alt_index] in {"", "."}:
        return None
    try:
        value = float(values[alt_index])
    except ValueError as exc:
        raise ValueError(f"invalid INFO/AF value: {values[alt_index]}") from exc
    if not math.isfinite(value) or value < 0 or value > 1:
        raise ValueError(f"INFO/AF must be finite and between 0 and 1, got {values[alt_index]}")
    return value


def genotype_af(format_text: str, sample_fields: list[str], alt_number: int) -> float | None:
    format_keys = format_text.split(":")
    if "GT" not in format_keys:
        return None
    gt_index = format_keys.index("GT")
    alt_count = 0
    called_alleles = 0
    for sample in sample_fields:
        values = sample.split(":")
        if gt_index >= len(values):
            continue
        genotype = values[gt_index]
        if genotype in {"", ".", "./.", ".|."}:
            continue
        for allele in re.split(r"[|/]", genotype):
            if allele == ".":
                continue
            try:
                allele_number = int(allele)
            except ValueError as exc:
                raise ValueError(f"invalid GT allele {allele!r}") from exc
            called_alleles += 1
            if allele_number == alt_number:
                alt_count += 1
    return alt_count / called_alleles if called_alleles else None


def observed_allele_af(paths: Iterable[Path]) -> dict[tuple[str, str, str, str], float]:
    """Read observed allele AFs, retaining the maximum for duplicate VCF records."""
    af_by_allele: dict[tuple[str, str, str, str], float] = {}
    for path in paths:
        with open_vcf(path) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    raise ValueError(f"malformed VCF record in {path}: {line.rstrip()}")
                chrom, pos, _, ref, alts, _, _, info_text = fields[:8]
                info = parse_info(info_text)
                samples = fields[9:] if len(fields) > 9 else []
                format_text = fields[8] if len(fields) > 8 else ""
                for alt_index, alt in enumerate(alts.split(",")):
                    if not is_small_variant(ref, alt):
                        continue
                    af = numeric_af(info, alt_index)
                    if af is None:
                        af = genotype_af(format_text, samples, alt_index + 1)
                    if af is None:
                        continue
                    allele = (chrom, pos, ref, alt)
                    af_by_allele[allele] = max(af, af_by_allele.get(allele, 0.0))
    return af_by_allele


def load_gene_whitelist(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    genes = {
        line.split("\t", 1)[0].strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    return genes


def format_float(value: float) -> str:
    return format(value, ".15g")


def write_q2_table(
    out_path: Path,
    alleles_by_gene: dict[str, set[tuple[str, str, str, str]]],
    observed_af: dict[tuple[str, str, str, str], float],
) -> None:
    rows: list[dict[str, str | int]] = []
    for gene in sorted(alleles_by_gene):
        clinvar_alleles = alleles_by_gene[gene]
        observed = [allele for allele in clinvar_alleles if observed_af.get(allele, 0.0) > 0.0]
        q = math.fsum(observed_af[allele] for allele in observed)
        if q > 1:
            raise ValueError(
                f"cumulative q exceeds 1 for gene {gene}: {format_float(q)}"
            )
        incidence = q * q
        # User-approved semantics: includes anyone with at least one P/LP allele,
        # including affected homozygotes; it is not heterozygote-only frequency.
        carrier_frequency = 1 - (1 - q) * (1 - q)
        rows.append(
            {
                "gene": gene,
                "n_plp_alleles_clinvar": len(clinvar_alleles),
                "n_plp_alleles_observed": len(observed),
                "q": format_float(q),
                "incidence": format_float(incidence),
                "carrier_freq": format_float(carrier_frequency),
            }
        )
    if not rows:
        raise ValueError("q2 aggregation produced no gene rows")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=(
            "carrier_freq = 1 - (1 - q)^2: frequency of individuals with at least one "
            "pathogenic allele, including affected homozygotes."
        ),
    )
    parser.add_argument("--clinvar", required=True, type=Path, help="ClinVar VCF (plain or gzip-compressed)")
    parser.add_argument("--vcf-glob", required=True, help="Glob identifying deterministic staged cohort VCFs")
    parser.add_argument("--out", required=True, type=Path, help="Destination q-squared TSV")
    parser.add_argument("--genes", type=Path, help="Optional one-gene-per-line whitelist")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if not args.clinvar.is_file():
            raise ValueError(f"ClinVar VCF does not exist: {args.clinvar}")
        vcf_paths = [Path(path) for path in sorted(glob.glob(args.vcf_glob)) if Path(path).is_file()]
        if not vcf_paths:
            raise ValueError(f"VCF glob {args.vcf_glob!r} matched no VCF files")
        if args.genes is not None and not args.genes.is_file():
            raise ValueError(f"gene whitelist does not exist: {args.genes}")
        alleles_by_gene = parse_clinvar_alleles(args.clinvar, load_gene_whitelist(args.genes))
        if not alleles_by_gene:
            raise ValueError(
                "ClinVar contains no eligible P/LP small variants assigned to genes"
            )
        write_q2_table(args.out, alleles_by_gene, observed_allele_af(vcf_paths))
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
