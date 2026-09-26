# src/annotate_dafs/io.py

"""VCF input/output and annotation helpers."""

from __future__ import annotations

from pathlib import Path

import vcfpy


def load_samples_from_file(sample_file: Path) -> list[str]:
    """Load sample IDs from a text file."""
    with sample_file.open() as handle:
        samples = [
            line.strip()
            for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        ]

    if not samples:
        raise ValueError(f"No samples found in {sample_file}.")

    return samples


def validate_samples(
    samples: set[str],
    available_samples: list[str],
) -> None:
    """Validate requested samples against samples present in the VCF."""
    missing = sorted(samples - set(available_samples))

    if missing:
        raise ValueError(
            "The following requested samples were not found in the VCF: "
            + ", ".join(missing)
        )


def add_info_field(
    header: vcfpy.Header,
    field_id: str,
    description: str,
    number: str = "A",
    field_type: str = "Float",
) -> None:
    """Add an INFO field if it does not already exist."""
    existing_fields = {
        record.mapping["ID"]
        for record in header.get_lines("INFO")
    }

    if field_id in existing_fields:
        return

    header.add_info_line(
        {
            "ID": field_id,
            "Number": number,
            "Type": field_type,
            "Description": description,
        }
    )


def parse_alt_alleles(
    variant: vcfpy.Record,
    alts_field: str | None,
) -> list[str]:
    """Retrieve alternate alleles from the VCF or a custom INFO field."""
    if alts_field is None:
        return [alt.value.upper() for alt in variant.ALT]

    alts = variant.INFO.get(alts_field)

    if alts is None:
        raise ValueError(
            f"Missing alternate allele field {alts_field!r} at "
            f"{variant.CHROM}:{variant.POS}."
        )

    if isinstance(alts, str):
        alts = [alts]

    return [str(alt).upper() for alt in alts]