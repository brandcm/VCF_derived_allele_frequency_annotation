# src/annotate_dafs/ancestral.py

"""Functions for retrieving and downloading ancestral sequences."""

from __future__ import annotations

import logging
import re
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path

import pysam
import vcfpy


LOGGER = logging.getLogger(__name__)


GRCH37_ANCESTRAL_URL = (
    "https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/"
    "homo_sapiens_ancestor_GRCh37_e71.tar.bz2"
)

GRCH38_ANCESTRAL_URL = (
    "https://ftp.ensembl.org/pub/release-{release}/fasta/ancestral_alleles/"
    "homo_sapiens_ancestor_GRCh38.tar.gz"
)


def _get_ancestral_fasta_url(
    assembly: str,
    release: int | None = None,
) -> str:
    """Return the Ensembl ancestral sequence archive URL."""
    if assembly == "GRCh37":
        if release is not None:
            raise ValueError(
                "An Ensembl release cannot be specified for GRCh37; "
                "the GRCh37 ancestral sequence is retrieved from release 75."
            )

        return GRCH37_ANCESTRAL_URL

    if assembly == "GRCh38":
        if release is None:
            raise ValueError(
                "An Ensembl release must be specified for GRCh38."
            )

        return GRCH38_ANCESTRAL_URL.format(release=release)

    raise ValueError(
        f"Unsupported assembly {assembly!r}; expected 'GRCh37' or 'GRCh38'."
    )


def _download_file(
    url: str,
    output: Path,
) -> None:
    """Download a file from a URL."""
    LOGGER.info("Downloading %s.", url)

    try:
        with urllib.request.urlopen(url) as response, output.open("wb") as handle:
            shutil.copyfileobj(response, handle)
    except Exception:
        if output.exists():
            output.unlink()
        raise


def _extract_ancestral_archive(
    archive: Path,
    output_directory: Path,
) -> None:
    """Extract an Ensembl ancestral sequence archive."""
    LOGGER.info("Extracting ancestral sequence archive.")

    with tarfile.open(archive, mode="r:*") as tar:
        tar.extractall(output_directory, filter="data")


def _normalize_ensembl_fasta_header(header: str) -> str:
    """Return the sequence name from an Ensembl ancestral FASTA header."""
    header = header.lstrip(">")

    match = re.match(
        r"ANCESTOR_for_(?:chromosome|supercontig):"
        r"[^:]+:([^:]+):\d+:\d+:\d+",
        header,
    )

    if match is None:
        raise ValueError(
            f"Unexpected Ensembl ancestral FASTA header: {header!r}."
        )

    return match.group(1)


def _combine_ancestral_fastas(
    fasta_paths: list[Path],
    output: Path,
) -> None:
    """Combine Ensembl ancestral FASTAs and normalize sequence names."""
    if not fasta_paths:
        raise ValueError(
            "No FASTA files were found in the Ensembl ancestral sequence archive."
        )

    LOGGER.info(
        "Combining %d ancestral sequence FASTA files.",
        len(fasta_paths),
    )

    with output.open("w") as output_handle:
        for fasta_path in sorted(fasta_paths):
            with fasta_path.open() as input_handle:
                for line in input_handle:
                    if line.startswith(">"):
                        sequence_name = _normalize_ensembl_fasta_header(
                            line.strip()
                        )
                        output_handle.write(f">{sequence_name}\n")
                    else:
                        output_handle.write(line)


def download_ancestral_fasta(
    assembly: str,
    output_directory: Path,
    release: int | None = None,
    force: bool = False,
) -> Path:
    """Download, combine, compress, and index an Ensembl ancestral FASTA."""
    output_directory = Path(output_directory)
    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    output = output_directory / f"homo_sapiens_ancestor_{assembly}.fa.gz"
    index = Path(f"{output}.fai")
    gzi_index = Path(f"{output}.gzi")

    if not force:
        existing_files = [
            path
            for path in (output, index, gzi_index)
            if path.exists()
        ]

        if existing_files:
            existing = ", ".join(str(path) for path in existing_files)
            raise FileExistsError(
                f"Output file(s) already exist: {existing}. "
                "Use --force to overwrite them."
            )

    else:
        for path in (output, index, gzi_index):
            if path.exists():
                path.unlink()

    url = _get_ancestral_fasta_url(
        assembly=assembly,
        release=release,
    )

    with tempfile.TemporaryDirectory() as temporary_directory:
        temporary_directory = Path(temporary_directory)

        archive = temporary_directory / Path(url).name
        extracted_directory = temporary_directory / "ancestral"
        combined_fasta = temporary_directory / "ancestral.fa"

        extracted_directory.mkdir()

        _download_file(
            url=url,
            output=archive,
        )

        _extract_ancestral_archive(
            archive=archive,
            output_directory=extracted_directory,
        )

        fasta_paths = [
            path
            for path in extracted_directory.rglob("*.fa")
            if path.is_file()
        ]

        _combine_ancestral_fastas(
            fasta_paths=fasta_paths,
            output=combined_fasta,
        )

        LOGGER.info("BGZF-compressing ancestral FASTA.")

        pysam.tabix_compress(
            str(combined_fasta),
            str(output),
            force=force,
        )

    LOGGER.info("Indexing ancestral FASTA.")
    pysam.faidx(str(output))

    LOGGER.info("Ancestral FASTA written to %s.", output)

    return output


def normalize_chromosome_name(chromosome: str) -> str:
    """Remove a leading 'chr' prefix from a chromosome name."""
    return re.sub(r"^chr", "", chromosome, flags=re.IGNORECASE)


def retrieve_ancestral_allele_from_fasta(
    fasta: pysam.FastaFile,
    chromosome: str,
    position: int,
    reference: str,
    exclude_low_confidence: bool = False,
) -> str | None:
    """Retrieve the ancestral allele from a FASTA file.

    FASTA-based ancestral allele retrieval currently supports SNVs only.
    """
    if len(reference) != 1:
        raise ValueError(
            "FASTA-based ancestral allele retrieval only supports SNVs; "
            f"found REF={reference} at {chromosome}:{position}."
        )

    contigs = set(fasta.references)

    if chromosome in contigs:
        fasta_chromosome = chromosome
    elif normalize_chromosome_name(chromosome) in contigs:
        fasta_chromosome = normalize_chromosome_name(chromosome)
    elif f"chr{normalize_chromosome_name(chromosome)}" in contigs:
        fasta_chromosome = f"chr{normalize_chromosome_name(chromosome)}"
    else:
        raise ValueError(
            f"Chromosome {chromosome!r} was not found in the ancestral FASTA."
        )

    ancestral_allele = fasta.fetch(
        fasta_chromosome,
        position - 1,
        position,
    )

    if ancestral_allele.upper() in {"N", "-", "."}:
        return None

    if ancestral_allele.islower() and exclude_low_confidence:
        return None

    ancestral_allele = ancestral_allele.upper()

    if ancestral_allele != reference.upper():
        LOGGER.debug(
            "Ancestral FASTA allele differs from VCF REF at %s:%d "
            "(VCF REF=%s, FASTA=%s).",
            chromosome,
            position,
            reference,
            ancestral_allele,
        )

    return ancestral_allele


def retrieve_ancestral_allele_from_info(
    variant: vcfpy.Record,
    aa_field: str,
    exclude_low_confidence: bool = False,
) -> str | None:
    """Retrieve the ancestral allele from a VCF INFO field."""
    ancestral_allele = variant.INFO.get(aa_field)

    if ancestral_allele is None:
        raise ValueError(
            f"Missing ancestral allele field {aa_field!r} at "
            f"{variant.CHROM}:{variant.POS}."
        )

    if isinstance(ancestral_allele, list):
        if not ancestral_allele:
            raise ValueError(
                f"Empty ancestral allele field {aa_field!r} at "
                f"{variant.CHROM}:{variant.POS}."
            )

        ancestral_allele = ancestral_allele[0]

    ancestral_allele = str(ancestral_allele).strip()

    if not ancestral_allele:
        raise ValueError(
            f"Empty ancestral allele field {aa_field!r} at "
            f"{variant.CHROM}:{variant.POS}."
        )

    if ancestral_allele.upper() in {"N", "-", "."}:
        return None

    if ancestral_allele.islower() and exclude_low_confidence:
        return None

    return ancestral_allele.upper()