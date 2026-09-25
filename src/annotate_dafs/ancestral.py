# src/annotate_dafs/ancestral.py

"""Functions for retrieving ancestral alleles."""

from __future__ import annotations

import logging
import re

import pysam
import vcfpy


LOGGER = logging.getLogger(__name__)


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