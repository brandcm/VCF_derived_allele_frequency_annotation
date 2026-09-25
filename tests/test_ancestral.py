# tests/test_ancestral.py

from pathlib import Path

import pysam
import pytest
import vcfpy

from annotate_dafs.ancestral import (
	normalize_chromosome_name,
	retrieve_ancestral_allele_from_fasta,
	retrieve_ancestral_allele_from_info,
)


DATA_DIR = Path(__file__).parent / "data"


def make_variant(
	chromosome: str = "chr22",
	position: int = 10,
	reference: str = "C",
	alternate: str = "A",
	info: dict | None = None,
) -> vcfpy.Record:
	"""Create a VCF record for testing."""
	return vcfpy.Record(
		CHROM=chromosome,
		POS=position,
		ID=["."],
		REF=reference,
		ALT=[
			vcfpy.Substitution(
				type_="SNV",
				value=alternate,
			)
		],
		QUAL=30,
		FILTER=["PASS"],
		INFO=info or {},
	)


def test_normalize_chromosome_name():
	"""Leading chr prefixes are removed case-insensitively."""
	assert normalize_chromosome_name("chr22") == "22"
	assert normalize_chromosome_name("CHR22") == "22"
	assert normalize_chromosome_name("22") == "22"


def test_retrieve_ancestral_allele_from_info():
	"""Ancestral alleles are retrieved and normalized from INFO."""
	variant = make_variant(info={"AA": "c"})

	ancestral_allele = retrieve_ancestral_allele_from_info(
		variant,
		"AA",
	)

	assert ancestral_allele == "C"


def test_retrieve_ancestral_allele_from_info_list():
	"""The first value is used when the INFO field contains a list."""
	variant = make_variant(info={"AA": ["a"]})

	ancestral_allele = retrieve_ancestral_allele_from_info(
		variant,
		"AA",
	)

	assert ancestral_allele == "A"


@pytest.mark.parametrize("ancestral_state", ["N", "n", "-", "."])
def test_retrieve_ancestral_allele_from_info_unresolved(ancestral_state):
	"""Unresolved ancestral states return None."""
	variant = make_variant(info={"AA": ancestral_state})

	ancestral_allele = retrieve_ancestral_allele_from_info(
		variant,
		"AA",
	)

	assert ancestral_allele is None


def test_retrieve_ancestral_allele_from_info_exclude_low_confidence():
	"""Lowercase ancestral alleles can be treated as unresolved."""
	variant = make_variant(info={"AA": "a"})

	ancestral_allele = retrieve_ancestral_allele_from_info(
		variant,
		"AA",
		exclude_low_confidence=True,
	)

	assert ancestral_allele is None


def test_retrieve_ancestral_allele_from_info_missing():
	"""A missing ancestral allele raises an informative error."""
	variant = make_variant()

	with pytest.raises(
		ValueError,
		match="Missing ancestral allele field",
	):
		retrieve_ancestral_allele_from_info(
			variant,
			"AA",
		)


def test_retrieve_ancestral_allele_from_info_empty():
	"""An empty ancestral allele raises an informative error."""
	variant = make_variant(info={"AA": ""})

	with pytest.raises(
		ValueError,
		match="Empty ancestral allele field",
	):
		retrieve_ancestral_allele_from_info(
			variant,
			"AA",
		)


def test_retrieve_ancestral_allele_from_fasta():
	"""The ancestral allele is retrieved from the FASTA."""
	fasta = pysam.FastaFile(DATA_DIR / "ancestral.fa")

	try:
		ancestral_allele = retrieve_ancestral_allele_from_fasta(
			fasta,
			"chr22",
			10,
			"C",
		)
	finally:
		fasta.close()

	assert ancestral_allele == "C"


def test_retrieve_ancestral_allele_from_fasta_without_chr_prefix():
	"""FASTA retrieval handles chromosome names without a chr prefix."""
	fasta = pysam.FastaFile(DATA_DIR / "ancestral.fa")

	try:
		ancestral_allele = retrieve_ancestral_allele_from_fasta(
			fasta,
			"22",
			10,
			"C",
		)
	finally:
		fasta.close()

	assert ancestral_allele == "C"


def test_retrieve_ancestral_allele_from_fasta_non_snv():
	"""FASTA retrieval rejects non-SNV reference alleles."""
	fasta = pysam.FastaFile(DATA_DIR / "ancestral.fa")

	try:
		with pytest.raises(
			ValueError,
			match="only supports SNVs",
		):
			retrieve_ancestral_allele_from_fasta(
				fasta,
				"chr22",
				10,
				"AT",
			)
	finally:
		fasta.close()

def test_retrieve_ancestral_allele_from_fasta_exclude_low_confidence(
	tmp_path,
):
	"""Lowercase FASTA alleles can be treated as unresolved."""
	fasta_path = tmp_path / "ancestral.fa"
	fasta_path.write_text(">chr22\naAAAAAAAA\n")
	pysam.faidx(str(fasta_path))

	fasta = pysam.FastaFile(fasta_path)

	try:
		ancestral_allele = retrieve_ancestral_allele_from_fasta(
			fasta,
			"chr22",
			1,
			"A",
			exclude_low_confidence=True,
		)
	finally:
		fasta.close()

	assert ancestral_allele is None


@pytest.mark.parametrize("ancestral_state", ["N", "n", "-", "."])
def test_retrieve_ancestral_allele_from_fasta_unresolved(
	tmp_path,
	ancestral_state,
):
	"""Unresolved FASTA states return None."""
	fasta_path = tmp_path / "ancestral.fa"
	fasta_path.write_text(f">chr22\n{ancestral_state}\n")
	pysam.faidx(str(fasta_path))

	fasta = pysam.FastaFile(fasta_path)

	try:
		ancestral_allele = retrieve_ancestral_allele_from_fasta(
			fasta,
			"chr22",
			1,
			"A",
		)
	finally:
		fasta.close()

	assert ancestral_allele is None