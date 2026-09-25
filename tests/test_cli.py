from pathlib import Path

import pysam
import pytest
import vcfpy

from annotate_dafs.cli import main


DATA_DIR = Path(__file__).parent / "data"


def read_records(path: Path) -> list[vcfpy.Record]:
	"""Read all records from a VCF file."""
	reader = vcfpy.Reader.from_path(path)
	return list(reader)


def test_cli_annotates_daf_from_info(tmp_path):
	"""CLI annotates DAF using ancestral alleles and INFO/AF."""
	input_vcf = DATA_DIR / "allele_frequencies.vcf"
	output_vcf = tmp_path / "output.vcf"

	exit_code = main(
		[
			"--aa-field",
			"AA",
			"--vcf",
			str(input_vcf),
			"--output",
			str(output_vcf),
		]
	)

	assert exit_code == 0

	reader = vcfpy.Reader.from_path(output_vcf)

	assert "DAF" in reader.header.info_ids()

	records = list(reader)

	assert records[0].INFO["DAF"] == [0.01]
	assert records[1].INFO["DAF"] == [0.1]
	assert records[2].INFO["DAF"] == [0.0]
	assert records[3].INFO["DAF"] == [0.0]
	assert records[4].INFO["DAF"] == [0.0, 0.03]
	assert records[5].INFO["DAF"] == [0.08, 0.09]


def test_cli_uses_selected_samples(tmp_path):
	"""CLI uses genotype-derived AF when samples are explicitly selected."""
	input_vcf = DATA_DIR / "genotypes.vcf"
	output_vcf = tmp_path / "output.vcf"

	exit_code = main(
		[
			"--aa-field",
			"AA",
			"--samples",
			"sample4",
			"--vcf",
			str(input_vcf),
			"--output",
			str(output_vcf),
		]
	)

	assert exit_code == 0

	records = read_records(output_vcf)

	assert records[0].INFO["DAF"] == [0.5]


def test_cli_uses_fasta(tmp_path):
	"""CLI annotates DAF using ancestral alleles from a FASTA."""
	input_vcf = DATA_DIR / "cli_input.vcf"
	ancestral_fasta = DATA_DIR / "ancestral.fa"
	output_vcf = tmp_path / "output.vcf"

	exit_code = main(
		[
			"--fasta",
			str(ancestral_fasta),
			"--vcf",
			str(input_vcf),
			"--output",
			str(output_vcf),
		]
	)

	assert exit_code == 0

	records = read_records(output_vcf)

	assert records[0].INFO["DAF"] == [1.0]


def test_cli_excludes_low_confidence_ancestral(tmp_path):
	"""CLI skips DAF annotation for lowercase ancestral alleles when requested."""
	input_vcf = DATA_DIR / "cli_input.vcf"
	ancestral_fasta = tmp_path / "ancestral.fa"
	output_vcf = tmp_path / "output.vcf"

	ancestral_fasta.write_text(">chr22\nAAAAAAAAAa\n")
	pysam.faidx(str(ancestral_fasta))

	exit_code = main(
		[
			"--fasta",
			str(ancestral_fasta),
			"--exclude-low-confidence-ancestral",
			"--vcf",
			str(input_vcf),
			"--output",
			str(output_vcf),
		]
	)

	assert exit_code == 0

	records = read_records(output_vcf)

	assert "DAF" not in records[0].INFO


def test_cli_adds_daf_header(tmp_path):
	"""CLI adds the DAF INFO definition to the output VCF header."""
	input_vcf = DATA_DIR / "cli_input.vcf"
	output_vcf = tmp_path / "output.vcf"

	exit_code = main(
		[
			"--aa-field",
			"AA",
			"--vcf",
			str(input_vcf),
			"--output",
			str(output_vcf),
		]
	)

	assert exit_code == 0

	reader = vcfpy.Reader.from_path(output_vcf)

	assert "DAF" in reader.header.info_ids()


def test_cli_uses_custom_alts_field(tmp_path):
	"""CLI uses a custom ALT field with the corresponding INFO AF field."""
	input_vcf = DATA_DIR / "cli_input.vcf"
	output_vcf = tmp_path / "custom_alts_output.vcf"

	exit_code = main(
		[
			"--aa-field",
			"AA",
			"--vcf",
			str(input_vcf),
			"--af-field",
			"CUSTOM_AF",
			"--alts-field",
			"CUSTOM_ALTS",
			"--output",
			str(output_vcf),
		]
	)

	assert exit_code == 0

	records = read_records(output_vcf)

	# VCF ALT is G, but CUSTOM_ALTS is T with CUSTOM_AF=0.4.
	# AA=A, so T is derived and should retain its frequency.
	assert records[0].INFO["DAF"] == pytest.approx([0.4])