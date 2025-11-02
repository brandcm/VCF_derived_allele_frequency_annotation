from collections import Counter
from pathlib import Path
import argparse
import pysam
import re
import vcfpy

def parse_args():
	parser = argparse.ArgumentParser("Annotate VCF with derived allele frequency (DAF) based on ancestral alleles.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--aa-field", dest="AA_field", type=str, default="AA", help="INFO field name with ancestral allele.")
	parser.add_argument("--fasta", type=str, help="Path to input FASTA file with ancestral alleles.")
	parser.add_argument("--vcf", type=str, required=True, help="Path to input VCF file for annotation. Requires AF INFO field or sample genotypes. DAF is calculated for all samples by default. Use the --samples or --sample-file options to specify individual samples.")
	parser.add_argument("--daf-field", dest="DAF_field", type=str, default="DAF", help="INFO field name for storing the derived allele frequency in the output VCF (default = DAF).")
	parser.add_argument("--daf-field-description", dest="DAF_field_description", type=str, default="Derived allele frequency.", help="INFO field description for derived allele frequencty (default = Derived allele frequency). Enclose in double quotes.")
	parser.add_argument("--af-field", dest="AF_field", type=str, default="AF", help="INFO field name to use for alternate allele frequencies (default = AF).")
	parser.add_argument("--alts-field", dest="ALTs_field", type=str, help="INFO field name to use for alternate alleles when not using VCF ALT.")
	parser.add_argument("--samples", type=str, nargs="+", help="Space-delimited list of sample names for which to calculate DAF.")
	parser.add_argument("--sample-file", dest="sample_file", type=str, help="Path to file with one sample name per line for which to calculate DAF.")
	parser.add_argument("--output", type=str, required=True, help="Path to output file. Will overwrite if it exists.")
	args = parser.parse_args()

	if not args.AA_field and not args.fasta:
		raise ValueError("Either --aa-field or --fasta must be provided.")

	if args.sample_file and args.samples:
		raise ValueError("Cannot use both --samples and --sample-file options. Choose one.")

	return args

def main():
	args = parse_args()

	samples = load_samples_from_file(args.sample_file) if args.sample_file else args.samples or []

	if samples:
		validate_samples(samples, args.vcf)

	with vcfpy.Reader.from_path(args.vcf) as reader:
		add_info_fields(reader.header, args.DAF_field, args.DAF_field_description, include_AA=bool(args.fasta))

		reference = pysam.FastaFile(args.fasta) if args.fasta else None

		with vcfpy.Writer.from_path(args.output, header=reader.header) as writer:
			for variant in reader:
				ancestral_allele = (
					retrieve_ancestral_allele_from_fasta(variant, reference)
					if reference
					else retrieve_ancestral_allele_from_info(variant, args.AA_field)
				)
				DAF = get_variant_derived_allele_frequency(variant, samples, ancestral_allele, args.AF_field, args.ALTs_field)
				if reference:
					variant.INFO['AA'] = ancestral_allele
				variant.INFO[args.DAF_field] = DAF
				writer.write_record(variant)

def load_samples_from_file(sample_file: str) -> list[str]:
	"""Load sample IDs from a text file.

	Reads a text file containing one sample name per line and returns
	the sample names as a list of strings.

	Args:
		sample_file (str):
			Path to the text file containing sample IDs.

	Returns:
		list[str]:
			A list of sample names, one for each line in the input file.

	Example:
		>>> samples = load_samples_from_file("samples.txt")
		>>> print(samples)
		['Sample1', 'Sample2', 'Sample3']
	"""
	return Path(sample_file).read_text().splitlines()

def validate_samples(sample_list: list[str], VCF_path: str) -> None:
	"""Validate that all requested samples are present in a VCF file.

	Compares the provided list of sample names against the samples defined
	in the VCF header. Raises a ValueError if any requested samples are
	missing.

	Args:
		sample_list (list[str]): 
			List of sample names to validate.
		VCF_path (str): 
			Path to the input VCF file.

	Raises:
		ValueError: 
			If one or more of the requested sample names are not found
			in the VCF header.

	Example:
		>>> validate_samples(['SampleA', 'SampleB'], 'variants.vcf.gz')
	"""
	with vcfpy.Reader.from_path(VCF_path) as reader:
		VCF_samples = set(reader.header.samples.names)
		invalid_samples = set(sample_list) - VCF_samples

		if invalid_samples:
			raise ValueError(
				f"The following samples are not found in the VCF: {', '.join(invalid_samples)}"
			)

def add_info_fields(header: "vcfpy.Header", DAF_field: str, DAF_description: str, include_AA: bool = False) -> None:
	"""Add derived allele frequency (DAF) and optional ancestral allele (AA) INFO fields to a VCF header.

	This function updates a VCF header in place by inserting new INFO definitions.
	The DAF field is always added. If `include_aa` is True, an additional INFO field
	for the ancestral allele is also included.

	Args:
		header (vcfpy.Header):
			The VCF header object to which INFO fields will be added.
		DAF_field (str):
			The ID of the derived allele frequency field to add.
		DAF_description (str):
			A description of the DAF field.
		include_AA (bool, optional):
			Whether to include an INFO field for the ancestral allele. Defaults to False.

	Returns:
		None: This function modifies the header object directly.

	Example:
		>>> reader = vcfpy.Reader.from_path("input.vcf")
		>>> add_info_fields(reader.header, "DAF", "Derived allele frequency", include_aa=True)
		>>> writer = vcfpy.Writer.from_path("output.vcf", reader.header)
	"""
	header.add_info_line(vcfpy.OrderedDict([
		('ID', DAF_field),
		('Number', '1'),
		('Type', 'Float'),
		('Description', DAF_description)
	]))
	if include_AA:
		header.add_info_line(vcfpy.OrderedDict([
			('ID', 'AA'),
			('Number', 'A'),
			('Type', 'String'),
			('Description', 'Ancestral allele call from FASTA')
		]))

def retrieve_ancestral_allele_from_fasta(variant: "vcfpy.Record", reference: "pysam.FastaFile") -> str:
	"""Retrieve the ancestral allele for a variant from a FASTA reference.

	Fetches the nucleotide at the variant's position from the provided
	reference FASTA file and returns it in uppercase.

	Args:
		variant (vcfpy.Record):
			Variant record containing chromosome and position information.
		reference (pysam.FastaFile):
			Open FASTA reference object from which to retrieve the allele.

	Returns:
		str:
			Uppercase ancestral allele corresponding to the reference base.

	Example:
		>>> ref = pysam.FastaFile("ancestral.fa")
		>>> allele = retrieve_ancestral_from_fasta(variant, ref)
		>>> print(allele)
		'A'
	"""
	return reference.fetch(variant.CHROM, variant.POS - 1, variant.POS).upper()


def retrieve_ancestral_allele_from_info(variant: "vcfpy.Record", AA_field: str) -> str | None:
	"""Retrieve the ancestral allele for a variant from the VCF INFO field.

	Extracts the ancestral allele from the specified INFO field in a VCF
	record. The value is expected to represent a single nucleotide variant
	(SNV) and should begin with the ancestral base.

	Args:
		variant (vcfpy.Record): 
			Variant record containing INFO annotations.
		AA_field (str): 
			Name of the INFO field containing the ancestral allele value.

	Returns:
		str | None: 
			Uppercase ancestral allele if present and valid, otherwise None.

	Example:
		>>> allele = retrieve_ancestral_from_info(variant, "AA")
		>>> print(allele)
		'G'
	"""
	AA_info = variant.INFO.get(AA_field)
	ancestral_allele = AA_info[0] if isinstance(AA_info, list) and AA_info else AA_info
	if not ancestral_allele or ancestral_allele in {'.', '-', 'N'}:
		return None
	return ancestral_allele[0].upper()

def get_variant_derived_allele_frequency(variant: vcfpy.Record, samples: list[str] | None = None, ancestral_allele: str | None = None, AF_field: str = "AF", ALTs_field: str | None = None) -> float | str:
	"""Calculates the derived allele frequency (DAF) for a single variant.

	Args:
		variant (vcfpy.Record): Variant record containing REF, ALT, INFO, and genotype calls.
		samples (list of str, optional): Subset of sample names to use for genotype-based frequency calculation.
		ancestral_allele (str, optional): Ancestral allele to use; if missing, returns '.'.
		AF_field (str): INFO field containing alternate allele frequencies.
		ALTs_field (str, optional): INFO field containing alternate alleles if different from variant.ALT.

	Returns:
		float or str: Derived allele frequency, or '.' if the ancestral allele is missing.

	Notes:
		- Uses AF field if available; otherwise calculates allele frequencies from genotype calls.
		- Missing genotypes (./.) are ignored; both phased (|) and unphased (/) genotypes are supported.
	"""
	if not ancestral_allele or ancestral_allele in {'.', '-', 'N'}:
		return '.'

	alt_AFs = variant.INFO.get(AF_field)
	alt_alleles = variant.INFO.get(ALTs_field) if ALTs_field else [alt.value for alt in variant.ALT]

	if alt_AFs:
		return compute_derived_allele_frequency(ancestral_allele, variant.REF, alt_AFs, alt_alleles, variant)

	calls = variant.calls if not samples else [call for call in variant.calls if call.sample in samples]
	N_alt_alleles = len(variant.ALT)
	alt_AFs = calculate_allele_frequencies(calls, N_alt_alleles)

	if alt_AFs is None:
		return '.'

	return compute_derived_allele_frequency(ancestral_allele, variant.REF, alt_AFs, alt_alleles, variant)

def compute_derived_allele_frequency(ancestral_allele: str, ref_allele: str, alt_AFs: list[float | str], alt_alleles: list[str], variant: vcfpy.Record | None = None) -> float:
	"""Computes the derived allele frequency (DAF) at a variant site.

	Args:
		ancestral_allele (str): The ancestral allele for the variant.
		ref_allele (str): The reference allele.
		alt_AFs (list of float or str): Frequencies of alternate alleles.
		alt_alleles (list of str): List of alternate alleles.
		variant (vcfpy.Record, optional): The VCF variant record, used only for error messages.

	Returns:
		float: Derived allele frequency. Returns 0.0 for monomorphic sites.

	Notes:
		- If the ancestral allele is equal to the reference, the DAF is the sum of ALT allele frequencies.
		- If the ancestral allele is one of the ALTs, the DAF is 1 minus the frequency of that ALT.
		- If the ancestral allele is not in REF or ALT, the DAF is 1.0 (fully derived).
		- Monomorphic sites (no ALT alleles or missing AF) return 0.0.
	"""
	try:
		alt_alleles_str = [a.value if hasattr(a, "value") else a for a in alt_alleles]

		if alt_alleles_str == ['.'] or not alt_AFs or all(a in {'.', None} for a in alt_AFs):
			return 0.0

		alt_AFs = [float(a) for a in alt_AFs]

		if ancestral_allele == ref_allele:
			return sum(alt_AFs)

		elif ancestral_allele in alt_alleles_str:
			idx = alt_alleles_str.index(ancestral_allele)
			ancestral_allele_AF = alt_AFs[idx] if idx < len(alt_AFs) else alt_AFs[-1]
			return 1 - ancestral_allele_AF

		else:
			return 1.0

	except Exception as e:
		if variant:
			print(f"Error processing variant: CHROM={variant.CHROM}, POS={variant.POS}, REF={variant.REF}, ALT={variant.ALT}")
			print(f"ancestral_allele={ancestral_allele}, alt_alleles={alt_alleles}, alt_afs={alt_AFs}")
		raise e

def calculate_allele_frequencies(calls: list["vcfpy.Call"], N_alt_alleles: int) -> list[float] | None:
	"""Calculate alternate allele frequencies from genotype calls.

	Computes per-allele frequencies based on genotype data for one variant.
	Missing genotypes (e.g., './.') are ignored in the frequency calculation.
	Frequencies are returned in the same order as the ALT alleles listed in
	the VCF record.

	Args:
		calls (list[vcfpy.Call]): 
			List of genotype calls from the VCF for all or selected samples.
		N_alt_alleles (int): 
			Number of alternate alleles at the site.

	Returns:
		list[float] | None: 
			A list of alternate-allele frequencies in the range [0, 1], or 
			``None`` if no non-missing genotypes are available.

	Example:
		>>> calls = [call1, call2, call3]  # vcfpy.Call objects with 'GT' fields
		>>> freqs = calculate_allele_frequencies(calls, 1)
		>>> print(freqs)
		[0.42]

	Notes:
		- Multiallelic sites return one frequency per ALT allele.
		- If all genotypes are missing, returns ``None``.
		- Missing alleles ('.') are not included in the calculation.
	"""
	gts = [call.data.get("GT", None) for call in calls]
	non_missing_gts = [gt for gt in gts if gt is not None and not re.fullmatch(r"\./\.|\.\|\.", gt)]

	if not non_missing_gts:
		return None

	allele_calls = [int(part) for gt in non_missing_gts for part in re.split("/|\\|", gt)]
	allele_counts = Counter(allele_calls)
	total_alleles = sum(allele_counts.values())

	return [allele_counts[i] / total_alleles for i in range(1, N_alt_alleles + 1)] if total_alleles > 0 else []

if __name__ == '__main__':
	main()