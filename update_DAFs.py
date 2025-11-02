from collections import Counter
from pathlib import Path
import argparse
import re
import vcfpy

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("--aa-field", dest="AA_field", type=str, default="AA", help="INFO field name with ancestral allele.")
	parser.add_argument("--vcf", type=str, required=True, help="Path to input VCF file for annotation. Requires sample genotypes. DAF is calculated for all samples by default. Use the --samples or --sample-file options to specify individual samples.")
	parser.add_argument("--update-field", dest="update_field", type=str, required=True, help="Name of the existing INFO field to update with recalculated DAF values.")
	parser.add_argument("--samples", type=str, nargs="+", help="Space-delimited list of sample names for which to update DAF.")
	parser.add_argument("--sample-file", dest="sample_file", type=str, help="Path to file with one sample name per line for which to update DAF.")
	parser.add_argument("--output", type=str, required=True, help="Path to output file. Will overwrite if it exists.")
	args = parser.parse_args()
	return args

def main():
	args = parse_args()

	if args.sample_file and args.samples:
		raise ValueError("Cannot use both --samples and --sample-file options. Choose one.")
	samples = load_samples_from_file(args.sample_file) if args.sample_file else args.samples or []

	with vcfpy.Reader.from_path(args.VCF) as reader:
		if args.update_field not in reader.header.info_ids():
			raise ValueError(f"INFO field '{args.update_field}' does not exist in the VCF file.")
		with vcfpy.Writer.from_path(args.output, header=reader.header) as writer:
			for variant in reader:
				DAF = recalculate_DAFs(variant, AA_field=args.AA_field, samples)
				variant.INFO[args.update_field] = DAF if DAF is not None else '.'
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

def recalculate_DAFs(variant: "vcfpy.Record", AA_field: str = "AA", samples: list[str] | None = None) -> float | str:
	"""Recalculate derived allele frequency (DAF) for a variant.

	Computes the derived allele frequency for a given variant, optionally
	restricting the calculation to a subset of samples. If the ancestral
	allele cannot be determined, returns '.'.

	Args:
		variant (vcfpy.Record): 
			Variant record from a VCF file, containing REF, ALT, INFO, and genotype call data.
		AA_field (str):
			Name of the INFO field containing the ancestral allele (e.g., 'AA').
		samples (list[str] | None): 
			Optional list of sample names to filter genotype calls. If None, all samples are used.

	Returns:
		float | str: 
			Derived allele frequency (DAF) as a float, or '.' if undetermined.

	Example:
		>>> DAF = recalculate_DAFs(variant, samples=['Sample1', 'Sample2'])
		>>> print(DAF)
		0.42
	"""
	ref_allele = variant.REF
	alt_alleles = [alt.value for alt in variant.ALT]

	ancestral_allele = retrieve_ancestral_allele_from_info(variant, AA_field)
	if not ancestral_allele:
		return '.'

	calls = [call for call in variant.calls if call.sample in samples] if samples else variant.calls

	alt_AFs = calculate_allele_frequencies(calls, len(alt_alleles)) or []

	if not alt_AFs:
		return '.'

	if ancestral_allele == ref_allele:
		return sum(alt_AFs)
	elif ancestral_allele in alt_alleles:
		idx = alt_alleles.index(ancestral_allele)
		return 1 - alt_AFs[idx]
	else:
		return 1.0

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