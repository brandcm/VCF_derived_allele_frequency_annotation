import argparse
import re
import vcfpy
from collections import Counter

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--vcf', type=str, required=True, help='Path to input VCF file for annotation. Requires AC and AN INFO field or sample genotypes. DAF is calculated for all samples by default. Use the --samples or --sample_file options to specify individual samples.')
	parser.add_argument('--update_field', type=str, required=True, help='Name of the INFO field to update DAF.')
	parser.add_argument('--samples', type=str, nargs='+', help='Space-delimited list of sample names for which to update DAF.')
	parser.add_argument('--sample_file', type=str, help='Path to file with one sample name per line for which to update DAF.')
	parser.add_argument('--output', type=str, required=True, help='Path to output file. Will overwrite if it exists.')
	args = parser.parse_args()
	return args

def main():
	args = parse_args()

	if args.sample_file and args.samples:
		raise ValueError("Cannot use both --samples and --sample_file options. Choose one.")
	samples = load_samples_from_file(args.sample_file) if args.sample_file else args.samples or []

	with vcfpy.Reader.from_path(args.vcf) as reader:
		if args.update_field not in reader.header.info_ids():
			raise ValueError(f"INFO field '{args.update_field}' does not exist in the VCF file.")
		with vcfpy.Writer.from_path(args.output, header=reader.header) as writer:
			for variant in reader:
				DAF = recalculate_DAFs(variant, samples, disable_AC_calc=bool(samples))
				variant.INFO[args.update_field] = DAF if DAF is not None else '.'
				writer.write_record(variant)

def load_samples_from_file(sample_file):
	"""
	Reads sample IDs from a text file (one per line).
	"""
	with open(sample_file, 'r') as file:
		return [line.strip() for line in file.readlines()]

def retrieve_ancestral_allele(variant):
	"""
	Extract and validate the ancestral allele from a VCF variant record.
	"""
	AA_info = variant.INFO.get('AA')
	ancestral_allele = AA_info[0] if isinstance(AA_info, list) and AA_info else AA_info
	if not ancestral_allele or ancestral_allele in {'.', '-', 'N'}:
		return None
	return ancestral_allele.upper()

def recalculate_DAFs(variant, samples=None, disable_AC_calc=False):
	"""
	Recalculates the derived allele frequency for a given site either after genotype filtering or indicating specific samples.

	Args:
		variant (vcfpy.Record): A variant record from a VCF file, containing REF, ALT, INFO, and genotype call data.
		samples (list[str] | None): Optional list of sample names to filter genotype calls. If None, all samples in the VCF are used.
		disable_AC_calc (bool): If True, disables the use of precomputed AC (allele count) and AN (allele number) fields from the INFO column for DAF calculation, forcing computation from genotype data.
 
	Returns:
		DAF (float | str): The derived allele frequency. If the ancestral allele cannot be determined, '.' is returned.
	"""
	
	ref_allele = variant.REF
	alt_alleles = [alt.value for alt in variant.ALT]

	ancestral_allele = retrieve_ancestral_allele(variant)
	if not ancestral_allele:
		return '.'

	calls = [call for call in variant.calls if call.sample in samples] if samples else variant.calls

	if 'AC' in variant.INFO and 'AN' in variant.INFO and not disable_AC_calc:
		alt_AFs = [AC / variant.INFO['AN'] for AC in variant.INFO['AC']]
	else:
		alt_AFs = calculate_allele_frequencies(calls, len(alt_alleles)) or []

	if not alt_AFs:
		return '.'

	if ancestral_allele == ref_allele:
		return sum(alt_AFs)
	elif ancestral_allele in alt_alleles:
		idx = alt_alleles.index(ancestral_allele)
		return 1 - alt_AFs[idx]
	else:
		return 1

def calculate_allele_frequencies(calls, N_alt_alleles):
	"""
	Calculates allele frequencies for one or multiple alternate alleles at a given position.

	Args:
		calls (list): List of genotype calls from the VCF.
		N_alt_alleles (int): Number of alternate alleles.

	Returns:
		list: A list of allele frequencies, where each index corresponds to an alternate allele.
	"""
	gts = [call.data.get('GT', None) for call in calls]
	non_missing_gts = [gt for gt in gts if not re.match(r'\.|\.\/\.', gt)]
	allele_calls = [int(part) for gt in non_missing_gts for part in re.split('/|\|', gt)]
	
	allele_counts = Counter(allele_calls)
	total_alleles = sum(allele_counts.values())

	return [allele_counts[i] / total_alleles for i in range(1, N_alt_alleles + 1)] if total_alleles > 0 else None

if __name__ == '__main__':
	main()