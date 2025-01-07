import argparse
import pysam
import re
import vcfpy
from collections import Counter

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--AA_field', type=str, default='AA', help='INFO field name with ancestral allele.')
	parser.add_argument('--fasta', type=str, help='Path to input FASTA file with ancestral alleles.')
	parser.add_argument('--vcf', type=str, required=True, help='Path to input VCF file for annotation. Requires AF INFO field or sample genotypes. DAF is calculated for all samples by default. Use the --samples or --sample_file options to specify individual samples.')
	parser.add_argument('--DAF_field', type=str, default='DAF', help='INFO field name for storing the derived allele frequency in the output VCF (default = DAF).')
	parser.add_argument('--DAF_field_description', type=str, default='Derived allele frequency.', help='INFO field description for derived allele frequencty (default = Derived allele frequency). Enclose in double quotes.')
	parser.add_argument('--AF_field', type=str, default='AF', help='INFO field name to use for alternate allele frequencies (default = AF).')
	parser.add_argument('--samples', type=str, nargs='+', help='Space-delimited list of sample names for which to calculate DAF.')
	parser.add_argument('--sample_file', type=str, help='Path to file with one sample name per line for which to calculate DAF.')
	parser.add_argument('--output', type=str, required=True, help='Path to output file. Will overwrite if it exists.')
	args = parser.parse_args()

	if not args.AA_field and not args.fasta:
		raise ValueError("Either --AA_field or --fasta must be provided.")
	return args

def main():
	args = parse_args()

	if args.sample_file and args.samples:
		raise ValueError("Cannot use both --samples and --sample_file options. Choose one.")
	samples = load_samples_from_file(args.sample_file) if args.sample_file else args.samples or []

	with vcfpy.Reader.from_path(args.vcf) as reader:
		if args.fasta:
			ancestral_allele_INFO = vcfpy.OrderedDict([('ID', 'AA'), ('Number', 'A'), ('Type', 'String'), ('Description', 'Ancestral allele call from Ensembl')])
			reader.header.add_info_line(ancestral_allele_INFO)
		DAF_INFO = vcfpy.OrderedDict([('ID', args.DAF_field), ('Number', '1'), ('Type', 'Float'), ('Description', args.DAF_field_description)])
		reader.header.add_info_line(DAF_INFO)

		with vcfpy.Writer.from_path(args.output, header=reader.header) as writer:
			if args.fasta:  # Only open the FASTA file if --fasta is provided
				with pysam.FastaFile(args.fasta) as reference:
					for variant in reader:
						ancestral_allele = retrieve_ancestral_from_fasta(variant, reference) if args.fasta else retrieve_ancestral_from_info(variant, args.AA_field)
						DAF = calculate_DAFs(variant, samples, ancestral_allele)
						variant.INFO['AA'] = ancestral_allele
						variant.INFO[args.DAF_field] = DAF
						writer.write_record(variant)
			else:
				for variant in reader:
					ancestral_allele = retrieve_ancestral_from_info(variant, args.AA_field)
					DAF = calculate_DAFs(variant, samples, ancestral_allele)
					variant.INFO[args.DAF_field] = DAF
					writer.write_record(variant)

def load_samples_from_file(sample_file):
	"""
	Reads sample IDs from a text file (one per line).
	"""
	with open(sample_file, 'r') as file:
		return [line.strip() for line in file.readlines()]

def retrieve_ancestral_from_fasta(variant, reference):
	"""
	Retrieve the ancestral allele for a given variant from a FASTA reference.
	"""
	return reference.fetch(variant.CHROM, variant.POS - 1, variant.POS).upper()


def retrieve_ancestral_from_info(variant, AA_field):
	"""
	Retrieves the ancestral allele for a given variant from an INFO field.
	"""
	AA_info = variant.INFO.get(AA_field)
	ancestral_allele = AA_info[0] if isinstance(AA_info, list) and AA_info else AA_info
	if not ancestral_allele or ancestral_allele in {'.', '-', 'N'}:
		return None
	return ancestral_allele.upper()

def calculate_DAFs(variant, samples=None, ancestral_allele=None, AF_field='AF'):
	"""
	Calculates the derived allele frequency for a given site.

	Args:
		variant (vcfpy.Record): A variant record from a VCF file, containing REF, ALT, INFO, and genotype call data.
		samples (list): List of sample names for which to calculate DAF.
		ancestral_allele (str): The ancestral allele to be used when calculating DAF.
		AF_field (str): INFO field name to fetch alternate allele frequencies.

	Returns:
		DAF (float | str): The derived allele frequency. If the ancestral allele cannot be determined, '.' is returned.
	"""
	ref_allele = variant.REF
	alt_alleles = [alt.value for alt in variant.ALT]

	if not ancestral_allele or ancestral_allele in {'.', '-', 'N'}:
		return '.'

	calls = [call for call in variant.calls if call.sample in samples] if samples else variant.calls

	if AF_field in variant.INFO:
		alt_AFs = [AF for AF in variant.INFO[AF_field]]
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

	return [allele_counts[i] / total_alleles for i in range(1, N_alt_alleles + 1)] if total_alleles > 0 else []

if __name__ == '__main__':
	main()