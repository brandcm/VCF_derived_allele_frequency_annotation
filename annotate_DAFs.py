import argparse
import pysam
import re
import vcfpy

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta', type=str, required=True, help='Path to input FASTA file with ancestral alleles')
	parser.add_argument('--vcf', type=str, required=True, help='Path to input VCF file for annotation.')
	parser.add_argument('--output', type=str, required=True, help='Path to output file. Will overwrite if it exists.')
	args = parser.parse_args()
	return args

def main():
	args = parse_args()

	with vcfpy.Reader.from_path(args.vcf) as reader:
		ancestral_allele_INFO = vcfpy.OrderedDict([('ID', 'AA'), ('Number', 'A'), ('Type', 'String'), ('Description', 'Ancestral allele call from Ensembl')])
		DAF_INFO = vcfpy.OrderedDict([('ID', 'DAF'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Derived allele frequency')])
		reader.header.add_info_line(ancestral_allele_INFO)
		reader.header.add_info_line(DAF_INFO)

		with vcfpy.Writer.from_path(args.output, header=reader.header) as writer, pysam.FastaFile(args.fasta) as reference:

			for variant in reader:
				ancestral_allele, DAF = calculate_DAFs(variant, reference)
				variant.INFO['AA'] = ancestral_allele
				variant.INFO['DAF'] = DAF
				writer.write_record(variant)

def calculate_DAFs(variant, reference):
	ancestral_allele = reference.fetch(variant.CHROM, variant.POS - 1, variant.POS)
	ref_allele = str(variant.REF[0])

	DAF = '.'

	if ancestral_allele == '.' or ancestral_allele == '-' or ancestral_allele == 'N':
		alt_alleles = [alt.value for alt in variant.ALT]
		DAF = '.'

	else:
		alt_alleles = [alt.value for alt in variant.ALT]

		if len(alt_alleles) == 0:
			if ancestral_allele.upper() == ref_allele:
				DAF = 0
			else:
				DAF = 1

		elif len(alt_alleles) == 1:
			alt_allele = variant.ALT[0].value

			if 'AF' in variant.INFO:
				alt_AF = variant.INFO['AF'][0]
			else:
				gts = [call.data.get('GT', None) for call in variant.calls]
				non_missing_gts = [gt for gt in gts if gt is not None and not re.match(r'\.|\.', gt)]
				allele_calls = [int(part) for gt in non_missing_gts for part in re.split('/|\|', gt)]
				alt_AF = sum(allele_calls)/len(allele_calls)

			if ancestral_allele.upper() == ref_allele:
				DAF = alt_AF
			elif ancestral_allele.upper() == alt_allele:
				DAF = 1 - alt_AF
			else:
				DAF = 1

		elif len(alt_alleles) > 1:
			if 'AF' in variant.INFO:
				alt_AFs = [AF for AF in variant.INFO['AF']]
				AFs_sum = sum(alt_AFs)

			else:
				gts = [call.data.get('GT', None) for call in variant.calls]
				non_missing_gts = [gt for gt in gts if gt is not None and not re.match(r'\.|\.', gt)]
				allele_calls = [int(part) for gt in non_missing_gts for part in re.split('/|\|', gt)]
				AFs_sum = sum(1 for allele in allele_calls if allele != 0) / len(allele_calls)
			
			if ref_allele == ancestral_allele.upper():
				DAF = AFs_sum

			elif ancestral_allele.upper() in alt_alleles:
				if 'AF' in variant.INFO:
					ancestral_allele_index = alt_alleles.index(ancestral_allele.upper())
					ancestral_alt_AF = alt_AFs[ancestral_allele_index]

				else:
					gts = [call.data.get('GT', None) for call in variant.calls]
					non_missing_gts = [gt for gt in gts if gt is not None and not re.match(r'\.|\.', gt)]
					allele_calls = [int(part) for gt in non_missing_gts for part in re.split('/|\|', gt)]
					ancestral_allele_index = alt_alleles.index(ancestral_allele.upper()) + 1
					ancestral_alt_AF = sum(1 for allele in allele_calls if allele == ancestral_allele_index) / len(allele_calls)
				
				DAF = 1 - ancestral_alt_AF

			elif ancestral_allele.upper() != ref_allele and ancestral_allele.upper() not in alt_alleles:
				DAF = 1

	return ancestral_allele, DAF

if __name__ == '__main__':
	main()