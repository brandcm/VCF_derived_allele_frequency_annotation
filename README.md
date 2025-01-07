# VCF Derived Allele Frequency Annotation
This repository contains scripts and supporting files to annotate the derived allele frequency of genomic positions in a VCF file. Calculation of the derived allele frequency, or DAF, requires ancestral allele information. The current script functions by using existing annotations or fetching the ancestral allele from a FASTA file with these data. Typically, I use the [ancestral allele sequences](http://www.ensembl.org/info/genome/compara/ancestral_sequences.html) from the latest version of Ensembl for VCFs with variants mapped to GRCh38/hg38. Two scripts are included in the [download_ancestral_sequence directory](https://github.com/brandcm/VCF_derived_allele_frequency_annotation/blob/main/download_ancestral_sequence) to retrieve either the GRCh37/hg19 or GRCh38/hg38 ancestral sequences. Note that the hg19 sequence was last provided in Ensembl Release 75. I recommend bgzipping and indexing the retrieved file after using the retrieval script: `bgzip -i homo_sapiens_ancestor_*.fa`. These files may need editing for your application depending on the chromosome notation in your VCF ('22' vs 'chr22'). Alternatively, you can edit the sequence headers in the FASTA file directly. Note that unplaced contig IDs do not include the "chrUn" prefix.

The DAF is calculated for sites where the ancestral allele is not missing or ambiguous (see ancestral allele sequence convention below). DAFs are calculated from either 1) allele frequency or 2) genotypes if allele frequency is not present. The resulting VCF contains one or two new annotations in the INFO field: 1) the ancestral allele call if not already provided and 2) the derived allele frequency for that genomic position. A simple command line prompt to annotate a VCF without ancestral allele calls would look like this:

```
python3 annotate_DAFs.py --fasta input.fa --vcf input.vcf --output out.vcf
```
or the ancestral alleles are annotated:
```
python3 annotate_DAFs.py --AA_field AA --vcf input.vcf --output out.vcf
```
## Optional Arguments
Optional arguments include specifying the new INFO field name (default: DAF) and INFO field description (default: Derived allele frequency) using the `--DAF_field` and `--DAF_field_description` options, respectively. Enclose any new description in double quotes on the command line.

One can also specify for which samples to calculate the DAF using the `--samples` or `--sample_file` options. The former takes a space-delimited string of sample names whereas the latter reads a text file where each line is a sample name. The script will use all samples when calculating a DAF by default.

Finally, one can specify the AF field to use to calculate DAFs using the `--AF_field` option. This is useful when multiple allele frequencies are annotated in a VCF, such as the superpopulation-specific allele frequencies in Thousand Genomes (AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF).

Script arguments are listed below.

```
--AA_field, type=str, default='AA', help=INFO field name with ancestral allele.

--fasta, type=str, help=Path to input FASTA file with ancestral alleles.

--vcf, type=str, required=True, help=Path to input VCF file for annotation. Requires AF INFO field or sample genotypes. DAF is calculated for all samples by default. Use the --samples or --sample_file options to specify individual samples.

--DAF_field, type=str, default='DAF', help=INFO field name for storing the derived allele frequency in the output VCF (default = DAF).

--DAF_field_description, type=str, default='Derived allele frequency.', help=INFO field description for derived allele frequencty (default = Derived allele frequency). Enclose in double quotes.

--AF_field, type=str, default='AF', help=INFO field name to use for alternate allele frequencies (default = AF).

--samples, type=str, nargs='+', help=Space-delimited list of sample names for which to calculate DAF.

--sample_file, type=str, help=Path to file with one sample name per line for which to calculate DAF.

--output, type=str, required=True, help=Path to output file. Will overwrite if it exists.
```

## Updating Existing DAFs
I have also included an additional script that one can use to update an existing DAF annotation, such as after genotype filtering, update_DAFs.py. The script requires the 1) input VCF (`--vcf`), 2) the annotation to update (`--update_field`), and 3) output VCF (`--output`). One can also use the `--samples` and `--sample_file` options above. If genotypes are not provided, this script will use allele count (AC) and allele number (AN) rather than allele frequency (AF) to update the DAF.

Script arguments are listed below.

```
--vcf, type=str, required=True, help=Path to input VCF file for annotation. Requires AF INFO field or sample genotypes. DAF is calculated for all samples by default. Use the --samples or --sample_file options to specify individual samples.

--update_field, type=str, required=True, help=Name of the INFO field to update DAF.

--samples, type=str, nargs='+', help=Space-delimited list of sample names for which to calculate DAF.

--sample_file, type=str, help=Path to file with one sample name per line for which to calculate DAF.

--output, type=str, required=True, help=Path to output file. Will overwrite if it exists.
```

## Requirements & Examples
The annotation script uses two libraries: [pysam](https://pysam.readthedocs.io/en/latest/api.html) and [vcfpy](https://vcfpy.readthedocs.io/en/stable/). I recommend creating a virtual environment with Python and these libraries to run this script. I've included two examples VCFs in the [example_VCFs directory](https://github.com/brandcm/VCF_derived_allele_frequency_annotation/tree/main/example_VCFs). One example has allele frequencies and the other has only genotypes. A Word document in that directory notes the derived allele frequencies that will be output by the program and summarizes the different derived allele frequency scenarios that I considered when writing the program.

Please reach out with any questions or comments: colin.brand@ucsf.edu.

&nbsp;

### Notes:
- Chromosome/contig names must match between the ancestral FASTA sequence and VCF
- Genotypes can be phased or unphased
- Missing genotypes are recognized and derived allele frequency is calculated from the sum of non-missing alleles
- VCFs are assumed to be "unsplit", i.e., multi-allelic positions are recorded on one rather than multiple lines
- Ancestral allele calls are largely missing in telomeric regions of the genome
- DAFs are calculated for both low- and high-confidence calls (see sequence convention below)
- DAFs are calculated for positions included for any positions included in a VCF that are fixed for the reference allele (DAF = 0 if the ancestral allele is the reference allele, otherwise DAF = 1)

&nbsp;

### Ancestral Allele Sequence Convention (per Ensembl):  
ACTG : high-confidence call, ancestral state supproted by the other two sequences  
actg : low-confindence call, ancestral state supported by one sequence only  
N    : failure, the ancestral state is not supported by any other sequence  
\-    : the extant species contains an insertion at this postion  
.    : no coverage in the alignment