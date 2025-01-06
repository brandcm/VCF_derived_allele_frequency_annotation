# VCF Derived Allele Frequency Annotation
This repository contains a Python script (annotate_DAFs.py) and supporting files to annotate the derived allele frequency of genomic positions in a VCF file. Calculation of the derived allele frequency, or DAF, requires ancestral allele information. The current script functions by fetching the ancestral allele from a FASTA file with these data. Typically, I use the [ancestral allele sequences](http://www.ensembl.org/info/genome/compara/ancestral_sequences.html) from the latest version of Ensembl for VCFs with variants mapped to GRCh38/hg38. Two scripts are included in the [download_ancestral_sequence directory](https://github.com/brandcm/VCF_derived_allele_frequency_annotation/blob/main/download_ancestral_sequence) to retrieve either the GRCh37/hg19 or GRCh38/hg38 ancestral sequences. Note that the hg19 sequence was last provided in Ensembl Release 75. I recommend bgzipping and indexing the retrieved file after using the retrieval script: `bgzip -i homo_sapiens_ancestor_*.fa`. These files may need editing for your application depending on the chromosome notation in your VCF ('22' vs 'chr22'). Alternatively, you can edit the sequence headers in the FASTA file directly. Note that unplaced contig IDs do not include the "chrUn" prefix.

The DAF is calculated for sites where the ancestral allele is not missing or ambiguous (see ancestral allele sequence convention below). DAFs are from either 1) allele frequency or 2) genotypes if allele frequency is not present. The resulting VCF contains two new annotations in the INFO field: 1) the ancestral allele call and 2) the derived allele frequency for that genomic position. The annotation script has three required arguments: 1) the path to the FASTA file with ancestral alleles, 2) the path to the input VCF for annotation, and 3) the path to the output file. Therefore, the command line prompt to annotate a VCF would look like this:

```
python3 annotate_DAFs.py --fasta input.fa --vcf input.vcf --out out.vcf
```

Optional arguments include specifying the new INFO field name (default: DAF) and INFO field description (default: Derived allele frequency) using the --DAF_field and --DAF_field_description options, respectively. Enclose any new description in double quotes on the command line.

One can also specify for which samples to calculate the DAF using the --samples or --sample_file options. The former takes a space-delimited string of sample names whereas the latter reads a text file where each line is a sample name. The script will use all samples when calculating a DAF by default.

I have also included an additional script that one can use to update an existing DAF annotation after genotype filtering: update_DAFs.py. The script requires the 1) input VCF (--vcf), 2) the annotation to update (--update_field), and 3) output VCF (--output). One can also use the --samples and --sample_file options above. If genotypes are not provided, this script will use allele count (AC) and allele number (AF) rather than allele frequency (AF) to update the DAF.

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