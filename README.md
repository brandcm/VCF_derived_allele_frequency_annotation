# VCF Derived Allele Frequency Annotation

`annotate-dafs` annotates VCF files with derived allele frequencies (DAFs) using ancestral allele information and alternate allele frequencies.

Ancestral alleles can be provided either:

* in an existing VCF INFO field using `--aa-field`, or
* from an ancestral FASTA sequence using `--fasta`.

Alternate allele frequencies can be obtained from:

1. an INFO field such as `AF` (default), using the VCF ALT alleles;
2. genotypes for explicitly selected samples using `--samples` or `--sample-file`; or
3. an INFO allele-frequency field paired with a corresponding INFO field containing the alternate alleles using `--af-field` and `--alts-field`.

For GRCh38/hg38 VCFs, I typically use the [ancestral allele sequences](http://www.ensembl.org/info/genome/compara/ancestral_sequences.html) provided by Ensembl. Two scripts in the `download_ancestral_sequence` directory can retrieve ancestral sequences for GRCh37/hg19 or GRCh38/hg38. The hg19 sequence was last provided in Ensembl Release 75.

After retrieving an ancestral FASTA, it can be bgzipped and indexed for use with `pysam`:

```bash
bgzip -i homo_sapiens_ancestor_*.fa
```

## Ancestral allele sequence convention

When using an ancestral FASTA from Ensembl, ancestral allele calls follow the conventions described by Ensembl:

| Character | Meaning |
|---|---|
| `A`, `C`, `T`, `G` | High-confidence ancestral allele call, supported by the other two sequences |
| `a`, `c`, `t`, `g` | Low-confidence ancestral allele call, supported by one sequence |
| `N` | Ancestral state is not supported by any other sequence |
| `-` | The extant species contains an insertion at this position |
| `.` | No coverage in the alignment |

Both high- and low-confidence ancestral allele calls are accepted when calculating DAF. Missing or unresolved ancestral states are not annotated.

## Installation

Create a Python environment and install the package:

```bash
conda create -n annotate-dafs python=3.12
conda activate annotate-dafs
pip install .
```

For development, install the package in editable mode:

```bash
pip install -e .
```

## Basic usage

Using an ancestral allele INFO field:

```bash
annotate-dafs \
	--aa-field AA \
	--vcf input.vcf \
	--output output.vcf
```

Using an ancestral FASTA:

```bash
annotate-dafs \
	--fasta ancestral.fa \
	--vcf input.vcf \
	--output output.vcf
```

## Allele frequency sources

By default, `annotate-dafs` uses the `AF` INFO field and the VCF ALT alleles.

A different INFO field can be selected with `--af-field`:

```bash
annotate-dafs \
    --aa-field AA \
    --af-field EAS_AF \
    --vcf input.vcf \
    --output output.vcf
```

### Calculate allele frequencies from selected samples

When `--samples` or `--sample-file` is provided, allele frequencies are calculated directly from the selected samples' genotypes rather than using an INFO AF field.

For example:

```bash
annotate-dafs \
	--aa-field AA \
	--samples sample1 sample2 sample3 \
	--vcf input.vcf \
	--output output.vcf
```

Or provide one sample ID per line in a file:

```bash
annotate-dafs \
	--aa-field AA \
	--sample-file samples.txt \
	--vcf input.vcf \
	--output output.vcf
```

Phased and unphased genotypes are supported, and missing alleles are excluded from the frequency calculation.

### Use alternate alleles from an INFO field

Some VCFs contain allele frequencies corresponding to a set of alternate alleles that differs from the VCF ALT field. In this case, `--alts-field` can be used together with the corresponding `--af-field`.

For example:

```bash
annotate-dafs \
	--aa-field AA \
	--af-field CUSTOM_AF \
	--alts-field CUSTOM_ALTS \
	--vcf input.vcf \
	--output output.vcf
```

The alleles in `CUSTOM_ALTS` and frequencies in `CUSTOM_AF` are paired by position and must contain the same number of values. The VCF REF allele is assumed to be the reference allele for both representations.

When samples are explicitly selected, genotype allele indices are interpreted against the VCF ALT field, so `--alts-field` is not used.

## DAF calculation

DAF is calculated by comparing the ancestral allele with the reference and alternate alleles at each variant.

- When the **ancestral allele is the REF allele**, the ALT allele frequencies are reported as DAF.
- When the **ancestral allele is an ALT allele**, that allele is assigned a DAF of 0, while the frequencies of the other ALT alleles are retained.
- When the **ancestral allele does not match the REF or any ALT allele**, all ALT alleles are assigned a DAF of 1.0.
- When the **ancestral allele is missing**, no DAF annotation is added.

For example, a variant with `REF=A`, `ALT=G`, and `AF=0.10` produces:

| Ancestral allele | DAF |
|---|---:|
| A | 0.10 |
| G | 0.00 |

For a multiallelic variant with `REF=A`, `ALT=G,T`, and `AF=0.10,0.20`:

| Ancestral allele | DAF |
|---|---|
| A | 0.10, 0.20 |
| G | 0.00, 0.20 |
| T | 0.10, 0.00 |

## Command-line options

| Option | Description |
|---|---|
| `--aa-field` | INFO field containing the ancestral allele. Required when `--fasta` is not used. |
| `--fasta` | FASTA file containing the ancestral sequence. Required when `--aa-field` is not used. |
| `--exclude-low-confidence-ancestral` | Treat lowercase ancestral alleles from a FASTA as unresolved and skip DAF annotation. |
| `--vcf` | Input VCF file to annotate. Required. |
| `--output` | Output VCF file. Required and must differ from the input VCF. |
| `--af-field` | INFO field containing alternate allele frequencies. Default: `AF`. |
| `--alts-field` | INFO field containing alternate alleles corresponding to `--af-field`. |
| `--samples` | Sample IDs whose genotypes should be used to calculate allele frequencies. |
| `--sample-file` | File containing one sample ID per line whose genotypes should be used to calculate allele frequencies. |
| `--daf-field` | INFO field name for the derived allele frequency annotation. Default: `DAF`. |
| `--daf-field-description` | Description for the DAF INFO field. Default: `Derived allele frequency`. |
| `--log-level` | Logging level: `DEBUG`, `INFO`, `WARNING`, or `ERROR`. Default: `INFO`. |

Use `annotate-dafs --help` to display the complete command-line interface and available options.

## Requirements and testing

`annotate-dafs` requires Python 3.10 or later and uses [pysam](https://pysam.readthedocs.io/) and [vcfpy](https://vcfpy.readthedocs.io/).

Tests can be run from the repository root with:

```bash
pytest -q
```
