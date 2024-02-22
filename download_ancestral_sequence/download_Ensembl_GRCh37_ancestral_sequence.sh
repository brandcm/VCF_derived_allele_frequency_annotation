# Colin M. Brand, University of California San Francisco, 02/21/2024

# make temporary directory
temp_dir=$(mktemp -d "./temp_dir.XXXXXX")

# check if tmp dir was created
if [[ ! "$temp_dir" || ! -d "$temp_dir" ]]; then
  echo "Could not create temp dir"
  exit 1
fi

# clean up function
function cleanup {      
  rm -rf "$temp_dir"
  echo "Deleted temporary working directory $temp_dir"
}

trap cleanup EXIT

# change into temp directory
cd "$temp_dir"

# get ancestral alleles
wget http://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2

# decompress directory
tar -xf homo_sapiens_ancestor_GRCh37_e71.tar.bz2

# cd into the decompressed directory
cd homo_sapiens_ancestor_GRCh37_e71

# add a header to the FASTA, necessary to include the first contig in the awk process below
echo "Homo_sapiens_Ancestor_GRCh37 FASTA" > header.txt

# concat the FASTAs and append to the FASTA header
cat *.fa > fastas.fa
cat header.txt fastas.fa > homo_sapiens_ancestor_GRCh37.fa

# rehead and rename FASTA
awk 'NR==FNR{A[$1]=$2; next} NF==2{$2=A[$2]; print ">" $2; next} 1' FS='\t' ../../rename_ancestral_GRCh37_seqs.txt FS='>' homo_sapiens_ancestor_GRCh37.fa | tail -n +2 > rehead_homo_sapiens_ancestor_GRCh37.fa
mv rehead_homo_sapiens_ancestor_GRCh37.fa homo_sapiens_ancestor_GRCh37.fa
mv homo_sapiens_ancestor_GRCh37.fa ../../
cd ../../