# SpliceAIh
Runs SpliceAI on phased variants

## Installation

Create conda environment
```
mamba create --name spliceaih pip python=3.7 spliceai=1.3.1 pyvcf \
pyensembl pandas tqdm -c bioconda -c anaconda -y
```

Install ENSEMBL release 108
```
mamba activate spliceaih
pyensembl install --release 108 --species homo_sapiens
mamba deactivate
```

Install SpliceAIh
```
git clone https://github.com/CherWeiYuan/SpliceAIh.git
cd SpliceAIh

mamba activate spliceaih
pip install .
mamba deactivate
```

## Test codes
```
mamba activate spliceaih
python3 test.py
mamba deactivate
```


## Download prerequisite files
Download and unzip genome fasta
```
mkdir -p genome
curl -o genome/hg38.fa.gz \
ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
cd genome
gunzip hg38.fa.gz
```


## Trial run on positive control
Run SpliceAIh on two CFTR variants. This is a good positive control because modifying both the T(n) tract and TG(n) tract in intron 8 leads to incorrect splicing of exon 9. The extent of mis-splicing is increased with less T in T(n) tract and more TG repeats in the TG(n) tract.
```
mamba activate spliceaih
mkdir -p sample_output
spliceaih \
--vcf_file test_data/CTFR_T5_TG13.vcf \
--max_read_length 300 \
--fasta genome/Homo_sapiens_assembly38.fasta \
--annotation_file annotations/grch38.txt \
--spliceai_dist 5000 \
--single_variants \
--outdir sample_output/TEST_CTFR_T5_TG13
mamba deactivate
```
