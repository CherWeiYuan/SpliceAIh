# SpliceAIh
Runs [SpliceAI](https://github.com/Illumina/SpliceAI) on phased variants generated from [WhatsHap](https://whatshap.readthedocs.io/en/latest/).

The output reports SpliceAI scores given the haplotype specified in the phased VCF.

WhatsHap can be replaced by any software that generate [pipe notation](https://whatshap.readthedocs.io/en/latest/guide.html#phasing-in-vcfs) in the VCF's GT field.

## How to run SpliceAIh
```
# Run WhatsHap
mamba activate whatshap
whatshap phase \
--indels \
-o sample.phased.vcf \
--reference genome/hg38.fa \
sample.in.vcf \
sample.in.bam
mamba deactivate
```

```
# Run SpliceAIh
mamba activate spliceaih
spliceaih \
--vcf_file sample.phased.vcf \
--max_read_length 200 \
--fasta genome/hg38.fa \
--annotation_file annotations/grch38.txt \
--spliceai_dist 5000 \
--single \
--outdir outdir/sample
mamba deactivate
```

## Installation

Option [1] Set up environment via conda
```
git clone https://github.com/CherWeiYuan/SpliceAIh.git
cd SpliceAIh

mamba create --name spliceaih pip python=3.7 spliceai=1.3.1 pyvcf \
pyensembl pandas tqdm -c bioconda -c anaconda -y

mamba activate spliceaih
pyensembl install --release 108 --species homo_sapiens
pip install .
mamba deactivate
```

Option [2] Install via yaml file
```
git clone https://github.com/CherWeiYuan/SpliceAIh.git
cd SpliceAIh

mamba create --name spliceaih pip python=3.7 -y
mamba env update --name spliceaih --file env/spliceaih.yaml

mamba activate spliceaih
pip install .
mamba deactivate
```

## Test codes
```
# Unzip test data
tar –xvzf test_data.tar.xz

# Run test script
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
