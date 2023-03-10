<img src="https://github.com/CherWeiYuan/SpliceAIh/blob/main/logo/Original%20on%20Transparent.png?raw=true" width=600 height=170>

Runs [SpliceAI](https://github.com/Illumina/SpliceAI) on phased variants generated from [WhatsHap](https://whatshap.readthedocs.io/en/latest/).

The output reports SpliceAI scores given the haplotype specified in the phased VCF.

WhatsHap can be replaced by any software that generate [pipe notation](https://whatshap.readthedocs.io/en/latest/guide.html#phasing-in-vcfs) in the VCF's GT field.

## How to run SpliceAIh
Run WhatsHap first to generate phased vcf
```
mamba activate whatshap
whatshap phase \
--indels \
-o sample.phased.vcf \
--reference genome/hg38.fa \
sample.in.vcf \
sample.in.bam
mamba deactivate
```

SpliceAIh takes the phased vcf as input and produces a tab-separated file (TSV) in the output directory
```
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
Note: Installation instructions assumes GRCh38/ hg38 reference genome is used to generate your variant call file (VCF)

Git clone and set up mamba environment
```
git clone https://github.com/CherWeiYuan/SpliceAIh.git
cd SpliceAIh

mamba create --name spliceaih pip python=3.7 spliceai=1.3.1 pyvcf=0.6.8 \
pyensembl=2.1.0 pandas tqdm -c bioconda -c anaconda -y

mamba activate spliceaih
pyensembl install --release 108 --species homo_sapiens
pip install .
mamba deactivate
```

Download and unzip GRCh38 genome fasta
```
mkdir -p genome
curl -o genome/hg38.fa.gz \
ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip genome/hg38.fa.gz
```


## Test codes
```
# Unzip test data
tar -xf test_data.tar.xz

# Run test script
mamba activate spliceaih
python3 test.py
mamba deactivate
```

Successful test result:
```
----------------------------------------------------------------------
Ran 43 tests in 3.677s

OK
```


## Trial run on positive control
Note: The below codes require test_data.tar.xz to be decompressed and genome fasta to be downloaded (see above)

Run SpliceAIh on two CFTR variants. This is a good positive control because modifying both the T(n) tract and TG(n) tract in intron 8 leads to incorrect splicing of exon 9. The extent of mis-splicing is increased with less T in T(n) tract and more TG repeats in the TG(n) tract. 

[[Source 1]](https://doi.org/10.1172/JCI639) [[Source 2]](https://doi.org/10.1093/hmg/8.13.2339) [[Source 3]](https://www.mdpi.com/2075-4418/11/2/168)
```
mamba activate spliceaih
mkdir -p sample_output
spliceaih \
--vcf_file test_data/CTFR_T5_TG13.vcf \
--max_read_length 300 \
--fasta genome/hg38.fa \
--annotation_file annotations/grch38.txt \
--spliceai_dist 5000 \
--single_variants \
--outdir sample_output/TEST_CTFR_T5_TG13
mamba deactivate
```
A sample output can be found in this repository in sample_output/TEST_CTFR_T5_TG13

## Output
Three rows from the sample output looks like this:
| variants  | landmark_pos | delta_scores | highest_score |
| ------------- | ------------- | ------------- | ------------- |
| ('chr7', 117548607, 'T', ['TGTGTG']) | 	117548607	| ['TGTGTG\|CFTR\|0.00\|0.04\|0.00\|0.00\|152\|34\|219\|211'] | 0.04 |
| ('chr7', 117548629, 'TTT', ['T']) | 117548629 |	['T\|CFTR\|0.01\|0.00\|0.00\|0.00\|130\|2908\|-16\|189'] | 0.01 |
| [('chr7', 117548607, 'T', ['TGTGTG'], '0\|1'), ('chr7', 117548629, 'TTT', ['T'], '0\|1')]	| 117548607	| ['TGTGTG\|CFTR\|0.00\|0.21\|0.00\|0.00\|150\|32\|333\|209']	| 0.21 |
| [('chr7', 117548607, 'T', ['TGTGTG'], '0\|1'), ('chr7', 117548629, 'TTT', ['T'], '0\|1')] | 117548629 | ['T\|CFTR\|0.01\|0.22\|0.00\|0.00\|130\|12\|-714\|189'] | 0.22  |

[Column 1] A list of (phased) variants. Each variant is annotated as (chromosome, genomic coordinate, reference allele, alternate alleles, genotype if variant is phased).

[Column 2] Landmark position refers to the position that is specified as wild-type in the reference genome and annotation; it will be considered as mutated during SpliceAI execution.

[Column 3] Delta score refers to the [SpliceAI](https://github.com/Illumina/SpliceAI) output of ```ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL```.

[Column 4] The highest score refers to the highest score amongst the four SpliceAI delta scores.

The first two rows in the table are computed when --single flag is raised, so single variants undergo SpliceAI. The scores are low (0.04 & 0.01). When the variants are considered together as phased (third and fourth row), the score is significantly higher (0.21 & 0.22).
