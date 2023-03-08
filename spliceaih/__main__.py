"""
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) CHER_WEI_YUAN, 07JUL2022 
License     : -
Maintainer  : E0031403@U.NUS.EDU
Portability : POSIX

This program runs SpliceAI on haplotypes provided by the input VCF
"""

from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from copy import deepcopy
import pandas as pd
import pyensembl
from spliceaih.aux import *
from spliceaih.utils import Annotator, get_delta_scores
from tqdm import tqdm

def parse_args():
    """Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    """
    description = "Performs SpliceAI on haplotypes"
    parser = ArgumentParser(description = description, 
                            formatter_class = RawTextHelpFormatter)
    parser.add_argument("-v", "--vcf_file",
                        type = str,
                        default = None,
                        help = "Name of variant call file")
    parser.add_argument("-l", "--max_read_length",
                        type = int,
                        default = 150,
                        help = "Max sequencing read length")
    parser.add_argument("-f", "--fasta",
                        type = str,
                        default = None,
                        required = True,
                        help = "Name of genome fasta file. " +\
                               "Use either GRCh37 or GRCh38")
    parser.add_argument("-a", "--annotation_file",
                        type = str,
                        required = True,
                        help = "grch38.txt or grch37.txt file downloaded " +\
                               " from https://github.com/Illumina/SpliceAI/tree/master/spliceai/annotations")   
    parser.add_argument("-d", "--spliceai_dist",
                        type = int,
                        default = 5000,
                        help = "Maximum distance between the variant and gained/lost splice site") 
    parser.add_argument("-s", "--single_variants", 
                        action = "store_true",
                        help = "Run SpliceAI on both single variants and haplotypes")
    parser.add_argument("-o", "--outdir",
                        type = str,
                        default = "Sample",
                        help = "Output directory")
    return parser.parse_args()

def main():
    args = parse_args()
    vcf_file = args.vcf_file
    genome_fasta = args.fasta
    max_read_length = args.max_read_length
    annotation_file = args.annotation_file
    spliceai_dist = args.spliceai_dist
    single_variants = args.single_variants
    outdir = args.outdir
    
    # Make output directory
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    # Make temp directory
    isExist = os.path.exists("spliceaih_temp")
    if not isExist:
        os.makedirs("spliceaih_temp")

    # Log
    init_logging(outdir)
    
    # Get haplotypes
    vcf_df = create_vcfdf(vcf_file)
    variant_blocks = find_variant_blocks(vcf_df, max_read_length)
    haplotype_blocks = split_variant_blocks(variant_blocks)

    # Free memory
    variant_blocks = None

    # Load genomic data
    genome_dict = parse_genome(genome_fasta)
    ensembl = pyensembl.EnsemblRelease(108)

    # Create output dataframe seed
    output_seed = []

    # Update single variant SpliceAI results in output dataframe seed
    logging.info("Running SpliceAI on single variants")
    if single_variants:
        output_seed = update_df_seed_single_variant(vcf_df, output_seed, 
                      genome_fasta, annotation_file, ensembl, spliceai_dist)

    # Update haplotype SpliceAI results in output dataframe seed
    logging.info("Running SpliceAI on haplotypes")
    pbar = tqdm(total = len(haplotype_blocks))
    for haplotype_block in haplotype_blocks:
        # If there is only one variant in the block, there is no phasing
        if len(haplotype_block) <= 1:
            pbar.update()
            continue
        chrom = []
        pos_list = [] # Positions are extracted as 1-based
        ref_list = []
        alt_list = []
        for item in haplotype_block:
            chrom.append(item[0])
            pos_list.append(item[1]) # Position extracted as 1-based
            ref_list.append(item[2])
            alt_list.append(item[3][0])
        chrom = list(set(chrom))
        check_input(chrom, pos_list, ref_list, alt_list)
        if len(chrom) == 1:
            chrom = chrom[0]
        else:
            sys.exit("Error: There are multiple chroms in one " +\
                     "variant block")
        output_seed = update_df_seed_haplotype(haplotype_block, output_seed, 
                        chrom, deepcopy(pos_list), ref_list, alt_list, 
                        annotation_file, genome_dict, ensembl, spliceai_dist)
        pbar.update()

    # Output as TSV
    out_df = pd.DataFrame(output_seed, columns = ["variants", 
                                                  "landmark_pos", 
                                                  "delta_scores", 
                                                  "highest_score"])
    out_df.to_csv(f"{outdir}/spliceaih_out.tsv", sep = "\t", index = False)

    # Clear previous files to prevent failure to overwrite
    for item in ("spliceaih_temp/temp_genome.fasta", 
                 "spliceaih_temp/temp_genome.fasta.fai",
                 f"spliceaih_temp/{annotation_file}.temp"):
        try:
            os.remove(item)
        except FileNotFoundError:
            pass

    logging.info("Complete")

if __name__ == "__main__":
    main()