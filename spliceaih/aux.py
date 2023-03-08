from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from copy import deepcopy
import logging
import os
import pandas as pd
import vcf
from spliceaih.utils import Annotator, get_delta_scores
import sys

def init_logging(out):
    """
    Initialise the logging facility. Write log statement indicating the 
    program has started, and also write out the command line from sys.argv
    and allow other logging statements throughout the program to be put in 
    the same file
    """
    logging.root.handlers = []
    logging.basicConfig(level = logging.INFO,
                        format = "%(asctime)s %(levelname)s - %(message)s",
                        datefmt = "%Y-%m-%dT%H:%M:%S%z",
                        handlers = [logging.FileHandler(out + "/log.txt", "w+"),
                                    logging.StreamHandler(sys.stdout)])
    logging.info("Program started")
    logging.info("Command line: %s", " ".join(sys.argv))

def create_vcfdf(vcf_file):
    """
    Create Pandas dataframe from VCF file
    
    Parameters
    ----------
    vcf_file: TYPE str. VCF file name
    
    Returns
    -------
    vcf_df: TYPE Pandas dataframe. Rows sorted by chrom then pos.
            Four columns | chrom (str), pos (int), ref (str), alt (list of str)

    Notes
    ------
    Cells in alt column contain a list of alternate alleles
    """
    # Open VCF file with pyvcf
    vcf_reader = vcf.Reader(open(vcf_file, "r"))
    sample_names = vcf_reader.samples
    if len(sample_names) > 1:
        sys.exit("Error: More than one sample names found in input VCF")
    else:
        sample_name = sample_names[0]
    # Loop through VCF and create dataframe
    vcfdf_seed = []
    for record in vcf_reader:
        chrom = str(record.CHROM)
        pos = int(record.POS) # VCF is 1-based; keep it that way
        ref = str(record.REF)
        alt = [str(x) for x in list(record.ALT) if x != "*"]
        genotype = str(record.genotype(sample_name).data.GT)
        vcfdf_seed.append([chrom, pos, ref, alt, genotype])
    vcf_df = pd.DataFrame(vcfdf_seed, columns = ["chrom", "pos", "ref", "alt", 
                                                 "genotype"])
    return vcf_df

def find_variant_blocks(vcf_df, max_read_length):
    """
    Find blocks of close variants (within spliceai_dist)

    Parameters
    ----------
    vcf_df: TYPE Pandas dataframe. Rows sorted by chrom then pos.
            Four columns | chrom (str), pos (int), ref (str), alt (list of str)
            
    read_length: TYPE int. Sequencing read length
    
    Returns
    -------
    blocks_list: TYPE list. List of variants located within read_length
                 For example, a blocks_list with two blocks is represented as:
    [
    [(chrom, pos, ref, alt, gt1), (chrom, pos, ref, alt, gt1)], # block 1
    [(chrom, pos, ref, alt, gt2), (chrom, pos, ref, alt, gt3), 
     (chrom, pos, ref, alt, gt2), # block 2
    ]
    
    Notes
    ------
    Sorted vcf_df is critical for correct function
    """
    # Sort vcf_df because order by chromosome and then position is required
    # to determine if neighboring variants should be placed in blocks
    vcf_df.sort_values(by = ["chrom", "pos"], inplace = True, ascending = True,
                      ignore_index = True)
    
    # Loop in sorted vcf_df to find variant blocks
    prev_chrom = None
    prev_pos = None
    blocks_list = []
    current_block = []
    for index, row in vcf_df.iterrows():
        chrom = vcf_df.loc[index, "chrom"]
        pos = vcf_df.loc[index, "pos"]
        ref = vcf_df.loc[index, "ref"]
        alt = vcf_df.at[index, "alt"]
        genotype = vcf_df.loc[index, "genotype"]
        # If within distance, add to current block
        if chrom == prev_chrom and abs(pos - prev_pos) <= max_read_length:
            current_block.append((chrom, pos, ref, alt, genotype))
        # If not within distance, flush current block (if not empty) into 
        # block list
        else:
            if len(current_block) == 1: 
                current_block = []
            elif len(current_block) > 1:
                blocks_list.append(current_block)
                current_block = []
            current_block.append((chrom, pos, ref, alt, genotype))
        prev_chrom = chrom
        prev_pos = pos
        
        # If last row is reached
        if index == len(vcf_df) - 1 and len(current_block) > 1:
            blocks_list.append(current_block)
    return blocks_list

def split_variant_blocks(blocks_list):
    """
    Split blocks of close variants according to their haplotypes

    Parameters
    ----------
    blocks_list: TYPE list. List of variants located within read_length
                 For example, a blocks_list with two blocks is represented as:
    [
    [(chrom, pos, ref, alt, gt1), (chrom, pos, ref, alt, gt1)], # block 1
    [(chrom, pos, ref, alt, gt1), (chrom, pos, ref, alt, gt2), 
     (chrom, pos, ref, alt, gt1)], # block 2
    ]

    Returns
    -------
    blocks_list: TYPE list. List of haplotypes, such as:
    [
    [(chrom, pos, ref, alt, gt1), (chrom, pos, ref, alt, gt1)], # block 1
    [(chrom, pos, ref, alt, gt2), (chrom, pos, ref, alt, gt3)],
    [(chrom, pos, ref, alt, gt2)]
     ]
    """
    haplotype_list = []
    for block in blocks_list:
        if len(block) == 1:
            continue
        block_haplotype = {}
        for variant in block:
            gt = str(variant[4])
            # If variant is not phased, ignore
            if "|" not in gt:
                continue
            # If variant is phased, add to dictionary block_haplotype
            else:
                chrom = variant[0]
                pos = variant[1]
                ref = variant[2]
                alt = variant[3]           # List of alternate alleles
                gt_split = variant[4].split("|") # List of haplotypes corresponding to alt
                if len(alt) != 1 and len(alt) != len(gt_split):
                    logging.warning("WARNING: Number of alt alleles do not " +\
                                    "correspond to number of haplotypes for " +\
                                    f"variant {variant}")
                for i in range(len(alt)):
                    try:
                        block_haplotype[gt_split[i]] += [(chrom, pos, ref, 
                                                          [alt[i]], gt)]
                    except KeyError:
                        block_haplotype[gt_split[i]] = [(chrom, pos, ref, 
                                                         [alt[i]], gt)]
        for key in block_haplotype.keys():
            haplotype_list.append(block_haplotype[key])
    return haplotype_list

def update_df_seed_single_variant(vcf_df, output_seed, genome_fasta, 
                                  annotation_file, ensembl, spliceai_dist):
    for index, row in vcf_df.iterrows():
        single_chrom = vcf_df.loc[index, "chrom"]
        single_pos = int(vcf_df.loc[index, "pos"])
        single_ref = vcf_df.loc[index, "ref"]
        single_alt = vcf_df.loc[index, "alt"]
        genes_nearby = ensembl.gene_names_at_locus(
            contig = single_chrom.replace("chr", ""), position = single_pos)
        for gene_name in genes_nearby:
            if gene_name == "":
                continue
            try:
                gene_index = find_gene_index(gene_name, annotation_file)
                exon_pos_dict = get_exon_positions(gene_name, annotation_file, 
                                                   gene_index)
            except GeneNotFoundError:
                logging.warning(f"Warning: {gene_name} is not found in " +\
                                "annotation file")
            delta_scores = get_delta_scores(vcfRecord(single_chrom,   # CHROM
                                                      single_pos,     # POS
                                                      single_ref,     # REF
                                                      single_alt),    # ALT
                                            Annotator(genome_fasta, 
                                                      annotation_file), 
                                                      spliceai_dist, mask = 1)
            highest_score = parse_spliceai(delta_scores, gene_name)
            output_seed.append(((single_chrom, single_pos, single_ref, single_alt), 
                                single_pos, delta_scores, highest_score))
    return output_seed

def update_df_seed_haplotype(haplotype_block, output_seed, 
                             chrom, pos_list, ref_list, alt_list, 
                             annotation_file, genome_dict, ensembl, 
                             spliceai_dist):
    """
    For each haplotype, run SpliceAI on each landmark position
    and update dataframe seed. The dataframe seed will be used to create an
    output dataframe at the end of the programme
    """
    for k in range(len(pos_list)):
        landmark_pos = int(pos_list[k])

        # Find genes nearby
        genes_nearby = ensembl.gene_names_at_locus(
            contig = chrom.replace("chr", ""), position = landmark_pos)
        for gene_name in genes_nearby:
            pos_list_copy = deepcopy(pos_list)
            if gene_name == "":
                continue
            try:
                gene_index = find_gene_index(gene_name, annotation_file)
                exon_pos_dict = get_exon_positions(gene_name, annotation_file, 
                                                   gene_index)
            except GeneNotFoundError:
                logging.warning(f"Warning: {gene_name} is not found in " +\
                                 "annotation file")
            TX_END = exon_pos_dict["TX_END"]
            EXON_START = exon_pos_dict["EXON_START"]
            EXON_END = exon_pos_dict["EXON_END"]

            # Positions are fed into function as 1-based
            try:
                new_landmark_pos = mutate_genome_and_annotation(
                    {chrom: genome_dict[chrom]}, # Careful not to modify genome_dict
                    landmark_pos, chrom, pos_list_copy, 
                    ref_list, alt_list, gene_index, annotation_file, 
                    TX_END, EXON_START, EXON_END)
            except AmbiguousDeletionError: 
                continue
            print(f"NEW LANDMARK: {new_landmark_pos}")
            delta_scores = get_delta_scores(
                vcfRecord(chrom,              # CHROM
                          new_landmark_pos,   # landmark POS after shift
                          ref_list[k],        # REF
                          [alt_list[k]]),     # ALT
                Annotator("spliceaih_temp/temp_genome.fasta", 
                          f"{annotation_file}.temp"), spliceai_dist, mask = 1)
            highest_score = parse_spliceai(delta_scores, gene_name)
            output_seed.append((haplotype_block, landmark_pos, delta_scores, 
                                highest_score))
            for item in ("spliceaih_temp/temp_genome.fasta", 
                         "spliceaih_temp/temp_genome.fasta.fai",
                        f"spliceaih_temp/{annotation_file}.temp"):
                try:
                    os.remove(item)
                except FileNotFoundError:
                    pass
    return output_seed

def parse_genome(genome_fasta):
    """
    Input
        genome_fasta: directory and name of fasta file
    Output
        genome_dict: dictionary of key (chromosome) and value (chromosomal 
                     sequence)
    Note 
        chrom.description is used instead of chrom.id to get full fasta header
        Otherwise, spaces in fasta header will lead to its truncation
    """
    ref_genome = SeqIO.parse(genome_fasta, "fasta")
    genome_dict = {}
    for contig in ref_genome:
        genome_dict[str(contig.name)] = str(contig.seq).upper()
    if genome_dict:
        pass
    else:
        raise ChromNotFoundError()
    return genome_dict

def snp_mutation(seq, pos, ref, alt):
    if pos == 0:
        raise ZeroIndexError()
    actual_ref = seq[pos-1]
    if ref == actual_ref:
        seq = seq[0: pos-1] + alt + seq[pos:]    
    else:
        raise UnexpectedRefError(f"Reference allele at position {pos + 1} is " +\
                                  "wrongly specified. Please check reference " +\
                                  "chromosome position.")
    return seq

def insertion_mutation(seq, pos, ref, alt):
    if pos == 0:
        raise ZeroIndexError()
    actual_ref = seq[pos-1]
    if ref == actual_ref:
        seq = seq[0: pos-1] + alt + seq[pos:]
    else:
        raise UnexpectedRefError(f"Reference allele at position {pos + 1} is " +\
                                  "wrongly specified. Please check reference " +\
                                  "chromosome position.")
    return seq

def deletion_mutation(seq, pos, ref, alt):
    if pos == 0:
        raise ZeroIndexError()
    num_del = len(ref) - len(alt)
    actual_ref = seq[pos-1: pos+num_del]
    if ref == actual_ref:
        seq = seq[0: pos] + seq[pos+num_del: ]
    elif alt == "": 
        seq = seq[0: pos-1] + seq[pos+num_del: ]
    else:
        raise UnexpectedRefError(f"Reference allele at position {pos + 1} is " +\
                                 "wrongly specified. Please check reference " +\
                                 "chromosome position.")
    return seq

def find_gene_index(gene_name, annotation_file):
    """
    Find index of gene in annotation dataframe
    """
    df = pd.read_csv(annotation_file, sep = "\t", header = 0,
                     names = ["NAME", "CHROM", "STRAND", "TX_START", "TX_END" ,
                              "EXON_START", "EXON_END"])
    gene_index = df.NAME[df.NAME == gene_name].index.to_list()
    if len(gene_index) == 0:
        raise GeneNotFoundError("Please ensure gene name can be found in " +\
                                "annotation file")
    elif len(gene_index) > 1:
        raise DuplicateAnnotationEntries("Please ensure annotation file does " +\
                                         "not have duplicate gene entries")
    elif len(gene_index) == 1:
        gene_index = int(gene_index[0])

    return gene_index

def process_annotation_file(annotation_file):
    """
    Create a temp annotation file for editing
    """
    try:
        os.remove(f"{annotation_file}.temp")
    except FileNotFoundError:
        pass
    temp = pd.read_csv(annotation_file, sep = "\t")
    temp.to_csv(f"{annotation_file}.temp", sep = "\t", index = False)
    pass

def get_exon_positions(gene_name, annotation_file, gene_index):
    """
    Parameters
    ----------
    gene_name : TYPE
        DESCRIPTION.
    annotation_file : TYPE
        DESCRIPTION.

    Returns
    -------
    out_dict : TYPE
        DESCRIPTION.
    """
    temp = pd.read_csv(f"{annotation_file}", sep = "\t")
    TX_END = temp.loc[gene_index, "TX_END"]
    EXON_START = [int(x) for x in temp.loc[gene_index, "EXON_START"].split(",") 
                  if x != '']
    EXON_END = [int(x) for x in temp.loc[gene_index, "EXON_END"].split(",") 
                if x != '']
    exon_pos_dict = {"TX_END": TX_END, 
                     "EXON_START": EXON_START, 
                     "EXON_END": EXON_END}
    return exon_pos_dict

def shift_annotation_position(current_pos, shift, mutation_type,
                              TX_END, EXON_START, EXON_END):
    """
    Change chromosomal positions of annotation due to insertion or deletion 
    mutations
    """
    if mutation_type == "insertion":
        shift = abs(shift)
    elif mutation_type == "deletion":
        shift = shift * -1
    else:
        raise UnexpectedMutationTypeError("Mutation type is neither " +\
                                          "'insertion' or 'deletion'")
    TX_END += shift
    
    for i in range(len(EXON_START)):
        if EXON_START[i] >= current_pos:
            EXON_START[i] += shift
    
    for i in range(len(EXON_END)):
        if EXON_END[i] >= current_pos:
            EXON_END[i] += shift 

    return [TX_END, EXON_START, EXON_END]     

def write_fasta(genome_dict, out):
    """
    Write genome dictionary as fasta file
    """
    record_list = []
    for key, value in genome_dict.items():
        record = SeqRecord(Seq(value), id = key, description = "")
        record_list.append(record)
    with open(out, "w") as output_handle:
        SeqIO.write(record_list, output_handle, "fasta")

def mutate_genome_and_annotation(genome_dict, landmark_pos, chrom, pos, ref, 
                                 alt, gene_index, annotation_file, TX_END, 
                                 EXON_START, EXON_END):
    """
    Returns mutated genome dictionary
    
    Ignore position at landmark because this position will be considered a
    variant position in SpliceAI and thereby needs to stay WT in ref. genome

    Positions in pos of input are 1-based
    """
    # Create a temp annotation file for editing
    process_annotation_file(annotation_file)

    # Raise error if index is zero in a 1-index system
    if 0 in pos:
        raise ZeroIndexError()
    
    # If there is only one position in pos, no mutation will be performed as
    # one variant position must retain WT status for SpliceAI
    if len(pos) == 1:
        raise NoExpectedChangeError()
    
    # If deletions delete a base that is supposed to be mutated,
    # move on to next variant block and leave a logging warning
    for i in range(len(pos) - 1):
        length_change = len(ref[i]) - len(alt[i])
        if length_change >= 1:
            if pos[i] + abs(length_change) >= pos[i + 1]:
                logging.info("WARNING: Deletion at position " +\
                             f"{chrom}:{pos[i]} deleted the base to be " +\
                             f"mutated at position {chrom}:{pos[i + 1]}")
                raise AmbiguousDeletionError()

    # Find index of landmark position in pos list
    if landmark_pos in pos:
        pivot = pos.index(landmark_pos)
    else:
        raise LandmarkNotFoundError()
    
    #pos = sorted(pos) # Sorting of pos but not reference and alternate allele 
                       # list will create bugs
    shift_seq = 0
    shift_annotation = 0
    new_landmark_pos = landmark_pos
    for i in range(len(pos)):
        # Adjust landmark pos for insertions and deletions before it
        if i == pivot and landmark_pos >= pos[i]:
            new_landmark_pos += shift_seq
        
        elif i != pivot:
            # Adjust current position by amount of shift done by the 
            # previous iteration
            pos[i] += shift_seq

            ## Modify reference sequence
            if len(ref[i]) == len(alt[i]):
                genome_dict[chrom] = snp_mutation(genome_dict[chrom], 
                                                  pos[i], ref[i], alt[i])
            # Insertion
            elif len(ref[i]) < len(alt[i]):
                genome_dict[chrom] = insertion_mutation(genome_dict[chrom], 
                                                        pos[i], ref[i], alt[i])
            # Deletion
            elif len(ref[i]) > len(alt[i]):
                genome_dict[chrom] = deletion_mutation(genome_dict[chrom], 
                                                       pos[i], ref[i], alt[i])
            else:
                raise UnexpectedMutation()

            shift_seq += len(alt[i]) - len(ref[i])

            ## Modify annotation file
            shift_annotation = abs(len(alt[i]) - len(ref[i]))
            
            # Insertion
            if len(ref[i]) < len(alt[i]):
                TX_END, EXON_START, EXON_END = shift_annotation_position(pos[i], 
                                                shift_annotation, "insertion", 
                                                TX_END, EXON_START, EXON_END)
            # Deletion
            elif len(ref[i]) > len(alt[i]):
                TX_END, EXON_START, EXON_END = shift_annotation_position(pos[i], 
                                                shift_annotation, "deletion", 
                                                TX_END, EXON_START, EXON_END)
            

    # Adjust TX_END, EXON_START, EXON_END in annotation file
    df = pd.read_csv(f"{annotation_file}.temp", sep = "\t", 
                     low_memory = False)
    
    exon_start = ''
    for i in EXON_START:
        exon_start += f"{i},"
    exon_end = ''
    for i in EXON_END:
        exon_end += f"{i},"

    df.loc[gene_index, "TX_END"] = TX_END
    df.loc[gene_index, "EXON_START"] = exon_start
    df.loc[gene_index, "EXON_END"] = exon_end

    # Clear previous files to prevent failure to overwrite
    for item in ("spliceaih_temp/temp_genome.fasta", 
                 "spliceaih_temp/temp_genome.fasta.fai"):
        try:
            os.remove(item)
        except FileNotFoundError:
            pass

    # Output files
    df.to_csv(f"{annotation_file}.temp", sep = "\t", index = False)
    write_fasta(genome_dict, "spliceaih_temp/temp_genome.fasta")
    return new_landmark_pos

class vcfRecord():
    def __init__(self, CHROM, POS, REF, ALT):
        self.chrom = CHROM
        self.pos = POS
        self.ref = REF
        self.alts = ALT

def str_to_numeric(string_number):
    """
    Convert probabilities from SpliceAI or Pangolin from string to numeric
    """
    # Some SpliceAI entries have "." instead of probabilities
    if string_number == ".":
        return 0
    try:
        return float(string_number)
    except:
        return int(string_number)

def parse_spliceai(delta_scores, gene_name):
    """
    Parse delta scores from SpliceAI

    Parameters
    ----------
    delta_scores : TYPE list
        DESCRIPTION. Formatted as a list of strings:
            ["ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL",
             "ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL"...]
            
    gene_name : TYPE str
        DESCRIPTION. Name of gene. Gene name should be provided if SpliceAI
                     is run with modified annotation because the annotation is
                     edited with one gene in mind.
    Returns
    -------
    out : TYPE int
        DESCRIPTION. Highest score in the score string from SpliceAI
    """
    out = []
    if len(delta_scores) == 0:
        return 0
    for score_string in delta_scores:
        record = score_string.split("|")
        if gene_name != None and record[1] != gene_name:
            out.append(0)
        out.append(max([str_to_numeric(record[2]), str_to_numeric(record[3]),   
                        str_to_numeric(record[4]), str_to_numeric(record[5])]))
    return max(out)

def check_input(chrom, pos_list, ref_list, alt_list):
    if len(pos_list) == len(ref_list) == len(alt_list):
        pass
    else:
        InputCheckError("Input consistency check: FAILED\n" +\
                 "Please ensure the number of items in pos, ref and alt " +\
                 "are the same")
    if type(pos_list) == type(ref_list) == len(alt_list) == list:
        pass
    else:
        InputCheckError("Input type check: FAILED\n" +\
            "Please ensure the input type of pos, ref and alt is a list")
    for lst in [pos_list, ref_list, alt_list]:
        if any(isinstance(i, list) for i in lst):
            InputCheckError("Input nestedness check: FAILED\n" +\
                     "Please ensure there are no nested lists within pos, " +\
                     "ref and alt")
    pass

# Error messages
class ChromNotFoundError(Exception):
    def __init__(self):
        super().__init__("User-supplied chromosome name is not found in fasta")    

class ZeroIndexError(Exception):
    def __init__(self):
        super().__init__("ERROR: Input positions is 1-indexed but position 0 was observed")

class UnexpectedRefError(Exception):
    pass

class GeneNotFoundError(Exception):
    pass

class DuplicateAnnotationEntries(Exception):
    pass

class LandmarkNotFoundError(Exception):
    def __init__(self):
        super().__init__("Landmark position is not found in the list of " +\
                         "positions to be mutated")

class NoExpectedChangeError(Exception):
    def __init__(self):
        super().__init__("There is only one position to be mutated and it is" +\
                         "the landmark position that must not be mutated. " +\
                         "Use snp_mutation(), insertion_mutation or " +\
                         "deletion_mutation() instead of mutate_genome()")
    
class AmbiguousDeletionError(Exception):
    pass

class UnexpectedMutation(Exception):
    def __init__(self):
        super().__init__("Mutation is not SNP, insertion or deletion")

class UnexpectedMutationTypeError(Exception):
    pass

class InputCheckError(Exception):
    pass