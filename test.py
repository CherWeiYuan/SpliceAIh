"""
Test for SpliceAIc tool
"""

from copy import deepcopy
from spliceaih.auxiliary import *
import unittest
import os

class TestReadVCF(unittest.TestCase):
    def test_basic_vcf(self):
        """
        Test reading basic VCF
        """
        result = find_variant_blocks(create_vcfdf("test_data/read_vcf_test_basic.vcf"), 
                                     max_read_length = 300)
        self.assertCountEqual(result, 
            [
            [("chr1", 93996983, "C", ["T"], "1|0"), ("chr1", 93997213, "G", ["A"], "1|0")]
            ])
    def test_long_vcf(self):
        """
        Test reading long VCF
        """
        result = find_variant_blocks(create_vcfdf("test_data/read_vcf_test_long.vcf"), 
                                     max_read_length = 300)
        self.maxDiff = None
        self.assertCountEqual(result, 
            [
            [("chr1", 93996983, "C", ["T"], "1|0"), ("chr1", 93997213, "G", ["A"], "1|0"),
             ("chr1", 93997351, "G", ["A"], "1|0"), ("chr1", 93997432, "A", ["AT"], "./."),
             ("chr1", 93997549, "G", ["A"], "1|0")], 
            [("chr1", 93998211, "C", ["T"], "0|1"), ("chr1", 93998212, "T", ["C"], "0|1"), 
             ("chr1", 93998239, "C", ["G"], "0|1")],
            [("chr1", 93998713, "ATTTATTTTAT", ["A"], "0|1"), ("chr1", 93998744, "T", ["G"], "0|1"),
             ("chr1", 93998884, "A", ["C"], "0/0"), ("chr1", 93998891, "A", ["C"], "0/0"),
             ("chr1", 93998896, "A", ["C"], "0/0"), ("chr1", 93998899, "C", ["T"], "0|1"),
             ("chr1", 93998905, "A", ["AT"], "./."), ("chr1", 93999047, "A", ["G"], "1/1"),
             ("chr1", 93999245, "C", ["T"], "0|1"), ("chr1", 93999328, "A", ["G"], "0|1"),
             ("chr1", 93999515, "T", ["C"], "0|1")]
            ])
    def test_edge_vcf(self):
        """
        Test reading edge VCF
        """
        result = find_variant_blocks(create_vcfdf("test_data/read_vcf_test_edge.vcf"), 
                                     max_read_length = 300)
        self.assertCountEqual(result, 
            [
            [("chr1", 0, "GGGGGG", ["G"], "0/1"), ("chr1", 12, "C", ["CTTTTT"], "1|0"),
             ("chr1", 100, "GATAGA", ["G"], "1|0")],
            [("chr3", 93997432, "ACCC", ["A"], "0|1"), ("chr3", 93997549, "G", ["A", "GCAT"], "1|0"),
             ("chr3", 93997600, "T", ["TTT", "TA"], "0|1")],
            ])

class TestSplitVariantBlock(unittest.TestCase):
    def test_basic_vcf_split(self):
        """
        Test reading basic VCF
        """
        result = split_variant_blocks(find_variant_blocks(
            create_vcfdf("test_data/read_vcf_test_basic.vcf"), 
            max_read_length = 300))
        self.assertCountEqual(result, 
            [
            [("chr1", 93996983, "C", ["T"], "1|0"), ("chr1", 93997213, "G", ["A"], "1|0")]
            ])
    def test_long_vcf_split(self):
        """
        Test reading long VCF
        """
        result = split_variant_blocks(find_variant_blocks(
            create_vcfdf("test_data/read_vcf_test_long.vcf"), 
            max_read_length = 300))
        self.maxDiff = None
        self.assertCountEqual(result, 
            [
            [("chr1", 93996983, "C", ["T"], "1|0"), ("chr1", 93997213, "G", ["A"], "1|0"),
             ("chr1", 93997351, "G", ["A"], "1|0"), ("chr1", 93997549, "G", ["A"], "1|0")], 
            [("chr1", 93998211, "C", ["T"], "0|1"), ("chr1", 93998212, "T", ["C"], "0|1"), 
             ("chr1", 93998239, "C", ["G"], "0|1")],
            [("chr1", 93998713, "ATTTATTTTAT", ["A"], "0|1"), ("chr1", 93998744, "T", ["G"], "0|1"),
             ("chr1", 93998899, "C", ["T"], "0|1"), ("chr1", 93999245, "C", ["T"], "0|1"), 
             ("chr1", 93999328, "A", ["G"], "0|1"), ("chr1", 93999515, "T", ["C"], "0|1")]
            ])
    def test_edge_vcf_split(self):
        """
        Test reading edge VCF
        """
        result = split_variant_blocks(find_variant_blocks(
            create_vcfdf("test_data/read_vcf_test_edge.vcf"), 
            max_read_length = 300))
        self.assertCountEqual(result, 
            [
            [("chr1", 12, "C", ["CTTTTT"], "1|0"), ("chr1", 100, "GATAGA", ["G"], "1|0")],   
            [("chr3", 93997432, "ACCC", ["A"], "0|1"), ("chr3", 93997549, "G", ["GCAT"], "1|0"),
             ("chr3", 93997600, "T", ["TTT"], "0|1")],
            [("chr3", 93997549, "G", ["A"], "1|0"), ("chr3", 93997600, "T", ["TA"], "0|1")]
            ])

class TestSNP(unittest.TestCase):
    def test_snp_mutation_1(self):
        """
        Test simple SNP mutation
        """
        result = snp_mutation("AAAATAAAAA", 5, "T", "A")
        self.assertEqual(result, "AAAAAAAAAA")

    def test_snp_mutation_2(self):
        """
        Test simple SNP mutation
        """
        result = snp_mutation("GAAAAAAAAA", 1, "G", "A")
        self.assertEqual(result, "AAAAAAAAAA")
        
    def test_snp_mutation_3(self):
        """
        Test simple SNP mutation
        """
        result = snp_mutation("AAAAAAAAAC", 10, "C", "A")
        self.assertEqual(result, "AAAAAAAAAA")

    def test_snp_sequential_mutation(self):
        """
        Test sequential SNP mutation
        """
        result = snp_mutation("AAAAAAAAAC", 10, "C", "A")
        result = snp_mutation(result, 10, "A", "T")
        result = snp_mutation(result, 1, "A", "C")
        self.assertEqual(result, "CAAAAAAAAT")

    def test_snp_zero_index_error(self):
        """
        Test if zero-index error is detected
        """
        self.assertRaises(ZeroIndexError, 
                          snp_mutation, "AAAAAAAAAC", 0, "C", "A")

    def test_snp_ref_error(self):
        """
        Test if wrong reference allele is detected
        """
        self.assertRaises(UnexpectedRefError, 
                          snp_mutation, "AAAAAAAAAA", 1, "C", "A")

class TestInsertion(unittest.TestCase):
    def test_insertion_mutation_1(self):
        """
        Test simple insertion
        """
        result = insertion_mutation("AAAAA", 5, "A", "ATTTTT")
        self.assertEqual(result, "AAAAATTTTT")

    def test_insertion_mutation_2(self):
        """
        Test simple insertion
        """
        result = insertion_mutation("AAAAA", 1, "A", "ATTTTT")
        self.assertEqual(result, "ATTTTTAAAA")

    def test_insertion_mutation_3(self):
        """
        Test simple insertion
        """
        result = insertion_mutation("AAAAA", 3, "A", "ATTTTT")
        self.assertEqual(result, "AAATTTTTAA")

    def test_insertion_sequential_mutation(self):
        """
        Test sequential insertion
        """
        result = insertion_mutation("AAAAA", 1, "A", "ATTTTT")
        result = insertion_mutation(result, 10, "A", "AGG")
        result = insertion_mutation(result, 12, "G", "GCCAT")
        result = insertion_mutation(result, 6, "T", "TATCG")
        self.assertEqual(result, "ATTTTTATCGAAAAGGCCAT")

    def test_insertion_zero_index_error(self):
        """
        Test if zero-index error is detected
        """
        self.assertRaises(ZeroIndexError, 
                          insertion_mutation, "AAAAAAAAAC", 0, "C", "CA")

    def test_insertion_ref_error(self):
        """
        Test if wrong reference allele is detected
        """
        self.assertRaises(UnexpectedRefError, 
                          insertion_mutation, "AAAAAAAAAC", 1, "C", "CA")

class TestDeletion(unittest.TestCase):
    def test_deletion_mutation_1(self):
        """
        Test simple deletion
        """
        result = deletion_mutation("AAATT", 1, "AAA", "A")
        self.assertEqual(result, "ATT")

    def test_deletion_mutation_2(self):
        """
        Test simple deletion
        """
        result = deletion_mutation("AAAAT", 4, "AT", "A")
        self.assertEqual(result, "AAAA")

    def test_deletion_mutation_3(self):
        """
        Test simple deletion
        """
        result = deletion_mutation("ATCGATCGATTGAC", 4, "GATCG", "G")
        self.assertEqual(result, "ATCGATTGAC")

    def test_deletion_mutation_4(self):
        """
        Test simple deletion
        """
        result = deletion_mutation("ATCGATCGATTGAC", 1, "ATCGATCGATTGAC", "A")
        self.assertEqual(result, "A")

    def test_deletion_zero_index_error(self):
        """
        Test if zero-index error is detected
        """
        self.assertRaises(ZeroIndexError, deletion_mutation, 
                          "ATCGATCGATTGAC", 0, "ATCGATCGATTGAC", "A")

    def test_deletion_ref_error(self):
        """
        Test if wrong reference allele is detected
        """
        self.assertRaises(UnexpectedRefError, deletion_mutation, 
                          "ATCGATCGATTGAC", 5, "CGAT", "C")

class TestGenomeMutation(unittest.TestCase):
    def test_two_snp(self):
        """
        Test two simple SNP mutation
        """
        genome_dict = {"chr1" : "AAAAA"}
        result = mutate_genome_and_annotation(
            genome_dict, landmark_pos = 3, chrom = "chr1", 
            pos = [1, 3, 5], ref = ["A", "A", "A"], alt = ["T", "T", "T"], 
            gene_index = 0, annotation_file = "test_data/grch38.txt", 
            TX_END = 1, EXON_START = [1], EXON_END = [2], out = "genome_dict") 
        self.assertEqual(result["chr1"], "TAAAT")

    def test_multiple_snp(self):
        """
        Test four simple SNP mutations
        """
        genome_dict = {"chr1" : "GATAACAAAAAAAAT"}
        result = mutate_genome_and_annotation(
            genome_dict, landmark_pos = 3, chrom = "chr1", pos = [1, 3, 6, 15], 
            ref = ["G", "T", "C", "T" ], alt = ["A", "A", "A", "A"], 
            gene_index = 0, annotation_file = "test_data/grch38.txt", 
            TX_END = 1, EXON_START = [1], EXON_END = [2], out = "genome_dict")
        self.assertEqual(result["chr1"], "AATAAAAAAAAAAAA")

    def test_two_deletions(self):
        """
        Test two simple deletion mutations
        """
        genome_dict = {"chr1" : "AAAAATTTTT"}
        result = mutate_genome_and_annotation(
            genome_dict, landmark_pos = 5, chrom = "chr1", pos = [1, 5, 8], 
            ref = ["AAA", "AT", "TTT"], alt = ["A", "A", "T"], gene_index = 0, 
            annotation_file = "test_data/grch38.txt", TX_END = 10, 
            EXON_START = [1, 4, 8], EXON_END = [2, 5, 10], out = "genome_dict")
        self.assertEqual(result["chr1"], "AAATTT")

    def test_multiple_deletions(self):
        """
        Test five deletion mutations
        """
        genome_dict = {"chr1" : "AATTTCCCCCCCCGGGGAAAAAAAAAAAA"}
        result = mutate_genome_and_annotation(
            genome_dict, landmark_pos = 1, chrom = "chr1", 
            pos = [1, 3, 6, 14, 18], 
            ref = ["AA", "TTT", "CCCCCCCC", "GGGG", "AAAAAAAAAAAA"], 
            alt = ["A", "T", "C", "G", "A"], 
            gene_index = 0, annotation_file = "test_data/grch38.txt", 
            TX_END = 10, EXON_START = [1, 4, 8], EXON_END = [2, 5, 10], 
            out = "genome_dict")
        self.assertEqual(result["chr1"], "AATCGA")

    def test_two_insertions(self):
        """
        Test two insertions
        """
        genome_dict = {"chr1" : "AATTTCCCCCCCCGGGGAAAAAAAAAAAA"}
        result = mutate_genome_and_annotation(
            genome_dict, landmark_pos = 5, chrom = "chr1", pos = [1, 5, 29], 
            ref = ["A", "T", "A"], alt = ["ATTT", "TA", "ATTTGCGGCG"], 
            gene_index = 0, annotation_file = "test_data/grch38.txt", 
            TX_END = 10, EXON_START = [1, 4, 8], EXON_END = [2, 5, 10], 
            out = "genome_dict")
        self.assertEqual(result["chr1"], "ATTTATTTCCCCCCCCGGGGAA" +\
                                         "AAAAAAAAAATTTGCGGCG")

    def test_multiple_insertions(self):
        """
        Test five insertions
        """
        genome_dict = {"chr1" : "TATAACAAGAAAAAG"}
        result = mutate_genome_and_annotation(
            genome_dict, landmark_pos = 6, chrom = "chr1", 
            pos = [1, 3, 6, 9, 15], 
            ref = ["T", "T", "C", "G", "G"], 
            alt = ["TTTT", "TAAA", "CCC", "GGG", "GGGGG"], 
            gene_index = 0, annotation_file = "test_data/grch38.txt", 
            TX_END = 10, EXON_START = [1, 4, 8], EXON_END = [2, 5, 10], 
            out = "genome_dict")
        self.assertEqual(result["chr1"], "TTTTATAAAAACAAGGGAAAAAGGGGG")
        
    def test_mutation_mixture(self):
        """
        Text a mixture of SNP, insertion and deletion mutations
        """
        genome_dict = {"chr1" : "AAAAATTTTTCCCCCGGGGG"}
        result = mutate_genome_and_annotation(
            genome_dict, landmark_pos = 6, chrom = "chr1", 
            pos = [1, 6, 11, 16, 20], 
            ref = ["AAAAA", "T", "C", "G", "G"], 
            alt = ["A", "G", "CCC", "GGGATAT", "GGGGG"], 
            gene_index = 0, annotation_file = "test_data/grch38.txt", 
            TX_END = 10, EXON_START = [1, 4, 8], EXON_END = [2, 5, 10], 
            out = "genome_dict")
        self.assertEqual(result["chr1"], "ATTTTTCCCCCCCGGGATATGGGGGGGG")
                               
    def test_ambiguous_deletion(self):
        """
        Test long deletion that deletes a base that is to be mutated in a later
        position
        """
        genome_dict = {"chr1" : "TATTCTATCTATGGTATGTA"}
        self.assertRaises(AmbiguousDeletionError, mutate_genome_and_annotation,
                          genome_dict, landmark_pos = 6, chrom = "chr1", 
                          pos = [1, 3, 6, 9, 20], 
                          ref = ["T", "TTCTA", "T", "CTATGGTAT", "A"], 
                          alt = ["TTTT", "T", "A", "C", "ATTT"], 
                          gene_index = 0, 
                          annotation_file = "test_data/grch38.txt", 
                          TX_END = 10, 
                          EXON_START = [1, 4, 8], 
                          EXON_END = [2, 5, 10], 
                          out = "genome_dict")
        
    def test_no_expected_change(self):
        """
        Test no expected change when landmark position is the same as the
        only position in pos list
        """
        genome_dict = {"chr1" : "AAAAA"}
        self.assertRaises(NoExpectedChangeError, mutate_genome_and_annotation,
                          genome_dict, landmark_pos = 1, chrom = "chr1", 
                          pos = [1], ref = ["A"], alt = ["ATTTT"], 
                          gene_index = 0, 
                          annotation_file = "test_data/grch38.txt", 
                          TX_END = 10, 
                          EXON_START = [1, 4, 8], 
                          EXON_END = [2, 5, 10], 
                          out = "genome_dict")    
        
    def test_landmark_not_found(self):
        """
        Test if landmark position is found in pos list
        """
        genome_dict = {"chr1" : "AAAAA"}
        self.assertRaises(LandmarkNotFoundError, mutate_genome_and_annotation,
                          genome_dict, 3, "chr1", 
                          [1, 2], ["A", "A"], ["T", "AAAAAA"], 
                          gene_index = 0, 
                          annotation_file = "test_data/grch38.txt", 
                          TX_END = 10, 
                          EXON_START = [1, 4, 8], 
                          EXON_END = [2, 5, 10], 
                          out = "genome_dict")

class TestGeneIndexSearch(unittest.TestCase):
    def test_find_OR4F5(self):
        result = find_gene_index("OR4F5", "test_data/grch38.txt")
        self.assertEqual(result, 0)

    def test_find_PC(self):
        result = find_gene_index("PC", "test_data/grch38.txt")
        self.assertEqual(result, 10548)

    def test_find_CDY1(self):
        result = find_gene_index("CDY1", "test_data/grch38.txt")
        self.assertEqual(result, 19304)

    def test_find_non_existent_gene(self):
        self.assertRaises(GeneNotFoundError, find_gene_index,
                          "TR0LLING", "test_data/grch38.txt")

class TestShiftPosition(unittest.TestCase):
    def test_insertion_shift(self):
        result = shift_annotation_position(current_pos = 5, 
                    shift = 2, mutation_type = "insertion", TX_END = 10, 
                    EXON_START = [1, 3, 5, 8], EXON_END = [2, 4, 6, 10])
        self.assertEqual(result, [12, [1, 3, 7, 10], [2, 4, 8, 12]])

    def test_deletion_shift(self):
        result = shift_annotation_position(current_pos = 5, 
                    shift = 1, mutation_type = "deletion", TX_END = 10, 
                    EXON_START = [1, 3, 5, 8], EXON_END = [2, 4, 6, 10])
        self.assertEqual(result, [9, [1, 3, 4, 7], [2, 4, 5, 9]])

class TestAnnotationEditing(unittest.TestCase):
    def test_SNP(self):
        """
        Test for no position change if all mutations are SNPs
        """
        # Read test data
        genome_dict = parse_genome("test_data/chr7.test.fasta")
        tx_end = 96322093
        exon_start_pos = [96120219,96121654,96121838,96131742,96146555,96170044,96171471,96184276,96184926,96189293,96189580,96191108,96193036,96208837,96234801,96277195,96296897,96321941]
        exon_end_pos = [96121377,96121745,96121997,96131881,96146696,96170125,96171524,96184435,96185011,96189378,96189674,96191247,96193183,96208977,96234917,96277338,96296951,96322093]
        
        # Run mutate_genome() to output a new annotation file: temp_anntest.tsv
        result = mutate_genome_and_annotation(genome_dict, 
                    landmark_pos = 96172431, chrom = "chr7",
                    pos = [96172431, 96172435], ref = ["A", "G"], 
                    alt = ["G", "A"], gene_index = 7245, 
                    annotation_file = "test_data/grch38.txt", 
                    TX_END = tx_end, EXON_START = deepcopy(exon_start_pos), 
                    EXON_END = deepcopy(exon_end_pos), out = "genome_dict")

        # Export and then load new output temp_anntest.tsv
        temp = pd.read_csv("test_data/grch38.txt.temp", sep = "\t", low_memory = False)
        
        # Check results in temp_anntest.tsv
        self.assertEqual(temp.loc[7245, "TX_END"], tx_end)
        result_exon_start = [int(x) for x in temp.loc[7245, "EXON_START"].split(",") 
                            if x != '']
        result_exon_end = [int(x) for x in temp.loc[7245, "EXON_END"].split(",") 
                          if x != ''] 
        self.assertEqual(result_exon_start, exon_start_pos)   
        self.assertEqual(result_exon_end, exon_end_pos)

    def test_deletion(self):
        """
        Test for position change if deletion occurred
        """
        # Read test data
        genome_dict = parse_genome("test_data/chr7.test.fasta")

        tx_end = 96322093
        exon_start_pos = [96120219,96121654,96121838,96131742,96146555,96170044,96171471,96184276,96184926,96189293,96189580,96191108,96193036,96208837,96234801,96277195,96296897,96321941]
        exon_end_pos = [96121377,96121745,96121997,96131881,96146696,96170125,96171524,96184435,96185011,96189378,96189674,96191247,96193183,96208977,96234917,96277338,96296951,96322093]

        expected_tx_end = 96322093 - 4
        expected_exon_start_pos = [x for x in exon_start_pos if x < 96172435] +\
                                [x - 4 for x in exon_start_pos if x >= 96172435]
        expected_exon_end_pos = [x for x in exon_end_pos if x < 96172435] +\
                                [x - 4 for x in exon_end_pos if x >= 96172435]

        # Run mutate_genome() to output a new annotation file: temp_mixtest.tsv
        result = mutate_genome_and_annotation(
            genome_dict, 
            landmark_pos = 96172431, 
            chrom = "chr7", pos = [96172431, 96172435], ref = ["A", "GTAGT"], 
            alt = ["G", "G"], gene_index = 7245, 
            annotation_file = "test_data/grch38.txt", TX_END = tx_end, 
            EXON_START = deepcopy(exon_start_pos), 
            EXON_END = deepcopy(exon_end_pos), out = "genome_dict")

        # Export and then load new output temp_mixtest.tsv
        temp = pd.read_csv("test_data/grch38.txt.temp", sep = "\t", low_memory = False)
        
        # Check results in temp_mixtest.tsv
        result_exon_start = [int(x) for x in temp.loc[7245, "EXON_START"].split(",") 
                             if x != '']
        result_exon_end = [int(x) for x in temp.loc[7245, "EXON_END"].split(",") 
                           if x != ''] 
        self.assertEqual(temp.loc[7245, "TX_END"], expected_tx_end)
        self.assertEqual(result_exon_start, expected_exon_start_pos)   
        self.assertEqual(result_exon_end, expected_exon_end_pos)

    def test_mixed(self):
        """
        Test for position change if SNP, insertion and deletion occurred
        """
        # Read test data
        genome_dict = parse_genome("test_data/chr7.test.fasta")

        tx_end = 96322093
        exon_start_pos = [96120219,96121654,96121838,96131742,96146555,96170044,96171471,96184276,96184926,96189293,96189580,96191108,96193036,96208837,96234801,96277195,96296897,96321941]
        exon_end_pos = [96121377,96121745,96121997,96131881,96146696,96170125,96171524,96184435,96185011,96189378,96189674,96191247,96193183,96208977,96234917,96277338,96296951,96322093]

        expected_tx_end = 96322093 - 3
        expected_exon_start_pos = [x for x in exon_start_pos if x < 96172435] +\
                                [x - 3 for x in exon_start_pos if x >= 96172435]
        expected_exon_end_pos = [x for x in exon_end_pos if x < 96172435] +\
                                [x - 3 for x in exon_end_pos if x >= 96172435]

        # Run mutate_genome() to output a new annotation file: temp_mixtest.tsv
        result = mutate_genome_and_annotation(genome_dict, 
                    landmark_pos = 96172431, chrom = "chr7", 
                    pos = [96172431, 96172435, 96172450, 96172451], 
                    ref = ["A", "GTAGT", "T", "G"], alt = ["G", "G", "TT", "A"], 
                    gene_index = 7245, annotation_file = "test_data/grch38.txt", 
                    TX_END = tx_end, EXON_START = deepcopy(exon_start_pos), 
                    EXON_END = deepcopy(exon_end_pos), out = "genome_dict")

        # Export and then load new output temp_mixtest.tsv
        temp = pd.read_csv("test_data/grch38.txt.temp", sep = "\t", low_memory = False)
        
        # Check results in temp_mixtest.tsv
        result_exon_start = [int(x) for x in temp.loc[7245, "EXON_START"].split(",") 
                             if x != '']
        result_exon_end = [int(x) for x in temp.loc[7245, "EXON_END"].split(",") 
                           if x != ''] 
        self.assertEqual(temp.loc[7245, "TX_END"], expected_tx_end)
        self.assertEqual(result_exon_start, expected_exon_start_pos)   
        self.assertEqual(result_exon_end, expected_exon_end_pos)

if __name__ == "__main__":
    unittest.main()
    try:
        os.remove("temp_test.tsv")
    except:
        pass


