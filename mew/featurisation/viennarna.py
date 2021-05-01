#!/usr/bin/env python

import RNA

BCD = "ATGGGCCCAAGTTCACTTAAAAAGGAGATCAACAATGAAAGCAATTTTCGTACTGAAACATCTTAATCATGCAGGGGAGGGTTTCTA"


def get_rna_coordinate(rna_coordinate, bcd):
    coordinate = rna_coordinate - len(bcd)
    if coordinate <= 0:
        coordinate -= 1

    return coordinate


def from_bppm_base_totals(sequence, bcd=''):
    """Return feature vector and mapping, with each feature the total probability of a base's involvement in pairing.

    :param sequence: str, DNA/RNA sequence of gene
    :param bcd: str, DNA/RNA sequence of upstream sequence

    :return feature_vector: list of float, with each float the total probability of a base's involvement in pairing
    :return feature_mapping: dict of {index: base position, ->}, with index and base position int. Index corresponds to
                             index in the feature vector
    """
    bppm = BPPM(sequence, bcd)
    feature_vector, feature_mapping = bppm.make_row_feature_vector()
    return feature_vector, feature_mapping


def from_bppm(sequence, bcd=''):
    """Return feature vector and mapping, with each feature the probability of two bases pairing.

    :param sequence: str, DNA/RNA sequence of gene
    :param bcd: str, DNA/RNA sequence of upstream sequence

    :return feature_vector: list of float, with each float the probability of two bases pairing
    :return feature_mapping: dict of {index: (base 1, base 2), ->}, with index, base 1 and base 2 int. Base 1 and
                             base 2 are the bases between which the pairing probability at feature_vector[index] was
                             calculated.
    """
    bppm = BPPM(sequence, bcd)
    feature_vector, feature_mapping = bppm.make_feature_vector()
    return feature_vector, feature_mapping


class BPPM:
    def __init__(self, sequence, bcd):
        self.full_sequence = bcd + sequence
        self.sequence = sequence
        self.bcd = bcd

        fc = RNA.fold_compound(self.full_sequence)
        self.propensity, self.ensemble_energy = fc.pf()
        self.bppm = fc.bpp()

    def make_feature_vector(self):
        feature_mapping = {}
        index = 0
        feature_vector = []
        for i in range(1, len(self.full_sequence)):
            coordinate_1 = get_rna_coordinate(i, self.bcd)
            for j in range(i + 1, len(self.full_sequence) + 1):
                coordinate_2 = get_rna_coordinate(j, self.bcd)
                feature_vector.append(self.bppm[i][j])
                feature_mapping[index] = (coordinate_1, coordinate_2)
                index += 1

        return feature_vector, feature_mapping

    def make_row_feature_vector(self):
        feature_mapping = {}
        index = 0
        feature_vector = []
        for i in range(1, len(self.full_sequence) + 1):
            bp_prob = sum(self.bppm[i])
            coordinate = get_rna_coordinate(i, self.bcd)
            feature_vector.append(bp_prob)
            feature_mapping[index] = coordinate
            index += 1

        return feature_vector, feature_mapping

    def print_bppm(self):
        for i in range(1, len(self.full_sequence)):
            for j in range(i + 1, len(self.full_sequence) + 1):
                print("pr(%d,%d) = %g" % (i - len(self.bcd), j - len(self.bcd), self.bppm[i][j]))



if __name__ == "__main__":
    bppm = BPPM("ATGGCATCATCAGAAGACGTAATTAAAGAATTTATGCGTTTTAAAGTTCGAATGGAGGGATCAGTTAACGGTCACGAATTCGAAATTGAAGGTGAGGGTGAAGGTCGTCCTTATGAGGGTACTCAGACGGCTAAGCTTAAGGTTACTAAGGGTGGTCCTCTCCCTTTCGCCTGGGATATTCTATCTCCTCAATTTCAGTATGGTTCTAAGGCCTATGTGAAGCATCCAGCCGATATTCCTGATTATCTGAAACTATCGTTTCCGGAGGGATTCAAGTGGGAGCGCGTGATGAACTTTGAAGATGGAGGCGTTGTTACAGTGACACAAGATTCATCTCTGCAGGATGGTGAGTTTATCTATAAGGTTAAGCTACGTGGGACGAATTTTCCTTCTGATGGTCCGGTTATGCAGAAGAAGACCATGGGTTGGGAGGCATCTACAGAGCGAATGTACCCGGAGGATGGTGCTCTAAAGGGCGAGATCAAGATGCGGCTGAAGCTGAAGGACGGTGGTCATTATGATGCTGAGGTGAAGACGACTTATATGGCTAAGAAGCCGGTGCAGCTGCCTGGTGCGTATAAGACTGACATAAAACTCGACATCACCTCCCACAACGAGGACTACACCATCGTCGAGCAGTACGAACGTGCCGAGGGACGCCACTCCACCGGCGCCTGA", BCD)
    print(bppm.make_row_feature_vector())
    bppm_2 = BPPM("ATGGCATCATCTGAAGACGTCATCAAAGAATTCATGCGCTTTAAAGTCCGCATGGAAGGATCAGTTAACGGGCATGAATTCGAAATCGAGGGGGAGGGCGAAGGTCGTCCCTATGAGGGAACACAAACAGCAAAGCTCAAGGTCACAAAGGGAGGGCCGCTACCATTCGCCTGGGATATCCTCTCCCCGCAGTTCCAGTACGGCTCCAAGGCGTACGTGAAGCACCCCGCGGACATCCCTGATTATCTGAAACTCTCCTTTCCGGAGGGCTTTAAATGGGAGCGAGTGATGAATTTCGAAGATGGGGGGGTAGTAACTGTTACTCAGGACTCGTCACTGCAGGATGGAGAATTTATTTACAAGGTAAAACTTCGGGGAACAAATTTTCCGTCGGATGGTCCAGTGATGCAGAAAAAGACCATGGGATGGGAAGCTTCCACAGAGCGCATGTACCCCGAGGATGGCGCGCTAAAGGGTGAAATTAAAATGCGACTTAAGCTTAAGGACGGCGGCCATTACGATGCCGAAGTAAAAACCACCTACATGGCCAAGAAACCTGTACAACTTCCCGGTGCTTATAAAACGGACATAAAGCTTGATATCACGTCGCATAATGAGGATTATACTATTGTGGAGCAGTATGAACGGGCCGAGGGTCGCCATTCTACTGGTGCGTGA", BCD)
    print(bppm_2.make_row_feature_vector())