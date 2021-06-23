#!/usr/bin/env python

import RNA


def get_rna_coordinate(rna_coordinate, five_utr):
    coordinate = rna_coordinate - len(five_utr) + 1
    if coordinate <= 0:
        coordinate -= 1

    return coordinate


def from_bppm_base_totals(sequence, five_utr='', terminator=''):
    """Return feature vector and mapping, with each feature the total probability of a base's involvement in pairing.

    :param sequence: str, DNA/RNA sequence of gene
    :param five_utr: str, DNA/RNA sequence of upstream sequence

    :return feature_vector: list of float, with each float the total probability of a base's involvement in pairing
    :return feature_mapping: dict of {index: base position, ->}, with index and base position int. Index corresponds to
                             index in the feature vector
    """

    bppm = BPPM(sequence, five_utr, terminator)
    feature_vector = bppm.make_row_feature_vector()
    return feature_vector


def from_bppm(sequence, five_utr='', terminator=''):
    """Return feature vector and mapping, with each feature the probability of two bases pairing.

    :param sequence: str, DNA/RNA sequence of gene
    :param five_utr: str, DNA/RNA sequence of upstream sequence

    :return feature_vector: list of float, with each float the probability of two bases pairing
    :return feature_mapping: dict of {index: (base 1, base 2), ->}, with index, base 1 and base 2 int. Base 1 and
                             base 2 are the bases between which the pairing probability at feature_vector[index] was
                             calculated.
    """
    bppm = BPPM(sequence, five_utr, terminator)
    feature_vector = bppm.make_feature_vector()
    return feature_vector


class BPPM:
    def __init__(self, sequence, five_utr, terminator):
        self.full_sequence = five_utr + sequence + terminator
        self.sequence = sequence
        self.five_utr = five_utr
        self.terminator = terminator

        print(f"Folding sequence starting with {sequence[:20]}...")

        fc = RNA.fold_compound(self.full_sequence)
        self.propensity, self.ensemble_energy = fc.pf()
        self.bppm = fc.bpp()

    def make_feature_vector(self):
        index = 0
        feature_vector = []
        for i in range(1, len(self.full_sequence)):
            for j in range(i + 1, len(self.full_sequence) + 1):
                feature_vector.append(self.bppm[i][j])
                index += 1

        return feature_vector

    def make_row_feature_vector(self):
        index = 0
        feature_vector = []
        for i in range(1, len(self.full_sequence) + 1):
            bp_prob = sum(self.bppm[i])
            feature_vector.append(bp_prob)
            index += 1

        return feature_vector

    def print_bppm(self):
        for i in range(1, len(self.full_sequence)):
            for j in range(i + 1, len(self.full_sequence) + 1):
                print("pr(%d,%d) = %g" % (i - len(self.five_utr), j - len(self.five_utr), self.bppm[i][j]))
