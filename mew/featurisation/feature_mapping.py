#!/usr/bin/env python

from mew.featurisation.one_hot import CodingSequence, INDEX_TO_BASE, INDEX_TO_CODON
from mew.featurisation.viennarna import get_rna_coordinate


def get_feature_mapping(representative_sequence, encoding_type, five_utr, terminator, length=None, utr_length=None):

    use_terminator = True

    if length != None:
        sequence = representative_sequence[:length]
        use_terminator = False
    else:
        sequence = representative_sequence

    if utr_length != None:
        five_utr = five_utr[len(five_utr) - utr_length:]

    if use_terminator:
        full_sequence = five_utr + sequence + terminator
    else:
        full_sequence = five_utr + sequence

    sequence = CodingSequence(sequence)

    feature_mapping = {}

    if encoding_type == 'one-hot-base':
        for i in range(len(sequence.sequence)):
            for j in range(4):
                feature_mapping[4 * i + j] = (i + 1, INDEX_TO_BASE[j])

    elif encoding_type == 'one-hot-codon':
        for i in range(len(sequence.codons)):
            for j in range(64):
                feature_mapping[64 * i + j] = (i + 1, INDEX_TO_CODON[j])

    elif encoding_type == 'one-hot-third-base':
        for i in range(int(len(sequence.sequence) / 3)):
            for j in range(4):
                feature_mapping[4 * i + j] = (3 * (i + 1), INDEX_TO_BASE[j])

    elif encoding_type == 'rna-bppm-totals':

        for i in range(len(full_sequence)):
            coordinate = get_rna_coordinate(i, five_utr)
            feature_mapping[i] = coordinate

    elif encoding_type == 'rna-bppm':

        index = 0
        for i in range(len(full_sequence)):
            coordinate_1 = get_rna_coordinate(i, five_utr)
            for j in range(i, len(full_sequence)):
                coordinate_2 = get_rna_coordinate(j, five_utr)
                feature_mapping[index] = (coordinate_1, coordinate_2)
                index += 1

    return feature_mapping


def get_feature_string_list(feature_mapping, encoding_type):
    feature_string_list = []
    feature_vector_length = len(list(feature_mapping.keys()))

    for i in range(feature_vector_length):
        feature_string_list.append(get_feature_string(i, feature_mapping, encoding_type))

    return feature_string_list


def get_feature_string(feature_index, feature_mapping, encoding_type):
    if encoding_type == 'one-hot-base' or encoding_type == 'one-hot-third-base':
        base_position, base_type = feature_mapping[feature_index]
        string = f'Base {base_position} is {base_type}'

    elif encoding_type == 'one-hot-codon':
        codon_position, codon_type = feature_mapping[feature_index]
        string = f'Codon {codon_position} is {codon_type}'

    elif encoding_type == 'rna-bppm-totals':
        base_position = feature_mapping[feature_index]
        string = f'Base pairing probability at base {base_position}'

    elif encoding_type == 'rna-bppm':
        base_position_1, base_position_2 = feature_mapping[feature_index]
        string = f'Base pairing probability between bases {base_position_1} and {base_position_2}'

    return string


ENCODING_TO_HEADER = {'one-hot-base': 'Base\tBase type\tFeature importance\n',
                      'one-hot-codon': 'Codon\tCodon type\tFeature importance\n',
                      'one-hot-third-base': 'Base\tBase type\tFeature importance\n',
                      'rna-bppm-totals': 'Base\tFeature importance\n',
                      'rna-bppm': 'Base 1\tBase 2\tFeature importance\n'}
