#!/usr/bin/env python

from mew.featurisation.one_hot import CodingSequence, INDEX_TO_BASE, INDEX_TO_CODON, get_coding_range
from mew.featurisation.viennarna import get_rna_coordinate


def get_feature_mapping(representative_sequence, encoding_type, five_utr, terminator, start_position=None, length=None, utr_length=None, coding_length=678):

    if start_position == None:
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

        elif encoding_type == 'rna-bppm-onehot-third':
            for i in range(len(full_sequence)):
                coordinate = get_rna_coordinate(i, five_utr)
                feature_mapping[i] = ('bpp', coordinate)

            for i in range(int(len(sequence.sequence) / 3)):
                for j in range(4):
                    feature_mapping[len(full_sequence) + 4 * i + j] = ('onehot3', 3 * (i + 1), INDEX_TO_BASE[j])

    else:
        sequence = get_coding_range(representative_sequence, start_position, coding_length=coding_length)

        if start_position < 0:
            five_utr = representative_sequence[:abs(start_position)]
        else:
            five_utr = ''

        terminator_start = len(representative_sequence) + start_position
        if terminator_start > coding_length:
            terminator = representative_sequence[coding_length - start_position:]
        else:
            terminator = ''

        feature_mapping = {}

        codon_start_position = max(0, start_position / 3)


        if encoding_type == 'one-hot-base':
            for i in range(len(sequence.sequence)):
                for j in range(4):
                    feature_mapping[4 * i + j] = (max(0, start_position) + i + 1, INDEX_TO_BASE[j])

        elif encoding_type == 'one-hot-codon':
            for i in range(len(sequence.codons)):
                for j in range(64):
                    feature_mapping[64 * i + j] = (codon_start_position + i + 1, INDEX_TO_CODON[j])

        elif encoding_type == 'one-hot-third-base':
            counter = 0
            for i in range(len(sequence.sequence)):
                if i % 3 == 2 - sequence.frame:
                    for j in range(4):
                        feature_mapping[4 * counter + j] = (max(0, start_position) + i + 1, INDEX_TO_BASE[j])
                    counter += 1

        elif encoding_type == 'rna-bppm-totals':

            for i in range(len(representative_sequence)):
                if start_position >= 0:
                    coordinate = i + start_position + 1
                else:
                    coordinate = get_rna_coordinate(i, five_utr)
                feature_mapping[i] = coordinate

        elif encoding_type == 'rna-bppm':


            index = 0
            for i in range(len(representative_sequence)):
                if start_position >= 0:
                    coordinate_1 = i + start_position + 1
                else:
                    coordinate_1 = get_rna_coordinate(i, five_utr)
                for j in range(i, len(representative_sequence)):
                    if start_position >= 0:
                        coordinate_2 = j + start_position + 1
                    else:
                        coordinate_2 = get_rna_coordinate(j, five_utr)
                    feature_mapping[index] = (coordinate_1, coordinate_2)
                    index += 1

        elif encoding_type == 'rna-bppm-onehot-third':
            for i in range(len(representative_sequence)):
                if start_position >= 0:
                    coordinate = i + start_position + 1
                else:
                    coordinate = get_rna_coordinate(i, five_utr)
                feature_mapping[i] = ('bpp', coordinate)

            counter = 0

            for i in range(len(sequence.sequence)):
                if i % 3 == 2 - sequence.frame:
                    for j in range(4):
                        feature_mapping[len(representative_sequence) + 4 * counter + j] = ('onehot3', max(0, start_position) + i + 1, INDEX_TO_BASE[j])
                    counter += 1

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

    elif encoding_type == 'rna-bppm-onehot-third':
        encoding_type = feature_mapping[feature_index][0]
        if encoding_type == 'bpp':
            base_position = feature_mapping[feature_index][1]
            string = f'Base pairing probability at base {base_position}'
        elif encoding_type == 'onehot3':
            base_position, base_type = feature_mapping[feature_index][1:]
            string = f'Base {base_position} is {base_type}'

    return string


ENCODING_TO_HEADER = {'one-hot-base': 'Base\tBase type\tFeature importance\n',
                      'one-hot-codon': 'Codon\tCodon type\tFeature importance\n',
                      'one-hot-third-base': 'Base\tBase type\tFeature importance\n',
                      'rna-bppm-totals': 'Base\tFeature importance\n',
                      'rna-bppm': 'Base 1\tBase 2\tFeature importance\n',
                      'rna-bppm-onehot-third': 'Encoding\tBase\tBase type\tFeature importance\n'}
