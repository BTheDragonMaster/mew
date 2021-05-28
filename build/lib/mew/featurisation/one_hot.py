#!/usr/bin/env python

from mew.data.codon_tables import AA_TO_CODONS, CODON_TO_AA



def do_onehot_encoding(sequence, type='one-hot-base'):
    sequence = CodingSequence(sequence)
    feature_vector = sequence.get_feature_vector(type)
    return feature_vector

def make_onehot_codon_mapping():
    codon_to_vector = {}
    codon_to_index = {}
    index_to_codon = {}
    index = 0

    for letter_1 in ['A', 'C', 'G', 'T']:
        for letter_2 in ['A', 'C', 'G', 'T']:
            for letter_3 in ['A', 'C', 'G', 'T']:
                codon = letter_1 + letter_2 + letter_3
                codon_to_vector[codon] = [0] * 64
                codon_to_index[codon] = index
                index_to_codon[index] = codon
                codon_to_vector[codon][index] = 1
                index += 1

    return codon_to_vector, codon_to_index, index_to_codon


CODON_TO_VECTOR, CODON_TO_INDEX, INDEX_TO_CODON = make_onehot_codon_mapping()

BASE_TO_VECTOR = {'A': [1, 0, 0, 0],
                  'C': [0, 1, 0, 0],
                  'G': [0, 0, 1, 0],
                  'T': [0, 0, 0, 1]}

BASE_TO_INDEX = {'A': 0,
                 'C': 1,
                 'G': 2,
                 'T': 3}

INDEX_TO_BASE = {0: 'A',
                 1: 'C',
                 2: 'G',
                 3: 'T'}



class CodingSequence:

    def __init__(self, sequence):
        self.sequence = sequence
        assert len(self.sequence) % 3 == 0
        
        self.codons = self.get_codons()

    def get_feature_vector(self, encoding_type='one-hot-base'):
        if encoding_type == 'one-hot-base':
            feature_vector= self.make_onehot_encoded_base_vector()
        elif encoding_type == 'one-hot-codon':
            feature_vector = self.make_onehot_encoded_codon_vector()
        elif encoding_type == 'one-hot-third-base':
            feature_vector = self.make_onehot_encoded_third_base_vector()

        return feature_vector

    def get_codons(self):
        codons = []

        for i in range(0, len(self.sequence), 3):
            codon = self.sequence[i:i + 3]
            codons.append(codon)

        return codons

    def make_onehot_encoded_third_base_vector(self):
        vector = []

        for i, base in enumerate(self.sequence):
            if i % 3 == 2:
                vector += BASE_TO_VECTOR[base]

        return vector


    def make_onehot_encoded_codon_vector(self):
        vector = []

        for i, codon in enumerate(self.codons):
            vector += CODON_TO_VECTOR[codon]

        return vector

    def make_onehot_encoded_base_vector(self):
        vector = []

        for i, base in enumerate(self.sequence):
            vector += BASE_TO_VECTOR[base]

        return vector

    def get_sliding_windows_codons(self, window_size):
        windows = []
        for i in range(0, len(self.codons) - window_size):
            codons = self.codons[i: i + window_size]
            windows.append(codons)

        return windows

    def get_sliding_windows_bases(self, window_size):

        windows = []
        for i in range(0, len(self.sequence) - window_size):
            bases = self.sequence[i: i + window_size]
            windows.append(bases)

        return windows

    def get_windows_codons(self, window_size):

        windows = []
        for i in range(0, len(self.codons), window_size):
            codons = self.codons[i: i + window_size]
            windows.append(codons)

        return windows

    def get_windows_bases(self, window_size):
        windows = []
        for i in range(0, len(self.sequence), window_size):
            bases = self.sequence[i: i + window_size]
            windows.append(bases)

        return windows
    
    
        
        
        

