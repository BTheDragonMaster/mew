#!/usr/bin/env python

AA_TO_CODONS = {"F": ["TTT","TTC"],
                "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
                "I": ["ATT", "ATC", "ATA"],
                "M": ["ATG"],
                "V": ["GTT", "GTC", "GTA", "GTG"],
                "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
                "P": ["CCT", "CCC", "CCA", "CCG"],
                "T": ["ACT", "ACC", "ACA", "ACG"],
                "A": ["GCT", "GCC", "GCA", "GCG"],
                "Y": ["TAT", "TAC"],
                "H": ["CAT", "CAC"],
                "Q": ["CAA", "CAG"],
                "N": ["AAT", "AAC"],
                "K": ["AAA", "AAG"],
                "D": ["GAT", "GAC"],
                "E": ["GAA", "GAG"],
                "C": ["TGT", "TGC"],
                "W": ["TGG"],
                "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
                "G": ["GGT", "GGC", "GGA", "GGG"],
                "*": ["TAA", "TAG", "TGA"]}


def reverse_dictionary(dictionary):
    """Return dict of {value: key, ->}

    Input:
    dictionary: dict of {key: [value, ->], ->}
    Output:
    reverse_dictionary: dict of {value: key, ->}

    """
    reverse_dictionary = {}

    for key, values in dictionary.items():
        for value in values:
            reverse_dictionary[value] = key

    return reverse_dictionary

CODON_TO_AA = reverse_dictionary(AA_TO_CODONS)