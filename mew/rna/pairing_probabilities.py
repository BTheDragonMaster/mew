#!/usr/bin/env python
from sys import argv

from mew.featurisation.viennarna import BPPM
from mew.parsers import read_sequences, read_utr_file


def get_pairing_probabilities(sequence_file, utr_file):
    well_to_sequence = read_sequences(sequence_file)
    five_utr, terminator = read_utr_file(utr_file)

    well_to_bpps = {}

    for well, sequence in well_to_sequence.items():
        bppm = BPPM(sequence, five_utr, terminator)
        bpps = bppm.make_row_feature_vector()
        well_to_bpps[well] = bpps

    return well_to_bpps


def write_pairing_probabilities(well_to_bpps, out_file):
    with open(out_file, 'w') as out:
        for well, bpps in well_to_bpps.items():
            out.write(well)
            for bpp in bpps:
                out.write(f',{bpp:.10f}')
            out.write('\n')


if __name__ == "__main__":
    well_to_bpps = get_pairing_probabilities(argv[1], argv[2])
    write_pairing_probabilities(well_to_bpps, argv[3])