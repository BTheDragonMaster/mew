#!/usr/bin/env python
import os
import argparse

from mew.parsers import read_sequences, read_pairing_probabilities, read_utr_file
from mew.featurisation.viennarna import get_rna_coordinate

def make_parser():
    parser = argparse.ArgumentParser(description="Divide sequences into windows and save them into files.")

    parser.add_argument('-s', '--sequence_data', required=True, type=str, help="Comma-separated file with sequence ID in the left column and sequence in the right column.")
    parser.add_argument('-w', '--window_size', required=True, type=int, help="Size of sequence windows to be analysed with machine learning.")
    parser.add_argument('-utr', '--utr_file', required=True, type=str, help="Two-line file, with the 5' UTR on the first line and the 3' UTR on the second line.")
    parser.add_argument('-bpp', '--base_pairing_data', required=True, type=str, help="Comma-separated file with sequence ID in the leftmost column and base pairing probabilities in the remaining columns.")
    parser.add_argument('-o', '--output', required=True, type=str, help="Output directory")
    return parser

def get_windows(sequence, window_size):
    windows = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        windows.append(window)
    return windows


def write_windows(sequence_file, base_pairing_file, window_size, out_dir, utr_file):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    window_out_dir = os.path.join(out_dir, 'sequence_windows')
    bpp_out_dir = os.path.join(out_dir, 'bpp_windows')

    if not os.path.exists(window_out_dir):
        os.mkdir(window_out_dir)
    if not os.path.exists(bpp_out_dir):
        os.mkdir(bpp_out_dir)

    well_to_sequence = read_sequences(sequence_file)
    well_to_probabilities = read_pairing_probabilities(base_pairing_file)
    five_utr, terminator = read_utr_file(utr_file)

    well_to_windows = {}
    well_to_bpp_windows = {}
    window_nr = 0

    for well, sequence in well_to_sequence.items():
        bpps = well_to_probabilities[well]
        full_sequence = five_utr + sequence + terminator
        windows = get_windows(full_sequence, window_size)
        bpp_windows = get_windows(bpps, window_size)

        if not window_nr:
            window_nr = len(windows)

        well_to_windows[well] = windows
        well_to_bpp_windows[well] = bpp_windows

    for i in range(window_nr):
        window_coordinate = i - len(five_utr)
        window_file = os.path.join(window_out_dir, f'window_{window_coordinate}.txt')
        bpp_file = os.path.join(bpp_out_dir, f'window_{window_coordinate}.txt')

        with open(window_file, 'w') as window_f:
            for well, windows in well_to_windows.items():
                window = windows[i]
                window_f.write(f'{well},{window}\n')

        with open(bpp_file, 'w') as bpp_f:
            for well, bpps_windows in well_to_bpp_windows.items():
                bpps_window = bpps_windows[i]
                bpp_f.write(f'{well}')
                for bpp in bpps_window:
                    bpp_f.write(f',{bpp:.10f}')
                bpp_f.write('\n')


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    write_windows(args.sequence_data, args.base_pairing_data, args.window_size, args.output, args.utr_file)

