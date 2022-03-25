#!/usr/bin/env python
import os
from sys import argv
import argparse
from sklearn.model_selection import train_test_split
from mew.parsers import read_sequences, read_flow, read_pairing_probabilities
from mew.writers import write_sequences, write_flow, write_bpp


def make_parser():
    parser = argparse.ArgumentParser(description="Split dataset into train and test set.")
    parser.add_argument('-s', '--sequences', required=True, type=str, help="Comma-separated file with sequence ID in the left column and sequence in the right column.")
    parser.add_argument('-e', '--expression_data', required=True, type=str, help="Comma-separated file with sequence ID in the left column and normalised expression in the right column.")
    parser.add_argument('-o', '--output_folder', default=os.path.join(os.getcwd(), 'cross-validation'), type=str, help="Output folder for saving cross-validation data.")
    parser.add_argument('-f', '--fraction', default=0.1, type=float, help="Fraction of data to be used as test set.")
    parser.add_argument('-b', '--base_pairing_data', default=None, type=str, help="Comma-separated file with sequence ID in the leftmost column and base pairing probabilities in the remaining columns. If given, base pairing probabilities don't have to be calculated from scratch.")

    return parser


def stratify_data(seq_ids, fraction):
    datasets = []
    for seq_id in seq_ids:
        dataset = seq_id[0]
        datasets.append(dataset)

    x_train, x_test = train_test_split(seq_ids, test_size=fraction, random_state=64810, stratify=datasets)

    return x_train, x_test


def split_train_test(sequences, expression_data, bpp_data, fraction, out_folder):
    well_to_sequence = read_sequences(sequences)
    wells = list(well_to_sequence.keys())
    train_wells, test_wells = stratify_data(wells, fraction)

    well_to_flow = read_flow(expression_data)
    well_to_bpp = read_pairing_probabilities(bpp_data)

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    out_train = os.path.join(out_folder, 'training')
    out_test = os.path.join(out_folder, 'test')

    training_sequences = os.path.join(out_train, 'sequence_data.csv')
    test_sequences = os.path.join(out_test, 'sequence_data.csv')

    training_flow = os.path.join(out_train, 'flow_data.csv')
    test_flow = os.path.join(out_test, 'flow_data.csv')

    training_bpp = os.path.join(out_train, 'bpp_data.csv')
    test_bpp = os.path.join(out_test, 'bpp_data.csv')

    train_well_to_sequence = {}
    test_well_to_sequence = {}

    train_well_to_flow = {}
    test_well_to_flow = {}

    train_well_to_bpp = {}
    test_well_to_bpp = {}

    for well in train_wells:
        sequence = well_to_sequence[well]
        flow = well_to_flow[well]
        bpp = well_to_bpp[well]

        train_well_to_sequence[well] = sequence
        train_well_to_flow[well] = flow
        train_well_to_bpp[well] = bpp

    for well in test_wells:
        sequence = well_to_sequence[well]
        flow = well_to_flow[well]
        bpp = well_to_bpp[well]

        test_well_to_sequence[well] = sequence
        test_well_to_flow[well] = flow
        test_well_to_bpp[well] = bpp

    write_sequences(train_well_to_sequence, training_sequences)
    write_sequences(test_well_to_sequence, test_sequences)
    write_flow(train_well_to_flow, training_flow)
    write_flow(test_well_to_flow, test_flow)
    write_bpp(train_well_to_bpp, training_bpp)
    write_bpp(test_well_to_bpp, test_bpp)


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()

    split_train_test(args.sequences, args.expression_data, args.base_pairing_data, args.fraction, args.output_folder)


