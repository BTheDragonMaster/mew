#!/usr/bin/env python

import argparse
import os
from multiprocessing import Pool

from mew.machine_learning.datasets import *


def make_parser():
    parser = argparse.ArgumentParser(description="Train and cross-validate one or more classifiers on sequence and flow data.")

    parser.add_argument('-s', '--sequence_data_dir', required=True, type=str, help="Comma-separated file with sequence ID in the left column and sequence in the right column.")
    parser.add_argument('-e', '--expression_data', required=True, type=str, help="Comma-separated file with sequence ID in the left column and normalised expression in the right column.")
    parser.add_argument('-o', '--output_folder', default=os.path.join(os.getcwd(), 'cross-validation'), type=str, help="Output folder for saving cross-validation data.")
    parser.add_argument('-f', '--fold', default=10, type=int, help="Supply x for x-fold cross-validation.")
    parser.add_argument('-n', '--name', required=True, type=str, help="Name for analysis")
    parser.add_argument('-b', '--base_pairing_data', default=None, type=str, help="Comma-separated file with sequence ID in the leftmost column and base pairing probabilities in the remaining columns. If given, base pairing probabilities don't have to be calculated from scratch.")
    parser.add_argument('-alpha', '--alpha', default=1.0, type=float, help="Alpha value for LASSO.")
    parser.add_argument('-cd', '--coding_length', default=678, type=int, help="Length of coding sequence. Default value represents length of mRFP")
    parser.add_argument('-t', '--threads', default=1, type=int, help="Number of threads.")

    ml_methods = parser.add_mutually_exclusive_group(required=True)
    ml_methods.add_argument('-la', '--lasso', action='store_true', help="Do cross-validation for LASSO regressor.")
    ml_methods.add_argument('-rf', '--random_forest', action='store_true', help="Do cross-validation for random forest regressor.")
    ml_methods.add_argument('-lr', '--linear_regression', action='store_true', help="Do cross-validation for linear regression.")
    ml_methods.add_argument('-nn', '--neural_net', action='store_true', help="Do cross-validation for neural net regressor.")

    featurisation_methods = parser.add_mutually_exclusive_group(required=True)
    featurisation_methods.add_argument('-ob', '--onehot_base', action='store_true', help="Use one-hot encoding of bases for featurisation.")
    featurisation_methods.add_argument('-o3', '--onehot_third_base', action='store_true', help="Use one-hot encoding of every third base for featurisation.")
    featurisation_methods.add_argument('-oc', '--onehot_codon', action='store_true', help="Use one-hot encoding of codons for featurisation.")
    featurisation_methods.add_argument('-bpp', '--base_pairing_probabilities', action='store_true', help="Use base-pairing probabilities derived from mRNA secondary structure for featurisation.")
    featurisation_methods.add_argument('-bppo3', '--base_pairing_and_onehot_third_base', action='store_true', help="Use base-pairing probabilities derived from mRNA secondary structure and one-hot encoding of every third base for featurisation.")

    return parser


def get_featurisation(args):

    featurisation = None

    if args.onehot_base:
        featurisation = 'one-hot-base'
    elif args.onehot_third_base:
        featurisation = 'one-hot-third-base'
    elif args.onehot_codon:
        featurisation = 'one-hot-codon'
    elif args.base_pairing_probabilities:
        featurisation = 'rna-bppm-totals'
    elif args.base_pairing_and_onehot_third_base:
        featurisation = 'rna-bppm-onehot-third'

    return featurisation


def get_algorithm(args):
    algorithm = None

    if args.lasso:
        algorithm = 'lasso'
    elif args.random_forest:
        algorithm = 'random_forest'
    elif args.linear_regression:
        algorithm = 'linear_regression'
    elif args.neural_net:
        algorithm = 'neural_net'

    return algorithm


def do_crossval_single(args, encoding, algorithm, sequence_file):
    window_number = int(sequence_file.split('.')[0].split('_')[1])
    full_path = os.path.join(args.sequence_data_dir, sequence_file)
    base_pairing_file = os.path.join(args.base_pairing_data, sequence_file)
    dataset = build_dataset(full_path, args.expression_data, f'{args.name}_window_{window_number}', encoding,
                            'crossval', bpps_file=base_pairing_file, fold=args.fold, algorithm=algorithm,
                            alpha=args.alpha, start_position=window_number, coding_length=args.coding_length)
    dataset.initialise_groups()
    dataset.populate_groups()

    output_folder = os.path.join(args.output_folder, f'{args.name}_window_{window_number}')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    dataset.do_crossval(output_folder, save_classifiers=False)
    dataset.plot_actual_vs_predicted(output_folder)
    dataset.write_actual_vs_predicted(output_folder)
    dataset.plot_feature_importances(output_folder)
    dataset.write_correlation_coefficients(output_folder)



def do_crossval_multi(args):
    encoding = get_featurisation(args)
    algorithm = get_algorithm(args)
    argument_list = []
    for sequence_file in os.listdir(args.sequence_data_dir):
        if sequence_file[-4:] == '.txt':
            arguments = (args, encoding, algorithm, sequence_file)
            argument_list.append(arguments)

    pools = Pool(args.threads)
    pools.starmap(do_crossval_single, argument_list)

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    do_crossval_multi(args)
