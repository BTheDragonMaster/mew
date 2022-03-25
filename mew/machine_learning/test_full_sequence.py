import argparse
import os
import joblib
from mew.machine_learning.datasets import build_dataset

def make_parser():
    parser = argparse.ArgumentParser(description="Train a classifier on sequence and flow data.")

    parser.add_argument('-s', '--sequence_data', required=True, type=str, help="Comma-separated file with sequence ID in the left column and sequence in the right column.")
    parser.add_argument('-e', '--expression_data', required=True, type=str, help="Comma-separated file with sequence ID in the left column and normalised expression in the right column.")
    parser.add_argument('-o', '--output_folder', required=True, type=str, help="Output directory for saving model.")
    parser.add_argument('-u', '--utr_data', type=str, default='', help="Two-line file, with the 5' UTR on the first line and the 3' UTR on the second line.")
    parser.add_argument('-n', '--name', required=True, type=str, help="Name for analysis")
    parser.add_argument('-b', '--base_pairing_data', default=None, type=str, help="Comma-separated file with sequence ID in the leftmost column and base pairing probabilities in the remaining columns. If given, base pairing probabilities don't have to be calculated from scratch.")
    parser.add_argument('-l', '--length', default=None, type=int, help="Number of bases in the sequence to be considered as features for machine learning.")
    parser.add_argument('-c', '--classifier', required=True, type=str, help="Directory to classifier")

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


def test_model(args):
    encoding = get_featurisation(args)
    if not os.path.exists(args.output_folder):
        os.mkdir(args.output_folder)
    dataset = build_dataset(args.sequence_data, args.expression_data, args.name, encoding,
                            'test', bpps_file=args.base_pairing_data)

    dataset.set_feature_vectors()
    classifier = joblib.load(args.classifier)

    representative_feature_vector = dataset.data_points[0].feature_vector
    if len(representative_feature_vector) != 0:
        output_folder = os.path.join(args.output_folder, args.name)
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        dataset.test_classifier(classifier)
        dataset.set_vectors_actual_and_predicted_flow()
        dataset.set_correlation_coefficients()
        dataset.plot_actual_vs_predicted(output_folder)
        dataset.write_actual_vs_predicted(output_folder)
        dataset.write_correlation_coefficients(output_folder)

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    test_model(args)