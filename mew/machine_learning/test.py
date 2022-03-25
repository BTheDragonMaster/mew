import argparse
import os
import joblib
from mew.machine_learning.datasets import build_dataset

def make_parser():
    parser = argparse.ArgumentParser(description="Train a classifier on sequence and flow data.")

    parser.add_argument('-s', '--sequence_data', required=True, type=str, help="Comma-separated file with sequence ID in the left column and sequence in the right column.")
    parser.add_argument('-e', '--expression_data', required=True, type=str, help="Comma-separated file with sequence ID in the left column and normalised expression in the right column.")
    parser.add_argument('-o', '--output_folder', required=True, type=str, help="Output directory for saving model.")
    parser.add_argument('-n', '--name', required=True, type=str, help="Name for analysis")
    parser.add_argument('-b', '--base_pairing_data', default=None, type=str, help="Comma-separated file with sequence ID in the leftmost column and base pairing probabilities in the remaining columns. If given, base pairing probabilities don't have to be calculated from scratch.")
    parser.add_argument('-l', '--length', default=None, type=int, help="Number of bases in the sequence to be considered as features for machine learning.")
    parser.add_argument('-alpha', '--alpha', default=1.0, type=float, help="Alpha value for LASSO.")
    parser.add_argument('-w', '--window', required=True, type=int, help="Number of analysed window")
    parser.add_argument('-c', '--classifier', required=True, type=str, help="Directory to classifier")

    ml_methods = parser.add_mutually_exclusive_group(required=True)
    ml_methods.add_argument('-la', '--lasso', action='store_true', help="Train LASSO regressor.")
    ml_methods.add_argument('-rf', '--random_forest', action='store_true', help="Train random forest regressor.")
    ml_methods.add_argument('-lr', '--linear_regression', action='store_true', help="Train linear regressor.")
    ml_methods.add_argument('-nn', '--neural_net', action='store_true', help="Train neural net regressor.")

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


def test_model(args):
    encoding = get_featurisation(args)
    algorithm = get_algorithm(args)
    if not os.path.exists(args.output_folder):
        os.mkdir(args.output_folder)
    dataset = build_dataset(args.sequence_data, args.expression_data, args.name, encoding,
                            'test', bpps_file=args.base_pairing_data, algorithm=algorithm,
                            start_position=args.window, alpha=args.alpha)

    dataset.set_feature_vectors()
    classifier = joblib.load(args.classifier)

    representative_feature_vector = dataset.data_points[0].feature_vector
    if len(representative_feature_vector) != 0:
        dataset.test_classifier(classifier)
        dataset.set_vectors_actual_and_predicted_flow()
        dataset.set_correlation_coefficients()
        dataset.plot_actual_vs_predicted(args.output_folder)
        dataset.write_actual_vs_predicted(args.output_folder)
        dataset.write_correlation_coefficients(args.output_folder)

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    test_model(args)