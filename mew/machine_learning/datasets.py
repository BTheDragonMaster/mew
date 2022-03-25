import random
import os

from joblib import dump
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression, Lasso
from sklearn.tree import export_graphviz

import matplotlib
matplotlib.use('TkAgg')
import pydotplus

from scipy.stats import pearsonr, spearmanr

from mew.parsers import read_flow, read_sequences, read_utr_file, read_pairing_probabilities
from mew.machine_learning.data_point import DataPoint
from mew.featurisation.feature_mapping import get_feature_mapping, get_feature_string_list, ENCODING_TO_HEADER
from mew.visualisation.plot import plot_actual_vs_predicted, plot_feature_importances


def build_dataset(sequence_file, flow_file, label, encoding, dataset_type, utr_file=None, bpps_file=None, length=None, utr_length=None, start_position=None, fold=10, algorithm="random_forest", alpha=1.0, coding_length=678):
    well_to_sequence = read_sequences(sequence_file)
    well_to_flow = read_flow(flow_file)

    if utr_file:
        five_utr, terminator = read_utr_file(utr_file)
    else:
        five_utr = ''
        terminator = ''

    if bpps_file:
        well_to_bpps = read_pairing_probabilities(bpps_file)

    else:
        well_to_bpps = {}

    representative_sequence = well_to_sequence[list(well_to_sequence.keys())[0]]

    if dataset_type == 'crossval':
        dataset = CrossvalSet(label, encoding, representative_sequence, five_utr=five_utr, terminator=terminator, fold=fold, algorithm=algorithm, length=length, utr_length=utr_length, alpha=alpha, start_position=start_position, coding_length=coding_length)
    elif dataset_type == 'training':
        dataset = TrainingSet(label, encoding, five_utr=five_utr, terminator=terminator, representative_sequence=representative_sequence, length=length, utr_length=utr_length, start_position=start_position, coding_length=coding_length)
    elif dataset_type == 'test':
        dataset = TestSet(label, encoding, five_utr=five_utr, terminator=terminator, representative_sequence=representative_sequence, length=length, utr_length=utr_length, start_position=start_position, coding_length=coding_length)

    for well, sequence in well_to_sequence.items():
        flow = well_to_flow[well]
        if bpps_file:
            bpps = well_to_bpps[well]
        else:
            bpps = None

        dataset.add_data_point(well, sequence, flow, bpps=bpps)

    return dataset


def get_feature_vectors(data_points):
    feature_vectors = []
    responses = []
    for data_point in data_points:
        feature_vectors.append(data_point.feature_vector)
        responses.append(data_point.flow)

    return feature_vectors, responses


class DataSet:
    def __init__(self, label, encoding, five_utr='', terminator='', representative_sequence=None, length=None, utr_length=None, data_points=None,
                 feature_mapping=None, start_position=None, coding_length=678):
        if not data_points:
            self.data_points = []
        else:
            self.data_points = data_points

        self.label = label
        self.encoding = encoding
        self.five_utr = five_utr
        self.terminator = terminator

        self.length = length
        self.utr_length = utr_length

        self.features = []
        self.responses = []
        self.start_position = start_position
        self.coding_length = coding_length

        if representative_sequence:

            self.feature_mapping = get_feature_mapping(representative_sequence, self.encoding, five_utr, terminator, length=self.length, utr_length=self.utr_length, start_position=self.start_position, coding_length=self.coding_length)
        else:
            if not feature_mapping:
                self.feature_mapping = {}
            else:
                self.feature_mapping = feature_mapping

        self.feature_strings = get_feature_string_list(self.feature_mapping, self.encoding)

        self.feature_importances = {}

        self.true_flow = []
        self.predicted_flow = []

    def __repr__(self):
        return ', '.join([x.__repr__() for x in self.data_points])

    def set_feature_vectors(self):
        self.features, self.responses = get_feature_vectors(self.data_points)

    def add_data_point(self, well, sequence, flow, bpps=None):
        self.data_points.append(DataPoint(well, sequence, flow, self.encoding, bpps=bpps, five_utr=self.five_utr, terminator=self.terminator, length=self.length, utr_length=self.utr_length, start_position=self.start_position, coding_length=self.coding_length))

    def set_vectors_actual_and_predicted_flow(self):

        for data_point in self.data_points:
            self.true_flow.append(data_point.flow)
            self.predicted_flow.append(data_point.predicted_flow)

    def set_correlation_coefficients(self):
        self.pearson, self.pearson_p = pearsonr(self.true_flow, self.predicted_flow)
        self.spearman, self.spearman_p = spearmanr(self.true_flow, self.predicted_flow)

    def plot_actual_vs_predicted(self, out_dir):
        plot_dir = os.path.join(out_dir, 'plots')
        if not os.path.exists(plot_dir):
            os.mkdir(plot_dir)

        plot_actual_vs_predicted(self, plot_dir)

    def write_actual_vs_predicted(self, out_dir):
        out_file = os.path.join(out_dir, self.label + '_actual_vs_predicted.txt')
        with open(out_file, 'w') as out:
            out.write('well\tactual\tpredicted\n')
            for data_point in self.data_points:
                out.write(f'{data_point.well}\t{data_point.flow}\t{data_point.predicted_flow}\n')

    def plot_feature_importances(self, out_dir):
        fi_dir = os.path.join(out_dir, 'feature_importances')
        plot_dir = os.path.join(out_dir, 'plots')
        plot_feature_importances(fi_dir, self.encoding, plot_dir)

    def write_correlation_coefficients(self, out_dir):
        correlation_dir = os.path.join(out_dir, 'correlation_coefficients.txt')
        with open(correlation_dir, 'w') as correlations:
            correlations.write(f'pearson,{self.pearson}\n')
            correlations.write(f'spearman,{self.spearman}\n')

class TestSet(DataSet):
    def __init__(self, label, encoding, five_utr='', terminator='', representative_sequence=None, length=None,
                 utr_length=None, data_points=None, feature_mapping=None, start_position=None, coding_length=678):
        super().__init__(label, encoding, five_utr=five_utr, terminator=terminator,
                         representative_sequence=representative_sequence,
                         length=length,
                         utr_length=utr_length,
                         data_points=data_points,
                         feature_mapping=feature_mapping,
                         start_position=start_position,
                         coding_length=coding_length)

    def test_classifier(self, classifier):
        predicted_responses = classifier.predict(self.features)

        for i, predicted_response in enumerate(predicted_responses):
            data_point = self.data_points[i]
            data_point.set_crossval_predicted_flow(predicted_response)

        return predicted_responses


class TrainingSet(DataSet):
    def __init__(self, label, encoding, five_utr='', terminator='', representative_sequence=None, length=None,
                 utr_length=None, data_points=None, feature_mapping=None, start_position=None, coding_length=678):
        super().__init__(label, encoding, five_utr=five_utr, terminator=terminator,
                         representative_sequence=representative_sequence,
                         length=length,
                         utr_length=utr_length,
                         data_points=data_points,
                         feature_mapping=feature_mapping,
                         start_position=start_position,
                         coding_length=coding_length)

        self.classifier = None

    def train_model(self, algorithm, out_dir, alpha=1.0):
        classifier_dir = os.path.join(out_dir, 'classifier')
        if algorithm == 'random_forest':
            self.train_random_forest(classifier_dir, save_classifier=True)
            #self.save_representative_tree(os.path.join(out_dir, 'crossval_trees'))

        elif algorithm == 'lasso':
            self.train_lasso(classifier_dir, alpha, save_classifier=True)

    def train_random_forest(self, classifier_dir, save_classifier=True):

        self.classifier = RandomForestRegressor()
        print(self.features)
        self.classifier.fit(self.features, self.responses)

        self.feature_importances = self.classifier.feature_importances_
        self.estimators = self.classifier.estimators_

        if not os.path.exists(classifier_dir):
            os.mkdir(classifier_dir)

        if save_classifier:
            dump(self.classifier, os.path.join(classifier_dir, self.label + '_rf.classifier'))

    def save_representative_tree(self, tree_dir):
        if not os.path.exists(tree_dir):
            os.mkdir(tree_dir)

        representative_tree = self.estimators[random.randint(0, 99)]
        tree_output_file = os.path.join(tree_dir, f'{self.label}_representative_tree.dot')
        export_graphviz(representative_tree,
                        out_file=tree_output_file,
                        feature_names=self.feature_strings,
                        rounded=True,
                        filled=True,
                        special_characters=True)

        with open(tree_output_file, 'r') as tree:

            graph = pydotplus.graph_from_dot_data(tree.read())
            graph.write_png(os.path.join(tree_dir, f'{self.label}_representative_tree.png'))

        os.remove(tree_output_file)


    def train_lasso(self, classifier_dir, alpha, save_classifier=True):
        self.classifier = Lasso(alpha=alpha, max_iter=10000)
        self.classifier.fit(self.features, self.responses)

        if not os.path.exists(classifier_dir):
            os.mkdir(classifier_dir)

        self.feature_importances = self.classifier.coef_

        if save_classifier:
            dump(self.classifier, os.path.join(classifier_dir, self.label + '_lasso.classifier'))

    def train_linear_regression(self, classifier_dir, save_classifier=True):
        self.classifier = LinearRegression()
        self.classifier.fit(self.features, self.responses)

        if not os.path.exists(classifier_dir):
            os.mkdir(classifier_dir)

        self.feature_importances = self.classifier.coef_

        if save_classifier:
            dump(self.classifier, os.path.join(classifier_dir, self.label + '_lr.classifier'))

    def train_neural_net(self, classifier_dir, save_classifier=True):
        raise NotImplementedError

    def write_feature_importances(self, feature_importances_dir):
        if not os.path.exists(feature_importances_dir):
            os.mkdir(feature_importances_dir)

        feature_importances_file = os.path.join(feature_importances_dir, f'{self.label}_feature_importances.txt')
        with open(feature_importances_file, 'w') as feature_importances:
            feature_importances.write(ENCODING_TO_HEADER[self.encoding])

            for i, importance in enumerate(self.feature_importances):
                mapping = self.feature_mapping[i]
                if self.encoding == 'one-hot-base' or self.encoding == 'one-hot-codon' or self.encoding == 'one-hot-third-base':
                    position, identity = mapping
                    feature_importances.write(f'{position}\t{identity}\t{importance:.3f}\n')
                elif self.encoding == 'rna-bppm-totals':
                    feature_importances.write(f'{mapping}\t{importance:.3f}\n')
                elif self.encoding == 'rna-bppm':
                    base_1, base_2 = mapping
                    feature_importances.write(f'{base_1}\t{base_2}\t{importance:.3f}\n')
                elif self.encoding == 'rna-bppm-onehot-third':
                    encoding_type = mapping[0]
                    if encoding_type == 'bpp':
                        position = mapping[1]
                        identity = 'N/A'

                    elif encoding_type == 'onehot3':
                        position = mapping[1]
                        identity = mapping[2]

                    feature_importances.write(f'{encoding_type}\t{position}\t{identity}\t{importance:.3f}\n')


class CrossvalSet(DataSet):

    def __init__(self, label, encoding, representative_sequence, five_utr='', terminator='', fold=10,
                 algorithm="random_forest", length=None, utr_length=None, alpha=1.0, start_position=None,
                 coding_length=678):
        super().__init__(label, encoding,
                         five_utr=five_utr,
                         terminator=terminator,
                         representative_sequence=representative_sequence,
                         length=length,
                         utr_length=utr_length,
                         start_position=start_position,
                         coding_length=coding_length)

        self.fold = fold

        self.groups = {}
        self.initialise_groups()

        self.feature_importances = {}
        self.algorithm = algorithm

        self.true_flow = []
        self.predicted_flow = []
        self.alpha = alpha

    def initialise_groups(self):
        for i in range(self.fold):
            self.groups[i] = []

    def populate_groups(self):

        random.shuffle(self.data_points)

        for i, data_point in enumerate(self.data_points):
            self.groups[i % self.fold].append(data_point)

    def get_train_and_test_sets(self, out_group):
        test_set = TestSet(f'test_crossval_{out_group}', self.encoding, feature_mapping=self.feature_mapping, five_utr=self.five_utr, terminator=self.terminator, length=self.length, utr_length=self.utr_length)
        training_set = TrainingSet(f'training_crossval_{out_group}', self.encoding, feature_mapping=self.feature_mapping, five_utr=self.five_utr, terminator=self.terminator, length=self.length, utr_length=self.utr_length)

        for group in range(self.fold):
            if group == out_group:
                for data_point in self.groups[group]:

                    test_set.data_points.append(data_point)
            else:
                for data_point in self.groups[group]:
                    training_set.data_points.append(data_point)

        test_set.set_feature_vectors()
        training_set.set_feature_vectors()

        return test_set, training_set

    def do_crossval(self, out_dir, save_classifiers=True):

        for group in range(self.fold):
            test_set, training_set = self.get_train_and_test_sets(group)

            classifier_dir = os.path.join(os.getcwd(), os.path.join(out_dir, 'crossval_classifiers'))

            if self.algorithm == 'random_forest':
                training_set.train_random_forest(classifier_dir, save_classifier=save_classifiers)
                if save_classifiers:
                    training_set.save_representative_tree(os.path.join(out_dir, 'crossval_trees'))

            elif self.algorithm == 'lasso':
                training_set.train_lasso(classifier_dir, self.alpha, save_classifier=save_classifiers)

            elif self.algorithm == 'neural_net':
                training_set.train_neural_net(classifier_dir, save_classifier=save_classifiers)

            training_set.write_feature_importances(os.path.join(out_dir, 'feature_importances'))

            predicted_responses = test_set.test_classifier(training_set.classifier)

            for i, predicted_response in enumerate(predicted_responses):
                data_point = self.groups[group][i]
                data_point.set_crossval_predicted_flow(predicted_response)

        self.set_vectors_actual_and_predicted_flow()
        self.set_correlation_coefficients()


    def set_vectors_actual_and_predicted_flow(self):

        for data_point in self.data_points:
            self.true_flow.append(data_point.flow)
            self.predicted_flow.append(data_point.predicted_flow)

    def set_correlation_coefficients(self):
        self.pearson, self.pearson_p = pearsonr(self.true_flow, self.predicted_flow)
        self.spearman, self.spearman_p = spearmanr(self.true_flow, self.predicted_flow)

    def plot_actual_vs_predicted(self, out_dir):
        plot_dir = os.path.join(out_dir, 'plots')
        if not os.path.exists(plot_dir):
            os.mkdir(plot_dir)

        plot_actual_vs_predicted(self, plot_dir)

    def write_actual_vs_predicted(self, out_dir):
        out_file = os.path.join(out_dir, self.label + '_actual_vs_predicted.txt')
        with open(out_file, 'w') as out:
            out.write('well\tactual\tpredicted\n')
            for data_point in self.data_points:
                out.write(f'{data_point.well}\t{data_point.flow}\t{data_point.predicted_flow}\n')

    def plot_feature_importances(self, out_dir):
        fi_dir = os.path.join(out_dir, 'feature_importances')
        plot_dir = os.path.join(out_dir, 'plots')
        plot_feature_importances(fi_dir, self.encoding, plot_dir)

    def write_correlation_coefficients(self, out_dir):
        correlation_dir = os.path.join(out_dir, 'correlation_coefficients.txt')
        with open(correlation_dir, 'w') as correlations:
            correlations.write(f'pearson,{self.pearson}\n')
            correlations.write(f'spearman,{self.spearman}\n')





