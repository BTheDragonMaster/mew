import random

from mew.machine_learning.data_point import DataPoint



class TrainingSet:
    pass


class CrossvalSet:
    def __init__(self, label, featurisation, fold=10):

        self.label = label
        self.featurisation = featurisation
        self.fold = fold

        self.data_points = []
        self.groups = {}
        self.make_groups()

        self.feature_importances = {}

    def initialise_groups(self):
        for i in range(self.fold):
            self.groups[i] = []

    def populate_groups(self):

        random.shuffle(self.data_points)

        for i, data_point in enumerate(self.data_points):
            self.groups[i % self.fold].append(data_point)

    def add_data_point(self, well, sequence, flow):
        self.data_points.append(DataPoint(well, sequence, flow))



    def make_training_test_set_crossval(self, test_group):
        self.test = []
        self.training = []

        for group in range(self.fold):
            if group == test_group:
                self.test += self.groups[group]

            else:
                self.training += self.groups[group]

        self.test_features, self.test_responses = make_feature_vectors(self.test)
        self.training_features, self.training_responses = make_feature_vectors(self.training)

    def train_classifier(self, classifier_dir):
        training_features, training_responses = make_feature_vectors(self.data_points)
        test_features, test_responses = make_feature_vectors(test_set.data_points)
        self.classifier = RandomForestRegressor()
        self.classifier.fit(training_features, training_responses)

        dump(self.classifier, classifier_dir + self.label + '_rf.classifier')

    def test_classifier(self, test_set):
        test_features, test_responses = make_feature_vectors(test_set.data_points)
        predicted_responses = self.classifier.predict(test_features)

        for i, predicted_response in enumerate(predicted_responses):
            data_point = test_set.data_points[i]
            data_point.set_crossval_predicted_flow(predicted_response)

    def do_crossval(self, tree_dir):
        for group in range(10):
            self.make_training_test_set_crossval(group)

            regressor = RandomForestRegressor()
            regressor.fit(self.training_features, self.training_responses)
            self.feature_importances[group] = regressor.feature_importances_
            representative_tree = regressor.estimators_[random.randint(0, 99)]
            tree_output_file = f'{tree_dir}{self.label}_representative_tree_{group}.dot'
            export_graphviz(representative_tree,
                            out_file=tree_output_file,
                            feature_names=self.feature_names,
                            rounded=True,
                            filled=True,
                            special_characters=True)
            with open(tree_output_file) as tree:
                graph = pydotplus.graph_from_dot_data(tree.read())
                graph.write_png(f'{tree_dir}{self.label}_representative_tree_{group}.png')

            predicted_responses = regressor.predict(self.test_features)
            for i, predicted_response in enumerate(predicted_responses):
                data_point = self.groups[group][i]
                data_point.set_crossval_predicted_flow(predicted_response)

    def write_feature_importances(self, out_dir):
        out_dir = f'{out_dir}_feature_importances.txt'
        out_file = open(out_dir, 'w')
        out_file.write('Codon_nr')

        for group in range(10):
            out_file.write('\t%d' % (group + 1))

        out_file.write('\n')

        codon_nr = 0
        for i, feature in enumerate(self.feature_importances[0]):

            if i % 12 == 0:
                codon_nr += 1

            out_file.write(str(codon_nr))

            for group in range(10):
                feature_importance = self.feature_importances[group][i]
                out_file.write('\t%.3f' % feature_importance)
            out_file.write('\n')

    def write_feature_importances_6(self, out_dir):
        out_file = open(out_dir, 'w')
        out_file.write('Codon_nr')

        for group in range(10):
            out_file.write('\t%d' % (group + 1))

        out_file.write('\n')

        codon_nr = 1

        aa_nr = 0
        current_aa = mRFP[aa_nr]
        nr_of_options = len(aa_dict[current_aa])
        option_nr = 0

        for i, feature in enumerate(self.feature_importances[0]):

            out_file.write(str(codon_nr))
            option_nr += 1

            if nr_of_options == option_nr:
                try:
                    codon_nr += 1
                    aa_nr += 1
                    current_aa = mRFP[aa_nr]
                    option_nr = 0
                    nr_of_options = len(aa_dict[current_aa])
                except IndexError:
                    pass

            for group in range(10):
                feature_importance = self.feature_importances[group][i]
                out_file.write('\t%.3f' % feature_importance)
            out_file.write('\n')

        out_file.close()