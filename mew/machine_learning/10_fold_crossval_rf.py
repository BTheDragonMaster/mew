from parsers import *
import random
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import export_graphviz

from sys import argv, path
import os
from joblib import dump, load
from statistics import mean

from pprint import pprint
import argparse
from one_hot_encoding import *
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import scipy
from scipy.stats import pearsonr
from IPython.display import Image
import pydotplus

mRFP = 'MASSEDVIKEFMRFKVRMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFQYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASTERMYPEDGALKGEIKMRLKLKDGGHYDAEVKTTYMAKKPVQLPGAYKTDIKLDITSHNEDYTIVEQYERAEGRHSTGA*'

class TrainingSet():
    def __init__(self, label):

        self.label = label
        self.data_points = []

        self.groups = {0: [],
                       1: [],
                       2: [],
                       3: [],
                       4: [],
                       5: [],
                       6: [],
                       7: [],
                       8: [],
                       9: []}

        self.feature_importances = {}
        self.names_vector_created = False

    def make_feature_names_vector(self, sequence):
        self.feature_names = []
        for i in range(len(sequence.sequence)):
            codon = 1 + i / 3
            base = 1 + i % 3
            self.feature_names.append('base %d of codon %d is A' % (base, codon))
            self.feature_names.append('base %d of codon %d is C' % (base, codon))
            self.feature_names.append('base %d of codon %d is G' % (base, codon))
            self.feature_names.append('base %d of codon %d is T' % (base, codon))
            
            

    def add_data_point(self, well, sequence, flow):
        self.data_points.append(DataPoint(well, sequence, flow))
        if not self.names_vector_created:
            self.make_feature_names_vector(sequence)
            self.names_vector_created = True


    def divide_into_groups(self):

        random.shuffle(self.data_points)

        for i, data_point in enumerate(self.data_points):
            self.groups[i % 10].append(data_point)
                         

    def make_training_test_set_crossval(self, out_group):
        self.test = []
        self.training = []
        
        for group in range(10):
            if group == out_group:
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
                            out_file = tree_output_file,
                            feature_names = self.feature_names,
                            rounded = True,
                            filled = True,
                            special_characters = True)
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

            

class DataPoint():
    def __init__(self, well, sequence, flow):
        self.flow = flow
        self.well = well
        self.sequence = sequence
        self.onehot_12 = sequence.onehot_encode_base_12()
        self.onehot = sequence.onehot_encode()
        self.onehot_6 = sequence.onehot_encode_6()

    def __repr__(self):
        return self.well

    def __hash__(self):
        return self.well

    def set_crossval_predicted_flow(self, flow):
        self.predicted_flow = flow
        self.difference = abs(self.predicted_flow - self.flow)

    
def make_training_set(flow_dict, translation_results):
    wells = list(flow_dict.keys())
    wells.sort()

    training_set = TrainingSet('all')
    low_training_set = TrainingSet('low')
    medium_training_set = TrainingSet('medium')
    high_training_set = TrainingSet('high')
    highlow_training_set = TrainingSet('highlow')


    for well in wells:
        flow = flow_dict[well]
        sequence = CodingSequence(translation_results.results[well].orf)#[:24])
        training_set.add_data_point(well, sequence, flow)

        if well.startswith('L'):
            low_training_set.add_data_point(well, sequence, flow)
            highlow_training_set.add_data_point(well, sequence, flow)
        elif well.startswith('M'):
            medium_training_set.add_data_point(well, sequence, flow)
        elif well.startswith('H'):
            high_training_set.add_data_point(well, sequence, flow)
            highlow_training_set.add_data_point(well, sequence, flow)



    training_set.divide_into_groups()
    low_training_set.divide_into_groups()
    medium_training_set.divide_into_groups()
    high_training_set.divide_into_groups()
    highlow_training_set.divide_into_groups()

    return training_set, low_training_set, medium_training_set, high_training_set, highlow_training_set

def plot_actual_vs_predicted(training_set, label, plot_dir):
    x = []
    y = []
    colours = []

    for data_point in training_set.data_points:
        x.append(data_point.flow)
        y.append(data_point.predicted_flow)
        if data_point.well.startswith('L'):
            colours.append('blue')
        elif data_point.well.startswith('M'):
            colours.append('orange')
        elif data_point.well.startswith('H'):
            colours.append('red')


    plt.xticks([])
    plt.yticks([])
    plt.xlabel('Actual')
    plt.ylabel('Predicted')

    plt.scatter(x, y, marker='o', color=colours)
    plt.savefig(plot_dir + label + '_actual_vs_predicted.png', dpi=500)
    plt.clf()

    correlation, p_val = pearsonr(x, y)
    print(label, correlation, p_val)

def write_actual_vs_predicted(training_set, label, out_dir):
    with open(out_dir + label + '_actual_vs_predicted.txt', 'w') as out_file:
        out_file.write('well\tactual\tpredicted\n')
        for data_point in training_set.data_points:
            out_file.write(f'{data_point.well}\t{data_point.flow}\t{data_point.predicted_flow}\n')




if __name__ == "__main__":
    flow_dir = argv[1]
    translation_dir = argv[2]
    feature_dir = argv[3]
    plot_dir = argv[4]
    out_dir = argv[5]
    tree_dir = argv[6]
    classifier_dir = argv[7]
    

    translation_results = TranslationResults(translation_dir)
    flow_dict = parse_mean_fluorescence(flow_dir)

    training_set, low_training_set, medium_training_set, high_training_set, highlow_training_set = make_training_set(flow_dict, translation_results)
    test_set, low_test_set, medium_test_set, high_test_set, highlow_test_set = make_training_set(flow_dict, translation_results)


 #   training_set.do_crossval(tree_dir)
 #   training_set.write_feature_importances(feature_dir + 'all_feature_importances.txt')

 #   plot_actual_vs_predicted(training_set, 'all', plot_dir)
 #   write_actual_vs_predicted(training_set, 'all', out_dir)

 #   training_set.train_classifier(classifier_dir)


 #   low_training_set.do_crossval(tree_dir)
 #   low_training_set.write_feature_importances(feature_dir + 'low_feature_importances.txt')

 #   plot_actual_vs_predicted(low_training_set, 'low', plot_dir)
 #   write_actual_vs_predicted(low_training_set, 'low', out_dir)

 #   low_training_set.train_classifier(classifier_dir)

 #   low_training_set.test_classifier(high_test_set)
 #   plot_actual_vs_predicted(high_test_set, 'high_from_low', plot_dir)
 #   write_actual_vs_predicted(high_test_set, 'high_from_low', out_dir)

 #   low_training_set.test_classifier(medium_test_set)
 #   plot_actual_vs_predicted(medium_test_set, 'medium_from_low', plot_dir)
 #   write_actual_vs_predicted(medium_test_set, 'medium_from_low', out_dir)



 #   medium_training_set.do_crossval(tree_dir)
 #   medium_training_set.write_feature_importances(feature_dir + 'medium_feature_importances.txt')
 #   plot_actual_vs_predicted(medium_training_set, 'medium', plot_dir)
 #   write_actual_vs_predicted(medium_training_set, 'medium', out_dir)

 #   medium_training_set.train_classifier(classifier_dir)

 #   medium_training_set.test_classifier(high_test_set)
 #   plot_actual_vs_predicted(high_test_set, 'high_from_medium', plot_dir)
 #   write_actual_vs_predicted(high_test_set, 'high_from_medium', out_dir)

 #   medium_training_set.test_classifier(low_test_set)
 #   plot_actual_vs_predicted(low_test_set, 'low_from_medium', plot_dir)
 #   write_actual_vs_predicted(low_test_set, 'low_from_medium', out_dir)



 #   high_training_set.do_crossval(tree_dir)
 #   high_training_set.write_feature_importances(feature_dir + 'high_feature_importances.txt')
 #   plot_actual_vs_predicted(high_training_set, 'high', plot_dir)
 #   write_actual_vs_predicted(high_training_set, 'high', out_dir)

 #   high_training_set.train_classifier(classifier_dir)

 #   high_training_set.test_classifier(medium_test_set)
 #   plot_actual_vs_predicted(medium_test_set, 'medium_from_high', plot_dir)
 #   write_actual_vs_predicted(medium_test_set, 'medium_from_high', out_dir)

 #   high_training_set.test_classifier(low_test_set)
 #   plot_actual_vs_predicted(low_test_set, 'low_from_high', plot_dir)
 #   write_actual_vs_predicted(low_test_set, 'low_from_high', out_dir)



 #   highlow_training_set.do_crossval(tree_dir)

 #   highlow_training_set.write_feature_importances(feature_dir + 'highlow_feature_importances.txt')
 #   plot_actual_vs_predicted(highlow_training_set, 'highlow', plot_dir)
 #   write_actual_vs_predicted(highlow_training_set, 'highlow', out_dir)

 #   highlow_training_set.train_classifier(classifier_dir)

 #   highlow_training_set.test_classifier(medium_test_set)
 #   plot_actual_vs_predicted(medium_test_set, 'medium_from_highlow', plot_dir)
 #   write_actual_vs_predicted(medium_test_set, 'medium_from_highlow', out_dir)










        
        

    
    




    
    
