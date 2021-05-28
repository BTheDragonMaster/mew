#!/usr/bin/env python

from sklearn.ensemble import RandomForestClassifier


from sys import argv, path
import os
from joblib import dump, load
from typing import Any, Dict, List, Tuple, Optional
import random
from pprint import pprint
import argparse
from one_hot_encoding import *

def define_arguments():
    parser = argparse.ArgumentParser(description = "Script to train a random \
forest.")
    parser.add_argument("--fasta", required = True, type = str, help = "Fasta \
directory.")
    parser.add_argument("--trees", type = int, default = 1000, help = "Number \
of trees in random forest.")
    parser.add_argument("--depth", type = int, default = 10, help = "Maximum \
tree depth.")
    parser.add_argument("--testtrain", default = None, help = "Directory of \
file identifying data.py point IDs as test or training data.py points.")
    parser.add_argument("--min_samples", default = 1, type = int,
                        help = "Minimum number of samples per leaf.")
    return parser

def get_sequence_and_flow_tuples(flow_dict, translation_results):
    seqflow_tuples = []
    wells = list(flow_dict.keys())
    wells.sort()

    for well in wells:
        flow = flow_dict[well]
        sequence = CodingSequence(translation_results.results[well].orf)
        seq_vector = sequence.onehot_encode_base_12()
        seqflow_tuples.append((seq_vector, flow))

    return seqflow_tuples
        
                                  

def train_RF(training_features, training_response, n_estimators,
             max_features = 'auto', max_depth = None, random_state = None,
             oob_score = True, bootstrap = True, min_samples_leaf = 1):
    forest = RandomForestClassifier(n_estimators = n_estimators,
                                    oob_score = oob_score,
                                    max_features = max_features,
                                    max_depth = max_depth,
                                    random_state = random_state,
                                    class_weight = "balanced_subsample",
                                    bootstrap = bootstrap,
                                    n_jobs = -1,
                                    min_samples_leaf = min_samples_leaf)
    forest.fit(training_features, training_response)
    return forest

def store_model(classifier, name):
    dump(classifier, '%s.joblib' % name)


if __name__ == "__main__":
    current_folder = os.path.dirname(os.path.abspath(__file__))
    parent_folder = os.path.dirname(current_folder)

#    properties = parse_15_properties("%s/Data_preparation/15_aa_properties.txt"
 #                                    % parent_folder)
    properties = parse_4_properties("%s/Data_preparation/aa_properties.txt"
                                    % parent_folder)
 #   features, response, seq_IDs = get_feature_matrix(argv[1], properties, 9)

    features_test, response_test, features_training, response_training = \
                   split_test_training_from_file(argv[1], properties, argv[2])


    chuck_dict = make_chuck_dict()
    
    new_response_training = []
    new_features_training = []
    for i, response in enumerate(response_training):
        if not chuck_dict[response] == 'reject':
            new_response_training.append(chuck_dict[response])
            new_features_training.append(features_training[i])
    response_training = new_response_training
    features_training = new_features_training

    new_response_test = []
    new_features_test = []
    for i, response in enumerate(response_test):
        if not chuck_dict[response] == 'reject':
            new_response_test.append(chuck_dict[response])
            new_features_test.append(features_test[i])
            
    response_test = new_response_test
    features_test = new_features_test
    
    forest = train_RF(features_training, response_training, 1000,
                      0.1, 9, 12)
    store_model(forest, 'RF_1000trees_0.1_depth9_4ft_groups')
    print(forest.oob_score_)
    print(test(forest, features_test, response_test))
    accuracy_dict = test_per_class(forest, features_test, response_test)

    print_accuracy_dict(accuracy_dict)

    

    
    
        #print("Forest size:", forest_size)
        #print("Out of bag:", forest.oob_score_)
        #print("Feature importances:", forest.feature_importances_)
        #print("Test:", test(forest, clean_test_features, clean_test_response))
        #print("Function:", forest.oob_decision_function_)
        

    
