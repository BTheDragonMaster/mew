#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt

from scipy.stats import pearsonr, spearmanr
import os
from sys import argv

from mew.parsers import read_feature_importances_from_dir

def plot_actual_vs_predicted(training_set, plot_dir):
    label = training_set.label

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
        elif data_point.well.startswith('N'):
            colours.append('purple')
        elif data_point.well.startswith('H'):
            colours.append('red')

    pearson_correlation, p_val = pearsonr(x, y)
    spearman_correlation, sp_val = spearmanr(x, y)


    #plt.xticks([])
    #plt.yticks([])
    plt.xlabel('Actual')
    plt.ylabel('Predicted')
    plt.title(f'Pearson: {pearson_correlation:.3f} (p={p_val:.3f}), Spearman: {spearman_correlation:.3f} (p={sp_val:.3f})')

   # plt.axis([0, 100000, 0, 100000])

    plt.scatter(x, y, marker='o', color=colours)
    plt.savefig(os.path.join(plot_dir, label + '_actual_vs_predicted.svg'))
    plt.clf()


def get_importance_vectors(average_dict, encoding):

    x = []
    a_importance_vector = []
    c_importance_vector = []
    g_importance_vector = []
    t_importance_vector = []

    for base, basetype_to_importance in average_dict.items():
        if encoding == 'one-hot-base':
            x.append(base)
        elif encoding == 'one-hot-third-base':
            x.append(int(base / 3))

        for basetype, importance in basetype_to_importance.items():
            if basetype == 'A':
                a_importance_vector.append(importance)
            elif basetype == 'C':
                c_importance_vector.append(importance)
            elif basetype == 'G':
                g_importance_vector.append(importance)
            elif basetype == 'T':
                t_importance_vector.append(importance)

    return x, a_importance_vector, c_importance_vector, g_importance_vector, t_importance_vector


def plot_importances(colour, x, y, label):

    plt.plot(x, y, color=colour, label=label)


def plot_feature_importances_onehot_base(fi_dir, encoding, out_dir):

    base_to_type_to_importance = read_feature_importances_from_dir(fi_dir, encoding)

    plt.figure(figsize=(27, 3))

    x, A, C, G, T = get_importance_vectors(base_to_type_to_importance, encoding)

    max_value = max(A + C + G + T)
    min_value = min(A + C + G + T)

    y_lim_top = max_value + 0.1 * max_value
    y_lim_bottom = min_value + 0.1 * min_value

    out_svg = os.path.join(out_dir, 'feature_importances.svg')

    plt.axis([1, max(x), y_lim_bottom, y_lim_top])

    if encoding == 'one-hot-base':
        x_label = 'Base position'
    elif encoding == 'one-hot-third-base':
        x_label = 'Codon position'

    plt.xlabel(x_label)
    plt.ylabel('Feature importance')

    plt.gcf().subplots_adjust(bottom=0.15)

    plot_importances('red', x, A, 'A')
    plot_importances('blue', x, T, 'T')
    plot_importances('green', x, G, 'G')
    plot_importances('orange', x, C, 'C')

    plt.legend(loc='upper right')
    plt.savefig(out_svg)
    plt.clf()


def plot_feature_importances_bpp(fi_dir, out_dir):

    base_to_importance = read_feature_importances_from_dir(fi_dir, 'rna-bppm-totals')

    plt.figure(figsize=(27, 3))

    x_axis = sorted(base_to_importance.keys())
    y_axis = []

    for base in x_axis:
        y_axis.append(base_to_importance[base])

    max_value = max(y_axis)
    min_value = min(y_axis)

    y_lim_top = max_value + 0.1 * max_value
    y_lim_bottom = min_value + 0.1 * min_value

    out_svg = os.path.join(out_dir, 'feature_importances.svg')

    plt.axis([x_axis[0], x_axis[-1], y_lim_bottom, y_lim_top])

    plt.xlabel('Base position')
    plt.ylabel('Feature importance')

    plt.gcf().subplots_adjust(bottom=0.15)

    plt.plot(x_axis, y_axis, color='grey')

    plt.savefig(out_svg)
    plt.clf()


def plot_feature_importances(fi_dir, encoding, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if encoding == 'rna-bppm-totals':
        plot_feature_importances_bpp(fi_dir, out_dir)
    elif encoding == 'one-hot-base' or encoding == 'one-hot-third-base':
        plot_feature_importances_onehot_base(fi_dir, encoding, out_dir)


if __name__ == "__main__":
    plot_feature_importances(argv[1], 'rna-bppm-totals', os.path.join(os.getcwd(), 'test'))
