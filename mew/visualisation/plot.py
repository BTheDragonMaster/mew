#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt

from scipy.stats import pearsonr, spearmanr
import os
from sys import argv

from mew.parsers import read_feature_importances_from_dir, read_correlations_from_dir


def plot_correlation_vs_window(window_dir, out_svg, window_size):
    window_to_correlations = read_correlations_from_dir(window_dir)

    plt.figure(figsize=(27, 3))

    x_axis = sorted(list(window_to_correlations.keys()))
    x = [y + (int(0.5 * window_size)) for y in x_axis]

    pearson = []
    spearman = []

    for window in x_axis:
        pearson.append(window_to_correlations[window]['pearson'])
        spearman.append(window_to_correlations[window]['spearman'])

    plot_importances('black', x, pearson, 'pearson')
    plot_importances('grey', x, spearman, 'spearman')


    max_value = max(spearman + pearson)
    min_value = min(spearman + pearson)

    y_lim_top = max(0, max_value + 0.1 * max_value)
    y_lim_bottom = min(0, min_value + 0.1 * min_value)

    plt.axis([x[0], x[-1], y_lim_bottom, y_lim_top])

    plt.xlabel("Base position")
    plt.ylabel('Correlation')

    plt.gcf().subplots_adjust(bottom=0.15)

    plt.legend(loc='upper right')
    plt.savefig(out_svg)
    plt.clf()
    plt.close()


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
            colours.append('orange')
        elif data_point.well.startswith('H'):
            colours.append('red')

    pearson_correlation, p_val = pearsonr(x, y)
    spearman_correlation, sp_val = spearmanr(x, y)


    plt.xlabel('Actual')
    plt.ylabel('Predicted')
    plt.title(f'Pearson: {pearson_correlation:.3f} (p={p_val:.3f}), Spearman: {spearman_correlation:.3f} (p={sp_val:.3f})')

    plt.scatter(x, y, marker='o', color=colours)
    plt.savefig(os.path.join(plot_dir, label + '_actual_vs_predicted.svg'))
    plt.clf()
    plt.close()


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

def get_importance_vectors_bpp_onehot(average_dict):

    x_onehot = []
    x_bpp = []

    a_importance_vector = []
    c_importance_vector = []
    g_importance_vector = []
    t_importance_vector = []
    bpp_importance_vector = []

    for base, basetype_to_importance in average_dict.items():
        x_bpp.append(base)
        if 'A' in basetype_to_importance:
            x_onehot.append(base)

        for basetype, importance in basetype_to_importance.items():

            if basetype == 'A':
                a_importance_vector.append(importance)
            elif basetype == 'C':
                c_importance_vector.append(importance)
            elif basetype == 'G':
                g_importance_vector.append(importance)
            elif basetype == 'T':
                t_importance_vector.append(importance)
            elif basetype == 'bpp':
                bpp_importance_vector.append(importance)

    return x_onehot, x_bpp, a_importance_vector, c_importance_vector, g_importance_vector, t_importance_vector, bpp_importance_vector


def plot_importances(colour, x, y, label):

    plt.plot(x, y, color=colour, label=label)


def subplot_importances(colour, x, y, label, ax):
    ax.plot(x, y, color=colour, label=label)


def plot_feature_importances_onehot_base(fi_dir, encoding, out_dir):

    base_to_type_to_importance = read_feature_importances_from_dir(fi_dir, encoding)

    plt.figure(figsize=(27, 3))

    x, A, C, G, T = get_importance_vectors(base_to_type_to_importance, encoding)

    max_value = max(A + C + G + T)
    min_value = min(A + C + G + T)

    y_lim_top = max(0, max_value + 0.1 * max_value)
    y_lim_bottom = min(0, min_value + 0.1 * min_value)

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
    plt.close()


def plot_feature_importances_bpp(fi_dir, out_dir):

    base_to_importance = read_feature_importances_from_dir(fi_dir, 'rna-bppm-totals')

    plt.figure(figsize=(27, 3))

    x_axis = sorted(base_to_importance.keys())
    y_axis = []

    for base in x_axis:
        y_axis.append(base_to_importance[base])

    max_value = max(y_axis)
    min_value = min(y_axis)

    y_lim_top = max(0, max_value + 0.1 * max_value)
    y_lim_bottom = min(0, min_value + 0.1 * min_value)

    out_svg = os.path.join(out_dir, 'feature_importances.svg')

    plt.axis([x_axis[0], x_axis[-1], y_lim_bottom, y_lim_top])

    plt.xlabel('Base position')
    plt.ylabel('Feature importance')

    plt.gcf().subplots_adjust(bottom=0.15)

    plt.plot(x_axis, y_axis, color='grey')

    plt.savefig(out_svg)
    plt.clf()
    plt.close()


def plot_feature_importances_bpp_onehot_base_third(fi_dir, out_dir):

    base_to_type_to_importance = read_feature_importances_from_dir(fi_dir, 'rna-bppm-onehot-third')

    plt.figure(figsize=(27, 3))

    x_onehot, x_bpp, A, C, G, T, bpp = get_importance_vectors_bpp_onehot(base_to_type_to_importance)

    max_value = max(A + C + G + T + bpp)
    min_value = min(A + C + G + T + bpp)

    y_lim_top = max(0, max_value + 0.1 * max_value)
    y_lim_bottom = min(0, min_value + 0.1 * min_value)

    x_axis = sorted(x_bpp)

    plt.axis([x_axis[0], x_axis[-1], y_lim_bottom, y_lim_top])

    plt.xlabel("Base position")
    plt.ylabel('Feature importance')

    plt.gcf().subplots_adjust(bottom=0.15)

    plot_importances('red', x_onehot, A, 'A')
    plot_importances('blue', x_onehot, T, 'T')
    plot_importances('green', x_onehot, G, 'G')
    plot_importances('orange', x_onehot, C, 'C')
    plot_importances('grey', x_bpp, bpp, 'mRNA')

    plt.legend(loc='upper right')

    out_svg = os.path.join(out_dir, 'feature_importances.svg')
    plt.savefig(out_svg)
    plt.clf()
    plt.close()


def subplot_feature_importances_onehot_base(fi_dir, encoding, label, ax):
    ax.set_title(label, fontsize=12)

    base_to_type_to_importance = read_feature_importances_from_dir(fi_dir, encoding)

    x, A, C, G, T = get_importance_vectors(base_to_type_to_importance, encoding)

    max_value = max(A + C + G + T)
    min_value = min(A + C + G + T)

    y_lim_top = max(0, max_value + 0.1 * max_value)
    y_lim_bottom = min(0, min_value + 0.1 * min_value)

    ax.set_ylim(y_lim_bottom, y_lim_top)
    ax.set_xlim(1, max(x))

    if encoding == 'one-hot-base':
        x_label = 'Base position'
    elif encoding == 'one-hot-third-base':
        x_label = 'Codon position'

    ax.set_xlabel(x_label)
    ax.set_ylabel('Feature importance')

    plt.gcf().subplots_adjust(bottom=0.15)

    subplot_importances('red', x, A, 'A', ax)
    subplot_importances('blue', x, T, 'T', ax)
    subplot_importances('green', x, G, 'G', ax)
    subplot_importances('orange', x, C, 'C', ax)

    ax.legend(loc='upper right')


def subplot_feature_importances_bpp_onehot_base_third(fi_dir, label, ax):
    ax.set_title(label, fontsize=12)

    base_to_type_to_importance = read_feature_importances_from_dir(fi_dir, 'rna-bppm-onehot-third')

    x_onehot, x_bpp, A, C, G, T, bpp = get_importance_vectors_bpp_onehot(base_to_type_to_importance)

    max_value = max(A + C + G + T + bpp)
    min_value = min(A + C + G + T + bpp)

    y_lim_top = max(0, max_value + 0.1 * max_value)
    y_lim_bottom = min(0, min_value + 0.1 * min_value)

    x_axis = sorted(x_bpp)

    ax.set_ylim(y_lim_bottom, y_lim_top)
    ax.set_xlim(x_axis[0], x_axis[-1])

    ax.set_xlabel('Base position')
    ax.set_ylabel('Feature importance')

    subplot_importances('red', x_onehot, A, 'A', ax)
    subplot_importances('blue', x_onehot, T, 'T', ax)
    subplot_importances('green', x_onehot, G, 'G', ax)
    subplot_importances('orange', x_onehot, C, 'C', ax)
    subplot_importances('grey', x_bpp, bpp, 'mRNA', ax)

    ax.legend(loc='upper right')


def subplot_feature_importances_bpp(fi_dir, label, ax):
    ax.set_title(label, fontsize=12)

    base_to_importance = read_feature_importances_from_dir(fi_dir, 'rna-bppm-totals')

    x_axis = sorted(base_to_importance.keys())
    y_axis = []

    for base in x_axis:
        y_axis.append(base_to_importance[base])

    max_value = max(y_axis)
    min_value = min(y_axis)

    y_lim_top = max(0, max_value + 0.1 * max_value)
    y_lim_bottom = min(0, min_value + 0.1 * min_value)

    ax.set_ylim(y_lim_bottom, y_lim_top)
    ax.set_xlim(x_axis[0], x_axis[-1])

    ax.set_xlabel('Base position')
    ax.set_ylabel('Feature importance')

    ax.plot(x_axis, y_axis, color='grey')


def plot_multiple(feature_importances_files, encodings, labels, figure_title, out_svg):

    fig, axes = plt.subplots(len(feature_importances_files), 1, figsize=(27, 3.5 * len(feature_importances_files)))

    plt.subplots_adjust(left=0.125, right=0.9, bottom=0.15, top=0.9, wspace=0.5, hspace=0.5)

    fig.suptitle(figure_title, fontsize=14)

    for i, fi_dir in enumerate(feature_importances_files):
        encoding = encodings[i]
        ax = axes[i]
        label = labels[i]

        if encoding == 'rna-bppm-totals':
            subplot_feature_importances_bpp(fi_dir, label, ax)
        elif encoding == 'one-hot-base' or encoding == 'one-hot-third-base':
            subplot_feature_importances_onehot_base(fi_dir, encoding, label, ax)
        elif encoding == 'rna-bppm-onehot-third':
            subplot_feature_importances_bpp_onehot_base_third(fi_dir, label, ax)

    plt.savefig(out_svg)
    plt.clf()
    plt.close()


def plot_feature_importances(fi_dir, encoding, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if encoding == 'rna-bppm-totals':
        plot_feature_importances_bpp(fi_dir, out_dir)
    elif encoding == 'one-hot-base' or encoding == 'one-hot-third-base':
        plot_feature_importances_onehot_base(fi_dir, encoding, out_dir)
    elif encoding == 'rna-bppm-onehot-third':
        plot_feature_importances_bpp_onehot_base_third(fi_dir, out_dir)


if __name__ == "__main__":
    plot_feature_importances(argv[1], 'rna-bppm-onehot-third', os.path.join(os.getcwd(), 'test'))
