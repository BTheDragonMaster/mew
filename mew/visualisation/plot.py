#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt, patches, rcParams
rcParams['font.family'] = 'palatino'

from scipy.stats import pearsonr, spearmanr
import os
from sys import argv

from mew.parsers import read_feature_importances_from_dir, read_correlations_from_dir


def plot_indicator_bar(smallest_x, max_x, ax):
    fragments = ["ATGGGCCCAAGTTCACTTAA", "AAAGGAGAT", "CAACAATGAAAGCAATTTTCGTACTGAAACATCTTAATCATGC", "AGGGGAGGGT", "TTCTA",
                 "ATGGCVAGYAGYGAAGAYGTBATYAARGAATTYATGCGYTTYAARGTBCGYATGGAAGGYAGYGTBAAYGGYCAYGAATTYGAAATYGAAGGYGAAGGYGAAGGYCGYCCGTAYGAAGGYACGCARACGGCVAARCTGAARGTBACGAARGGYGGYCCGCTGCCGTTCGCCTGGGAYATYCTGAGYCCGCARTTYCARTAYGGYAGYAARGCVTAYGTBAARCAYCCGGCVGAYATYCCGGAYTAYCTGAAACTGAGYTTYCCGGAAGGYTTYAARTGGGAACGYGTBATGAAYTTYGAAGAYGGYGGYGTBGTBACGGTBACGCARGAYAGYAGYCTGCARGAYGGYGAATTYATYTAYAARGTBAARCTGCGYGGYACGAAYTTYCCGAGYGAYGGYCCGGTBATGCARAARAARACCATGGGYTGGGAAGCVAGYACGGAACGYATGTAYCCGGAAGAYGGYGCVCTGAARGGYGAAATYAARATGCGYCTGAARCTGAARGACGGYGGYCAYTAYGAYGCVGAAGTBAARACGACGTAYATGGCVAARAARCCGGTBCARCTGCCGGGYGCVTAYAARACGGACATAAARCTGGAYATYACGAGYCAYAAYGAAGAYTAYACGATYGTBGAACARTAYGAACGYGCVGAAGGYCGYCAYAGYACGGGYGCVTGA",
                 "TGCCGACTCAGTTGCTGCTTCTACTGGGCG", "CCCCGCTTCGGCGGGGTTTTTTT"]
    labels = ["", "RBS1", "BCD", "RBS2", "", "CDS", "UTR", "TER"]
    colors = ["#A0B2BB".lower(), "#F38424".lower(), "#A0B2BB".lower(), "#F38424".lower(), "#A0B2BB".lower(), "#0497DB".lower(), "#A0B2BB".lower(), "#FF1F5D".lower()]

    bcd_length = len(fragments[0]) + len(fragments[1]) + len(fragments[2]) + len(fragments[3]) + len(fragments[4])

    labels_and_ranges = []

    total_length = -bcd_length

    for i, fragment in enumerate(fragments):

        print(total_length)
        label = labels[i]
        fragment_range = (total_length, total_length + len(fragment))
        labels_and_ranges.append((label, fragment_range))
        total_length += len(fragment)

    for i, label_and_range in enumerate(labels_and_ranges):
        label, fragment_range = label_and_range
        color = colors[i]
        x_min = max(smallest_x, fragment_range[0])
        x_max = min(fragment_range[1], max_x)

        if not x_max <= smallest_x and not x_min >= max_x:
            text_x = int((x_min + x_max) / 2)

            if x_max - x_min < 8:
                label = ''
                rotation = 'vertical'
            elif x_max - x_min < 20:
                rotation = 'vertical'
            else:
                rotation = 'horizontal'

            if rotation == 'horizontal':

                ax.text(text_x, 0.77, label, fontdict={'family': "palatino", 'size': 18, 'weight': 'bold'}, horizontalalignment='center',
                        verticalalignment='center', rotation=rotation)

            y_min = 0.65
            y_max = 0.95

            if i == 0:
                x_min += 1.0
            elif i == 7:
                x_max -= 1.0
            ax.axvspan(x_min, x_max, y_min, y_max, facecolor=color, fill=True, edgecolor='black', label=label, linewidth=2.0)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)


def plot_multiple_correlations(window_dirs, window_sizes, labels, out_svg, title, correlation_type="pearson"):
    figure, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 2]}, figsize=(27, 5))


    plt.subplots_adjust(left=0.125, right=0.9, bottom=0.15, top=0.9, wspace=0.5, hspace=0.5)
    colours = ["#FC0003".lower(), "#0096DB".lower(), "#F17B14".lower(), "#823D8B".lower(), "pink", "black", "grey", "cyan"]
    max_xs = []
    min_xs = []
    max_vals = []
    min_vals = []
    ax = axes[0]


    for i, window_dir in enumerate(window_dirs):
        window_size = window_sizes[i]
        label = labels[i]

        colour = colours[i % len(colours)]

        window_to_correlations = read_correlations_from_dir(window_dir)

        x_axis = sorted(list(window_to_correlations.keys()))
        print(window_size)

        x = [y + (int(0.5 * window_size)) for y in x_axis]
        print(x[0], x[-1])

        pearson = []
        spearman = []

        for window in x_axis:
            pearson.append(window_to_correlations[window]['pearson'])
            spearman.append(window_to_correlations[window]['spearman'])

        if correlation_type == "pearson":
            subplot_importances(colour, x, pearson, label, ax)
            max_value = max(pearson)
            min_value = min(pearson)
        elif correlation_type == "spearman":
            subplot_importances(colour, x, spearman, label, ax)
            max_value = max(spearman)
            min_value = min(spearman)

        max_vals.append(max_value)
        min_vals.append(min_value)

        max_xs.append(x[-1])
        min_xs.append(x[0])




    max_value = max(max_vals)
    min_value = min(min_vals)
    max_x = max(max_xs)
    min_x = min(min_xs)

    indicator_bar = axes[1]
    indicator_bar.get_xaxis().set_visible(False)
    indicator_bar.get_yaxis().set_visible(False)
    plot_indicator_bar(min_x, max_x, indicator_bar)

    y_lim_top = max(0, max_value + 0.1 * max_value)
    y_lim_bottom = min(0, min_value + 0.1 * min_value)

    ax.set_ylim(y_lim_bottom, y_lim_top)
    ax.set_xlim(min_x, max_x)

    indicator_bar.set_ylim(0, 1)
    indicator_bar.set_xlim(min_x, max_x)

    ax.set_xlabel("Base position", fontdict={'size': 20})
    ax.set_ylabel('Correlation', fontdict={'size': 16})

    plt.gcf().subplots_adjust(bottom=0.15)

    figure.suptitle(title, fontsize=22)
    ax.legend(loc='upper right')
    plt.savefig(out_svg)
    plt.clf()
    plt.close()


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
            colours.append('#019ADF'.lower())
        elif data_point.well.startswith('M'):
            colours.append('#FFC622'.lower())
        elif data_point.well.startswith('N'):
            colours.append('#FFC622'.lower())
        elif data_point.well.startswith('H'):
            colours.append('#FF1F5D'.lower())

    pearson_correlation, p_val = pearsonr(x, y)
    spearman_correlation, sp_val = spearmanr(x, y)


    plt.xlabel('Actual', fontsize=16)
    plt.ylabel('Predicted', fontsize=16)
    plt.title(f'Pearson: {pearson_correlation:.3f} (p={p_val:.3f})\nSpearman: {spearman_correlation:.3f} (p={sp_val:.3f})', fontsize=16)

    plt.scatter(x, y, marker='o', color=colours)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
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

    plot_importances('#FF1656'.lower(), x, A, 'A')
    plot_importances('#0194DD'.lower(), x, T, 'T')
    plot_importances('#06B658'.lower(), x, G, 'G')
    plot_importances('#FFDA6B'.lower(), x, C, 'C')

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

    plot_importances('#FF1656'.lower(), x_onehot, A, 'A')
    plot_importances('#0194DD'.lower(), x_onehot, T, 'T')
    plot_importances('#06B658'.lower(), x_onehot, G, 'G')
    plot_importances('#FFDA6B'.lower(), x_onehot, C, 'C')
    plot_importances('grey', x_bpp, bpp, 'mRNA')

    plt.legend(loc='upper right')

    out_svg = os.path.join(out_dir, 'feature_importances.svg')
    plt.savefig(out_svg)
    plt.clf()
    plt.close()


def subplot_feature_importances_onehot_base(fi_dir, encoding, label, ax):
    ax.set_title(label, fontsize=20)

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

    ax.set_xlabel(x_label, fontsize=18)
    ax.set_ylabel('Feature importance', fontsize=16)

    plt.gcf().subplots_adjust(bottom=0.15)

    subplot_importances('#FF1656'.lower(), x, A, 'A', ax)
    subplot_importances('#0194DD'.lower(), x, T, 'T', ax)
    subplot_importances('#06B658'.lower(), x, G, 'G', ax)
    subplot_importances('#FFDA6B'.lower(), x, C, 'C', ax)

    ax.legend(loc='upper right')


def subplot_feature_importances_bpp_onehot_base_third(fi_dir, label, ax):
    ax.set_title(label, fontsize=20)

    base_to_type_to_importance = read_feature_importances_from_dir(fi_dir, 'rna-bppm-onehot-third')

    x_onehot, x_bpp, A, C, G, T, bpp = get_importance_vectors_bpp_onehot(base_to_type_to_importance)

    max_value = max(A + C + G + T + bpp)
    min_value = min(A + C + G + T + bpp)

    y_lim_top = max(0, max_value + 0.1 * max_value)
    y_lim_bottom = min(0, min_value + 0.1 * min_value)

    x_axis = sorted(x_bpp)

    ax.set_ylim(y_lim_bottom, y_lim_top)
    ax.set_xlim(x_axis[0], x_axis[-1])

    ax.set_xlabel('Base position', fontsize=18)
    ax.set_ylabel('Feature importance', fontsize=16)

    subplot_importances('#FF1656'.lower(), x_onehot, A, 'A', ax)
    subplot_importances('#0194DD'.lower(), x_onehot, T, 'T', ax)
    subplot_importances('#06B658'.lower(), x_onehot, G, 'G', ax)
    subplot_importances('#FFDA6B'.lower(), x_onehot, C, 'C', ax)
    subplot_importances('grey', x_bpp, bpp, 'mRNA', ax)

    ax.legend(loc='upper right')


def subplot_feature_importances_bpp(fi_dir, label, ax):
    ax.set_title(label, fontsize=20)

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

    ax.set_xlabel('Base position', fontsize=18)
    ax.set_ylabel('Feature importance', fontsize=16)

    ax.plot(x_axis, y_axis, color='black')


def plot_multiple(feature_importances_files, encodings, labels, figure_title, out_svg):

    fig, axes = plt.subplots(len(feature_importances_files) + 1, 1, gridspec_kw={'height_ratios': [4.0] * len(feature_importances_files) + [2]},
                             figsize=(27, 3.0 * len(feature_importances_files) + 4))

    plt.subplots_adjust(left=0.125, right=0.9, bottom=0.15, top=0.9, wspace=0.5, hspace=0.5)

    fig.suptitle(figure_title, fontsize=14)

    min_x = 10000
    max_x = -1000

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

        x_axis = ax.get_xaxis()
        min_x_candidate, max_x_candidate = x_axis.get_data_interval()
        if min_x_candidate < min_x:
            min_x = min_x_candidate
        if max_x_candidate > max_x:
            max_x = max_x_candidate

    indicator_bar = axes[-1]

    indicator_bar.get_xaxis().set_visible(False)
    indicator_bar.get_yaxis().set_visible(False)
    plot_indicator_bar(min_x, max_x, indicator_bar)

    indicator_bar.set_ylim(0, 1)
    indicator_bar.set_xlim(min_x, max_x)

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
