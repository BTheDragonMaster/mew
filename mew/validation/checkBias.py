import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib import patches
from sys import argv
from mew.parsers import read_sequences
import numpy as np
import os

MmRFP = 'ATGGCNTCNTCNGARGAYGTNATHAARGARTTYATGCGNTTYAARGTNCGNATGGARGGNTCNGTNAAYGGNCAYGARTTYGAaATHGARGGNGARGGNGARGGNCGNCCNTAYGARGGNACNCARACNGCNAARCTNAARGTNACNAARGGNGGNCCNCTNCCNTTcGCcTGGGAYATHCTNTCNCCNCARTTYCARTAYGGNTCNAARGCNTAYGTNAARCAYCCNGCNGAYATHCCNGAYTAYCTgAAaCTNTCNTTYCCNGARGGNTTYAARTGGGARCGNGTNATGAAYTTYGARGAYGGNGGNGTNGTNACNGTNACNCARGAYTCNTCNCTgCARGAYGGNGARTTYATHTAYAARGTNAARCTNCGNGGNACNAAYTTYCCNTCNGAYGGNCCNGTNATGCARAARAARACcATGGGNTGGGARGCNTCNACNGARCGNATGTAYCCNGARGAYGGNGCNCTNAARGGNGARATHAARATGCGNCTNAARCTNAARGAcGGNGGNCAYTAYGAYGCNGARGTNAARACNACNTAYATGGCNAARAARCCNGTNCARCTNCCNGGNGCNTAYAARACNGAcATaAARCTNGAYATHACNTCNCAYAAYGARGAYTAYACNATHGTNGARCARTAYGARCGNGCNGARGGNCGNCAYTCNACNGGNGCN'
HmRFP = 'ATGGCVAGYAGYGAAGAYGTBATYAARGAATTYATGCGYTTYAARGTBCGYATGGAAGGYAGYGTBAAYGGYCAYGAATTYGAAATYGAAGGYGAAGGYGAAGGYCGYCCGTAYGAAGGYACGCARACGGCVAARCTGAARGTBACGAARGGYGGYCCGCTGCCGTTCGCCTGGGAYATYCTGAGYCCGCARTTYCARTAYGGYAGYAARGCVTAYGTBAARCAYCCGGCVGAYATYCCGGAYTAYCTGAAACTGAGYTTYCCGGAAGGYTTYAARTGGGAACGYGTBATGAAYTTYGAAGAYGGYGGYGTBGTBACGGTBACGCARGAYAGYAGYCTGCARGAYGGYGAATTYATYTAYAARGTBAARCTGCGYGGYACGAAYTTYCCGAGYGAYGGYCCGGTBATGCARAARAARACCATGGGYTGGGAAGCVAGYACGGAACGYATGTAYCCGGAAGAYGGYGCVCTGAARGGYGAAATYAARATGCGYCTGAARCTGAARGACGGYGGYCAYTAYGAYGCVGAAGTBAARACGACGTAYATGGCVAARAARCCGGTBCARCTGCCGGGYGCVTAYAARACGGACATAAARCTGGAYATYACGAGYCAYAAYGAAGAYTAYACGATYGTBGAACARTAYGAACGYGCVGAAGGYCGYCAYAGYACGGGYGCV'
LmRFP = 'ATGGCWTCNTCNGAGGACGTMATAAAGGAGTTCATGAGRTTCAAGGTMAGRATGGAGGGRTCNGTMAATGGRCACGAGTTCGAAATAGAGGGRGAGGGRGAGGGRAGRCCHTACGAGGGRACWCAAACWGCWAAGCTHAAGGTMACWAAGGGRGGRCCHCTHCCHTTCGCCTGGGACATACTHTCNCCHCAATTCCAATACGGRTCNAAGGCWTACGTMAAGCACCCHGCWGACATACCHGACTACCTGAAACTHTCNTTCCCHGAGGGRTTCAAGTGGGAGAGRGTMATGAATTTCGAGGACGGRGGRGTMGTMACWGTMACWCAAGACTCNTCNCTGCAAGACGGRGAGTTCATATACAAGGTMAAGCTHAGRGGRACWAATTTCCCHTCNGACGGRCCHGTMATGCAAAAGAAGACCATGGGRTGGGAGGCWTCNACWGAGAGRATGTACCCHGAGGACGGRGCWCTHAAGGGRGAGATAAAGATGAGRCTHAAGCTHAAGGACGGRGGRCACTACGACGCWGAGGTMAAGACWACWTACATGGCWAAGAAGCCHGTMCAACTHCCHGGRGCWTACAAGACWGACATAAAGCTHGACATAACWTCNCACAATGAGGACTACACWATAGTMGAGCAATACGAGAGRGCWGAGGGRAGRCACTCNACWGGRGCW'

MmRFP = MmRFP.upper()
HmRFP = HmRFP.upper()
LmRFP = LmRFP.upper()


def split_data(well_to_sequence):
    well_to_sequence_low = {}
    well_to_sequence_medium_1 = {}
    well_to_sequence_medium_2 = {}
    well_to_sequence_high = {}

    for well, sequence in well_to_sequence.items():

        if well.startswith('L'):
            well_to_sequence_low[well] = sequence

        elif well.startswith('M'):
            well_to_sequence_medium_1[well] = sequence

        elif well.startswith('N'):
            well_to_sequence_medium_2[well] = sequence

        elif well.startswith('H'):
            well_to_sequence_high[well] = sequence

    return (well_to_sequence_low,
            well_to_sequence_medium_1,
            well_to_sequence_medium_2,
            well_to_sequence_high)


DEGENERATE_DICT = {'A': ['A'],
                   'C': ['C'],
                   'T': ['T'],
                   'G': ['G'],
                   'R': ['A', 'G'],
                   'Y': ['C', 'T'],
                   'S': ['G', 'C'],
                   'W': ['A', 'T'],
                   'K': ['G', 'T'],
                   'M': ['A', 'C'],
                   'B': ['C', 'G', 'T'],
                   'D': ['A', 'G', 'T'],
                   'H': ['A', 'C', 'T'],
                   'V': ['A', 'C', 'G'],
                   'N': ['A', 'G', 'C', 'T']}


def get_base_counts(mRFP_sequence, sequences):
    index_to_base_to_count = {}

    for i, base in enumerate(mRFP_sequence):
        index_to_base_to_count[i + 1] = {}
        
        for base_2 in DEGENERATE_DICT[base]:
            index_to_base_to_count[i + 1][base_2] = 0

    wrong_bases = 0
    wrong_sequences = set()
            
    for i, base in enumerate(mRFP_sequence):
        for j, sequence in enumerate(sequences):
            base_2 = sequence[i]
            
            try:
                index_to_base_to_count[i + 1][base_2] += 1
            except KeyError:
                wrong_bases += 1
                wrong_sequences.add(sequence)
                print("Base 2:", base_2)
                print("Actual base:", mRFP_sequence[i])
                print("Base index:", i)
                print("Sequence index:", j)
                pass

    print("Number of removed sequences:")
    print(len(wrong_sequences))
    print("Number of added sequences:")
    print(len(sequences))

    return index_to_base_to_count


def make_bias_dict(mrfp_seq, index_to_base_to_count, seq_nr):
    bias_dict = {}

    for i, base in enumerate(mrfp_seq):
        bases = list(index_to_base_to_count[i + 1].keys())
        
        bias_dict[i + 1] = {}

        expected_frequency = 1.0 / len(bases)

        for base_2 in bases:

            count = index_to_base_to_count[i + 1][base_2]
            observed_frequency = float(count) / seq_nr
            difference = observed_frequency - expected_frequency

            bias_dict[i + 1][base_2] = difference

    return bias_dict

        
def plot_biases(mRFP, sequences, label, out_folder):
    index_to_base_to_count = get_base_counts(mRFP, sequences)
    bias_dict = make_bias_dict(mRFP, index_to_base_to_count, len(sequences))

    plot(mRFP, label, bias_dict, out_folder)


def plot_all_biases(mrfps, sequence_sets, labels, out_folder):
    bias_dicts = []
    for i, mrfp in enumerate(mrfps):
        sequences = sequence_sets[i]
        index_to_base_to_count = get_base_counts(mrfp, sequences)
        bias_dict = make_bias_dict(mrfp, index_to_base_to_count, len(sequences))
        bias_dicts.append(bias_dict)

    plot_all(mrfps, labels, bias_dicts, out_folder)



# Take negative and positive data apart and cumulate
def get_cumulated_array(data, **kwargs):
    cum = data.clip(**kwargs)
    cum = np.cumsum(cum, axis=0)
    d = np.zeros(np.shape(data))
    d[1:] = cum[:-1]
    return d

def make_patches(colours):
    legend_patches = []
    for colour in colours:
        patch = patches.Patch(edgecolor='black', facecolor=colour)
        legend_patches.append(patch)

    return legend_patches


def plot(mRFP, label, bias_dict, out_folder):

    a = []
    c = []
    g = []
    t = []

    for position, base_to_difference in bias_dict.items():
        for base in ['A', 'C', 'G', 'T']:
            if base not in base_to_difference:
                if base == 'A':
                    a.append(0.0)
                elif base == 'C':
                    c.append(0.0)
                elif base == 'G':
                    g.append(0.0)
                elif base == 'T':
                    t.append(0.0)
            else:
                difference = base_to_difference[base]

                if base == 'A':
                    a.append(difference)
                elif base == 'C':
                    c.append(difference)
                elif base == 'G':
                    g.append(difference)
                elif base == 'T':
                    t.append(difference)

    data = np.array([a, c, g, t])

    data_shape = np.shape(data)

    cumulated_data = get_cumulated_array(data, min=0)
    cumulated_data_neg = get_cumulated_array(data, max=0)

    # Re-merge negative and positive data.
    row_mask = (data<0)
    cumulated_data[row_mask] = cumulated_data_neg[row_mask]
    data_stack = cumulated_data

    cols = ["red", "orange", "green", "blue"]

    plt.rcParams.update({'font.size': 8,
                         'figure.figsize': (85, 5)})

    fig = plt.figure()
    ax = plt.subplot(111)

    print(data_shape)

    for i in range(4):
        ax.bar(np.arange(data_shape[1]), data[i], width=1, bottom=data_stack[i], color=cols[i])

    ax.set_title(f"Difference between observed and expected frequency for bases in {label}.", fontsize=28)
    ax.set_xticks(np.arange(len(mRFP)))
    ax.set_xticklabels(list(mRFP))

    ax.legend(make_patches(cols), ['A', 'C', 'G', 'T'])
    ax.set_ylabel('observed - expected frequency', fontsize=14)
    ax.set_xlim(-2, len(mRFP) + 2)

    plt.savefig(os.path.join(out_folder, label + '.svg'))
    plt.clf()


def plot_all(mrfps, labels, bias_dicts, out_folder):
    plt.rcParams.update({'font.size': 8,
                         'figure.figsize': (85, 20)})

    fig, axes = plt.subplots(3, 1)

    fig.suptitle(f"Difference between observed and expected base frequencies.", fontsize=20)

    for i, ax in enumerate(axes):
        label = labels[i]
        bias_dict = bias_dicts[i]
        mRFP = mrfps[i]
        a = []
        c = []
        g = []
        t = []

        for position, base_to_difference in bias_dict.items():
            for base in ['A', 'C', 'G', 'T']:
                if base not in base_to_difference:
                    if base == 'A':
                        a.append(0.0)
                    elif base == 'C':
                        c.append(0.0)
                    elif base == 'G':
                        g.append(0.0)
                    elif base == 'T':
                        t.append(0.0)
                else:
                    difference = base_to_difference[base]

                    if base == 'A':
                        a.append(difference)
                    elif base == 'C':
                        c.append(difference)
                    elif base == 'G':
                        g.append(difference)
                    elif base == 'T':
                        t.append(difference)

        data = np.array([a, c, g, t])

        data_shape = np.shape(data)

        cumulated_data = get_cumulated_array(data, min=0)
        cumulated_data_neg = get_cumulated_array(data, max=0)

        # Re-merge negative and positive data.
        row_mask = (data < 0)
        cumulated_data[row_mask] = cumulated_data_neg[row_mask]
        data_stack = cumulated_data

        cols = ["red", "orange", "green", "blue"]



        print(data_shape)

        for i in range(4):
            ax.bar(np.arange(data_shape[1]), data[i], width=1, bottom=data_stack[i], color=cols[i])

        ax.set_title(label, fontsize=16)
        ax.set_xticks(np.arange(len(mRFP)))
        ax.set_xticklabels(list(mRFP))

        ax.legend(make_patches(cols), ['A', 'C', 'G', 'T'])
        ax.set_ylabel('observed - expected frequency', fontsize=14)
        ax.set_xlim(-2, len(mRFP) + 2)

    plt.savefig(os.path.join(out_folder, 'biases.svg'))
    plt.clf()
        
        

if __name__ == "__main__":
    well_to_sequence = read_sequences(argv[1])
    out_folder = argv[2]
    well_to_low, well_to_medium_1, well_to_medium_2, well_to_high = split_data(well_to_sequence)

    low_sequences = list(well_to_low.values())
    medium_1_sequences = list(well_to_medium_1.values())
    medium_2_sequences = list(well_to_medium_2.values())
    medium_sequences = medium_1_sequences + medium_2_sequences
    high_sequences = list(well_to_high.values())

    plot_biases(MmRFP, medium_1_sequences, 'CAI medium 1', out_folder)
    plot_biases(MmRFP, medium_2_sequences, 'CAI medium 2', out_folder)
    plot_biases(MmRFP, medium_sequences, 'CAI medium', out_folder)
    plot_biases(LmRFP, low_sequences, 'CAI low', out_folder)
    plot_biases(HmRFP, high_sequences, 'CAI high', out_folder)

    mrfps = [LmRFP, MmRFP, HmRFP]
    sequence_sets = [low_sequences, medium_sequences, high_sequences]
    labels = ["Low CAI", "Medium CAI", "High CAI"]
    plot_all_biases(mrfps, sequence_sets, labels, out_folder)