#!/usr/bin/env python
import os


BCD = "ATGGGCCCAAGTTCACTTAAAAAGGAGATCAACAATGAAAGCAATTTTCGTACTGAAACATCTTAATCATGCAGGGGAGGGTTTCTA"
TERMINATOR = "tgccgactcagttgctgcttctactgggcgccccgcttcggcggggttttttt".upper()

def parse_crrfp_data(data_file, out_folder):
    flow_dict = {}
    sequence_dict = {}
    cai_dict = {}

    with open(data_file, 'r') as data:
        data.readline()
        for line in data:
            line = line.strip()
            line_data = line.split(',')

            dataset = line_data[0]
            plate_id = line_data[1]
            flow = line_data[8]
            sequence = line_data[11]
            cai = line_data[12]

            dataset_id = dataset.split('-')[1]
            if dataset_id == 'M2':
                dataset_id = 'N'

            sample_id = dataset_id + plate_id
            flow = float(flow)
            cai = float(cai)

            flow_dict[sample_id] = flow
            sequence_dict[sample_id] = sequence
            cai_dict[sample_id] = cai

    flow_file = out_folder + 'flow_data_raw.csv'
    sequence_file = out_folder + 'sequence_data_raw.csv'
    cai_file = out_folder + 'cai_data_raw.csv'

    with open(flow_file, 'w') as flow_out:
        for sample_id, flow in flow_dict.items():
            flow_out.write(f'{sample_id},{flow}\n')

    with open(sequence_file, 'w') as seq_out:
        for sample_id, sequence in sequence_dict.items():
            seq_out.write(f'{sample_id},{sequence}\n')

    with open(cai_file, 'w') as cai_out:
        for sample_id, cai in cai_dict.items():
            cai_out.write(f'{sample_id},{cai}\n')







def parse_rnalfold_output(rnalfold_output_folder, well_to_orf):
    well_to_local_structures = {}

    for file_name in os.listdir(rnalfold_output_folder):
        if file_name[-4:] == '.txt':
            well = file_name.split('_')[0]
            sequence = BCD + well_to_orf[well] + TERMINATOR
            rnalfold_data = RNALFoldOutput(well, sequence)
            file_location = rnalfold_output_folder + file_name

            with open(file_location, 'r') as rnalfold:
                for line in rnalfold:
                    line = line.strip()
                    if len(line.split()) == 3:
                        structure, mfe, start_pos = line.split()
                        mfe = float(mfe[1:-1].strip())

                    elif len(line.split()) == 4:

                        structure, junk, mfe, start_pos = line.split()
                        mfe = float(mfe[:-1].strip())

                    start_pos = int(start_pos)
                    if not 'A' in structure:
                        rnalfold_data.add_local_structure(sequence, structure, mfe, start_pos)

            well_to_local_structures[well] = rnalfold_data

    return well_to_local_structures


class RNALFoldOutput:
    def __init__(self, well, sequence):
        self.well = well
        self.sequence = sequence
        self.local_structures = []

    def add_local_structure(self, sequence, structure, mfe, start):
        length = len(structure)
        end = start + length - 1
        sequence_fragment = sequence[start - 1: end]
        local_structure = LocalStructure(sequence_fragment, structure, mfe, start, end)
        self.local_structures.append(local_structure)

    def get_local_structures_from_position(self, position):
        local_structures = []
        for local_structure in self.local_structures:
            if (local_structure.start - 1) <= position <= (local_structure.end - 1):
                local_structures.append(local_structure)
        return local_structures




class LocalStructure:
    def __init__(self, sequence, structure, mfe, start, end):
        self.sequence = sequence
        self.structure = structure
        self.mfe = mfe
        self.start = start
        self.end = end


def parse_structure_dir(structure_dir):
    well_to_structure = {}
    with open(structure_dir, 'r') as structure_file:
        for line in structure_file:
            line = line.strip()
            well, structure = line.split('\t')
            well_to_structure[well] = structure

    return well_to_structure

def parse_translation_dir(translation_dir):
    well_to_orf = {}

    with open(translation_dir, 'r') as translation_file:
        for line in translation_file:
            line = line.strip()
            well, orf = line.split(',')
            well_to_orf[well] = orf

    return well_to_orf


def parse_feature_importance(fi_dir):
    fi_file = open(fi_dir, 'r')
    fi_file.readline()
    feature_nr = 0
    feature_dict = {}
    
    for line in fi_file:
        line = line.strip()
        codon = line.split('\t')[0]
        if not codon in feature_dict:
            feature_dict[codon] = {}

        feature_dict[codon][feature_nr % 12] = line.split('\t')[1:]
        for i, feature_imp in enumerate(feature_dict[codon][feature_nr % 12]):
            feature_dict[codon][feature_nr % 12][i] = float(feature_imp)
        
        feature_nr += 1
    fi_file.close()
    return feature_dict

def parse_structure_feature_importance(fi_dir):
    fi_file = open(fi_dir, 'r')
    fi_file.readline()
    feature_dict = {}

    for line in fi_file:
        line = line.strip()
        base_nr = int(line.split('\t')[0])
        importances = list(map(float, line.split('\t')[1:]))
        feature_dict[base_nr] = sum(importances) / len(importances)

    fi_file.close()
    return feature_dict


def read_feature_importances(fi_file, encoding, base_to_importance = None):
    if not base_to_importance:
        base_to_importance = {}

    with open(fi_file, 'r') as fi:
        fi.readline()
        for line in fi:
            line = line.strip()
            if line:
                if encoding == 'rna-bppm-totals':
                    base, importance = line.split('\t')
                    base = int(base)
                    importance = float(importance)

                    if not base in base_to_importance:
                        base_to_importance[base] = []

                    base_to_importance[base].append(importance)
                elif encoding == 'one-hot-base' or encoding == 'one-hot-third-base':
                    base, base_type, importance = line.split('\t')
                    base = int(base)
                    importance = float(importance)

                    if not base in base_to_importance:
                        base_to_importance[base] = {}

                    if not base_type in base_to_importance[base]:
                        base_to_importance[base][base_type] = []

                    base_to_importance[base][base_type].append(importance)

                elif encoding == 'one-hot-codon':
                    codon, codon_type, importance = line.split('\t')
                    codon = int(codon)
                    importance = float(importance)

                    if not codon in base_to_importance:
                        base_to_importance[codon] = {}

                    if not codon_type in base_to_importance[codon]:
                        base_to_importance[codon][codon_type] = []

                    base_to_importance[codon][codon_type].append(importance)

                elif encoding == 'rna-bppm-onehot-third':
                    encoding_type, base, base_type, importance = line.split('\t')

                    base = int(base)
                    importance = float(importance)

                    if not base in base_to_importance:
                        base_to_importance[base] = {}

                    if encoding_type == 'bpp':
                        if not 'bpp' in base_to_importance[base]:
                            base_to_importance[base]['bpp'] = []

                        base_to_importance[base]['bpp'].append(importance)

                    elif encoding_type == 'onehot3':
                        if not base_type in base_to_importance[base]:
                            base_to_importance[base][base_type] = []

                        base_to_importance[base][base_type].append(importance)

    return base_to_importance

def read_correlations(correlations_file):
    with open(correlations_file, 'r') as correlations:
        pearson = correlations.readline().split(',')[1]
        spearman = correlations.readline().split(',')[1]

    return pearson, spearman

def read_correlations_from_dir(window_dir):
    window_to_correlations = {}
    for directory in os.listdir(window_dir):
        if directory[:7] == 'windows':
            full_path = os.path.join(window_dir, directory)
            coefficients_file = os.path.join(full_path, "correlation_coefficients.txt")
            pearson, spearman = read_correlations(coefficients_file)
            window = int(directory.split('_')[-1])
            window_to_correlations[window] = {}
            window_to_correlations[window]['pearson'] = float(pearson)
            window_to_correlations[window]['spearman'] = float(spearman)

    return window_to_correlations





def read_feature_importances_from_dir(fi_dir, encoding):
    for fi_name in os.listdir(fi_dir):
        if fi_name[-4:] == '.txt':
            fi_file = os.path.join(fi_dir, fi_name)
            base_to_importance = read_feature_importances(fi_file, encoding)

    if encoding == 'one-hot-codon' or encoding == 'one-hot-base' or encoding == 'one-hot-third-base' or encoding == 'rna-bppm-onehot-third':
        for base, type_to_importance in base_to_importance.items():
            for base_type, importances in type_to_importance.items():
                base_to_importance[base][base_type] = sum(importances) / float(len(importances))

    elif encoding == 'rna-bppm-totals':
        for base, importances in base_to_importance.items():
            base_to_importance[base] = sum(importances) / float(len(importances))

    return base_to_importance


def parse_feature_importance_per_base(fi_dir):
    fi_file = open(fi_dir, 'r')
    fi_file.readline()
    base_nr = 0
    feature_dict = {}

    for line in fi_file:

        if not base_nr in feature_dict:
            feature_dict[base_nr] = {}

        if base_nr % 4 == 0:
            feature_dict[base_nr]['A'] = list(map(float, line.split('\t')[1:]))

        elif base_nr % 4 == 1:
            feature_dict[base_nr]['C'] = list(map(float, line.split('\t')[1:]))

        elif base_nr % 4 == 2:
            feature_dict[base_nr]['G'] = list(map(float, line.split('\t')[1:]))

        elif base_nr % 4 == 3:
            feature_dict[base_nr]['T'] = list(map(float, line.split('\t')[1:]))


        base_nr += 1
    fi_file.close()
    return feature_dict

def parse_fasta(fasta_dir):
    with open(fasta_dir, 'r') as fasta_file:
        fasta_dict = {}
        sequence = []
        for line in fasta_file:
            line = line.strip()

            if line.startswith(">"):
                if sequence:
                    fasta_dict[ID] = ''.join(sequence)

                    ID = line[1:]
                    sequence = []
                else:
                    ID = line[1:]

            else:
                sequence.append(line)

        fasta_dict[ID] = ''.join(sequence)
        fasta_file.close()
    return fasta_dict



def parse_feature_importance_6(fi_dir):
    fi_file = open(fi_dir, 'r')
    fi_file.readline()
    feature_nr = 0
    feature_dict = {}

    counter = -1
    
    for line in fi_file:
        counter += 1
        line = line.strip()
        codon = line.split('\t')[0]
        if not codon in feature_dict:
            feature_dict[codon] = {}
            counter = -1

        

        feature_dict[codon][counter] = line.split('\t')[1:]
        for i, feature_imp in enumerate(feature_dict[codon][counter]):
            feature_dict[codon][counter][i] = float(feature_imp)
        
        feature_nr += 1
    fi_file.close()
    return feature_dict



class AIResult():
    def __init__(self, line):
        self.parse_line(line)

    def parse_line(self, line):
        line_info = line.split(',')
        self.well, self.orf, self.window_size = line_info[0], line_info[1],\
                                                int(line_info[2])
        

        
        self.scores = {}
        for i, score in enumerate(line_info[3:]):
            score = float(score)
            self.scores[i + (0.5 * self.window_size) + 0.5] = score
                                  
        
class AIResults():
    def __init__(self, ai_dir):
        self.ai_dir = ai_dir
        self.parse_ai()

    def parse_ai(self):
        ai_file = open(self.ai_dir)
        ai_file.readline()
        self.results = {}
        for line in ai_file:
            line = line.strip()
            if line:
                result = AIResult(line)
                self.results[result.well] = result

        ai_file.close()


def read_pairing_probabilities(bpp_file):
    well_to_bpps = {}

    with open(bpp_file, 'r') as bpp_data:
        for line in bpp_data:
            line = line.strip()
            if line:
                well_and_bpps = line.split(',')
                well = well_and_bpps[0]
                bpps = [float(x) for x in well_and_bpps[1:]]

                well_to_bpps[well] = bpps

    return well_to_bpps


def read_sequences(sequence_file):
    well_to_sequence = {}

    with open(sequence_file, 'r') as sequence_data:
        for line in sequence_data:
            line = line.strip()
            if line:
                well, sequence = line.split(',')
                well_to_sequence[well] = sequence

    return well_to_sequence


def read_flow(flow_file):
    well_to_flow = {}

    with open(flow_file, 'r') as flow_data:
        for line in flow_data:
            line = line.strip()
            if line:
                well, flow = line.split(',')
                well_to_flow[well] = float(flow)

    return well_to_flow


def read_utr_file(utr_file):
    with open(utr_file, 'r') as utr_data:
        five_utr = utr_data.readline().strip()
        terminator = utr_data.readline().strip()

    return five_utr, terminator

class TranslationResults:
    def __init__(self, translation_dir):
        self.dir = translation_dir
        self.parse_translation_results()

    def parse_translation_results(self):
        self.results = {}
        translation_file = open(self.dir, 'r')
        for line in translation_file:
            line = line.strip()
            translation_result = TranslationResult(line)
            if translation_result.well in self.results:
                print(translation_result.well)
            self.results[translation_result.well] = translation_result

        translation_file.close()

    def get_sequences(self):
        sequences = []
        for result in self.results:
            sequences.append(self.results[result].orf)

        return sequences

    def get_medium_sequences(self):
        sequences = []
        for result in self.results:
            if result.startswith('M'):
                sequences.append(self.results[result].orf)

        return sequences

    def get_medium_2_sequences(self):
        sequences = []


    def get_high_sequences(self):
        sequences = []
        for result in self.results:
            if result.startswith('H'):
                sequences.append(self.results[result].orf)

        return sequences

    def get_low_sequences(self):
        sequences = []
        for result in self.results:
            if result.startswith('L'):
                sequences.append(self.results[result].orf)

        return sequences



class TranslationResult:
    def __init__(self, line):
        self.parse_translation_result(line)

    def parse_translation_result(self, line):
        
        self.well, self.orf = line.split(',')
