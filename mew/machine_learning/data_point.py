from mew.featurisation.viennarna import from_bppm_base_totals
from mew.featurisation.one_hot import do_onehot_encoding, get_coding_range


class DataPoint:

    def __init__(self, well, sequence, flow, encoding, bpps=None, five_utr='', terminator='', length=None, utr_length=None, start_position=None, coding_length=678):
        self.flow = flow
        self.well = well
        self.sequence = sequence
        self.five_utr = five_utr
        self.terminator = terminator
        self.encoding = encoding
        self.bpps = bpps
        self.start_position = start_position
        self.coding_length = coding_length

        self.length = length
        self.utr_length = utr_length

        self.feature_vector = []
        self.feature_mapping = {}

        self.encode()

    def __repr__(self):
        return self.well

    def __hash__(self):
        return hash(self.well)

    def encode(self):

        if self.start_position != None:

            sequence_for_rna = get_coding_range(self.sequence, self.start_position, coding_length=self.coding_length).sequence

            if self.start_position < 0:
                self.five_utr = self.sequence[:abs(self.start_position)]
            else:
                self.five_utr = ''

            terminator_start = len(self.sequence) + self.start_position
            if terminator_start > self.coding_length:
                self.terminator = self.sequence[self.coding_length - self.start_position:]
            else:
                self.terminator = ''

        else:
            sequence_for_rna = self.sequence

        if self.encoding == 'rna-bppm-totals':
            if self.bpps:
                self.feature_vector = self.bpps

            else:
                self.feature_vector = from_bppm_base_totals(sequence_for_rna, self.five_utr, self.terminator)

            if self.length != None:
                self.feature_vector = self.feature_vector[:len(self.five_utr) + self.length]

            if self.utr_length != None:
                self.feature_vector = self.feature_vector[len(self.five_utr) - self.utr_length:]

        elif self.encoding == 'one-hot-base':
            if self.length != None:
                self.feature_vector = do_onehot_encoding(self.sequence[:self.length], type=self.encoding, start_position=self.start_position, coding_length=self.coding_length)
            else:
                self.feature_vector = do_onehot_encoding(self.sequence, type=self.encoding, start_position=self.start_position)

        elif self.encoding == 'one-hot-third-base':
            if self.length != None:
                self.feature_vector = do_onehot_encoding(self.sequence[:self.length], type=self.encoding, start_position=self.start_position, coding_length=self.coding_length)
            else:
                self.feature_vector = do_onehot_encoding(self.sequence, type=self.encoding, start_position=self.start_position, coding_length=self.coding_length)

        elif self.encoding == 'one-hot-codon':
            if self.length != None:
                self.feature_vector = do_onehot_encoding(self.sequence[:self.length], type=self.encoding, start_position=self.start_position, coding_length=self.coding_length)
            else:
                self.feature_vector = do_onehot_encoding(self.sequence, type=self.encoding, start_position=self.start_position, coding_length=self.coding_length)

        elif self.encoding == 'rna-bppm-onehot-third':
            if self.bpps:
                feature_vector_bpp = self.bpps

            else:
                feature_vector_bpp = from_bppm_base_totals(sequence_for_rna, self.five_utr, self.terminator)

            if self.length != None:
                feature_vector_bpp = feature_vector_bpp[:len(self.five_utr) + self.length]
                feature_vector_onehot_third = do_onehot_encoding(self.sequence[:self.length], type='one-hot-third-base', start_position=self.start_position, coding_length=self.coding_length)

            else:
                feature_vector_onehot_third = do_onehot_encoding(self.sequence, type='one-hot-third-base', start_position=self.start_position, coding_length=self.coding_length)

            if self.utr_length != None:
                feature_vector_bpp = feature_vector_bpp[len(self.five_utr) - self.utr_length:]

            self.feature_vector = feature_vector_bpp + feature_vector_onehot_third


    def set_crossval_predicted_flow(self, flow):
        self.predicted_flow = flow
        self.difference = abs(self.predicted_flow - self.flow)