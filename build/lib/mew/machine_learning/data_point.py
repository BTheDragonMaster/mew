from mew.featurisation.viennarna import from_bppm_base_totals
from mew.featurisation.one_hot import do_onehot_encoding


class DataPoint:

    def __init__(self, well, sequence, flow, encoding, bpps=None, five_utr='', terminator='', length=None, utr_length=None):
        self.flow = flow
        self.well = well
        self.sequence = sequence
        self.five_utr = five_utr
        self.terminator = terminator
        self.encoding = encoding
        self.bpps = bpps

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
        if self.encoding == 'rna-bppm-totals':
            if self.bpps:
                self.feature_vector = self.bpps

            else:
                self.feature_vector = from_bppm_base_totals(self.sequence, self.five_utr, self.terminator)

            if self.length != None:
                self.feature_vector = self.feature_vector[:len(self.five_utr) + self.length]

            if self.utr_length != None:
                self.feature_vector = self.feature_vector[len(self.five_utr) - self.utr_length:]

        elif self.encoding == 'one-hot-base':
            if self.length != None:
                self.feature_vector = do_onehot_encoding(self.sequence[:self.length], type=self.encoding)
            else:
                self.feature_vector = do_onehot_encoding(self.sequence, type=self.encoding)

        elif self.encoding == 'one-hot-third-base':
            if self.length != None:
                self.feature_vector = do_onehot_encoding(self.sequence[:self.length], type=self.encoding)
            else:
                self.feature_vector = do_onehot_encoding(self.sequence, type=self.encoding)

        elif self.encoding == 'one-hot-codon':
            if self.length != None:
                self.feature_vector = do_onehot_encoding(self.sequence[:self.length], type=self.encoding)
            else:
                self.feature_vector = do_onehot_encoding(self.sequence, type=self.encoding)


    def set_crossval_predicted_flow(self, flow):
        self.predicted_flow = flow
        self.difference = abs(self.predicted_flow - self.flow)