from mew.featurisation.viennarna import from_bppm_base_totals
from mew.featurisation.one_hot import *

class DataPoint:

    def __init__(self, well, sequence, flow, encoding, bcd=''):
        self.flow = flow
        self.well = well
        self.sequence = sequence
        self.bcd = bcd
        self.encoding = encoding

        self.feature_vector = []
        self.feature_mapping = {}

        self.encode()

    def __repr__(self):
        return self.well

    def __hash__(self):
        return self.well

    def encode(self):
        if self.encoding == 'rna-bppm-totals':
            self.feature_vector, self.feature_mapping = from_bppm_base_totals(self.sequence, bcd)
        elif self.encoding == 'one-hot':
            self.feature_vector, self.feature_mapping = None, None


    def set_crossval_predicted_flow(self, flow):
        self.predicted_flow = flow
        self.difference = abs(self.predicted_flow - self.flow)