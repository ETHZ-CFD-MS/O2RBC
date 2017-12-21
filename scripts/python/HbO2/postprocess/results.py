"""
Interfaces for simulation results
"""

from HbO2.parse.readhemoglobin import HemoglobinOnSegmentsResultReader


class CaseResults(object):
    """
    General interface for case results
    """

    def __init__(self, case_path):
        super(CaseResults, self).__init__()
        self.case_path = case_path


class HemoglobinOnSegmentsResults(CaseResults):
    """
    Interface for results produced by HemoglobinOnSegmentsPostProcessor
    """

    def __init__(self, case_path):
        super(HemoglobinOnSegmentsResults, self).__init__(case_path)
        self.reader = HemoglobinOnSegmentsResultReader(case_path)
        self.result_dict = self.reader.result_dict()

    def __getitem__(self, item):
        return self.result_dict[item]

    def result_names(self):
        return self.result_dict.keys()

    def hb_mean_upstream(self):
        return self['hbMeanUp']

    def hb_mean_downstream(self):
        return self['hbMeanDown']

    def hb_mean_on_edge(self):
        return 0.5*(self.hb_mean_upstream() + self.hb_mean_downstream())

    def hb_mean_drop(self):
        return self.hb_mean_upstream() - self.hb_mean_downstream()

    def rbc_flow(self):
        return self['meanRBCFlow']

    def rbc_velocity(self):
        return self['meanVelocity']

    def hematocrit(self):
        return self['meanHT']

    def mean_arrival_transit_time(self):
        return self['meanArrivalTT']

    def mean_arrival_path_length(self):
        return self['meanArrivalPath']

    def edge_ids(self):
        return self['segment']


class MultiCaseResults(object):
    """
    Interface to access the same set of results from different simulations

    Attributes:
        results (list): list of instances of CaseResults
        combine_op (func): function that take a single iterable as an argument
                           and combines its elements as desired by the user
    """

    def __init__(self, results, combine_op):
        super(MultiCaseResults, self).__init__()
        self.results = results
        self.combine_op = combine_op

    def __getattr__(self, item):
        if all([callable(getattr(result, item)) for result in self.results]):
            return self.combine_result_functions(item)
        else:
            raise AttributeError("'{:s}' object has no attribute {:s}".format(type(self), item))

    def topological_tissue_radii(self):
        return self.combine_op([result.topological_tissue_radii()
                                for result in self.results])

    def combine_result_functions(self, func_name):
        def combined_func():
            return self.combine_op([getattr(result, func_name)()
                                    for result in self.results])
        return combined_func


