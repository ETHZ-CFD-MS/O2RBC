import numpy as np
import os

from parse.readfile import load_csv_to_array, load_numeric_to_dictionary


class HemoglobinOnSegmentsResultReader(object):
    """
    Read postprocessing results for hemoglobin saturation on segments and provide
    access to data.
    """
    segment_key = 'segment'

    def __init__(self, case_path):
        self.case_path = case_path
        from HbO2.postprocess.hemoglobingraph import HemoglobinOnSegmentsPostProcessor
        self.result_filename = HemoglobinOnSegmentsPostProcessor.output_file_name
        self.hemoglobin_graph_results = {}
        self._load_hemoglobin_graph_results()

    def edge_ids(self):
        return self.hemoglobin_graph_results[self.segment_key]

    def result_dict(self):
        return self.hemoglobin_graph_results

    def _load_hemoglobin_graph_results(self):
        data = load_numeric_to_dictionary(open(os.path.join(self.case_path, self.result_filename)))
        sorted_ids = np.argsort(data[self.segment_key])
        data[self.segment_key] = data[self.segment_key].astype(int)
        for key in data:
            data[key] = data[key][sorted_ids]
        self.hemoglobin_graph_results = data


class TissueVolumesReader(object):
    """
    Read topological tissue volumes stored in file.
    """

    def __init__(self, case_path, topological_volume_file):
        self.case_path = case_path
        self.topological_volumes_file = topological_volume_file
        self.eids_topological = []
        self.topological_tissue_volumes_list = []
        self._load_tissue_volumes()

    def topological_tissue_volumes(self):
        return np.array(self.topological_tissue_volumes_list)

    def edge_ids(self):
        return self.eids_topological

    def _load_tissue_volumes(self):
        data = load_csv_to_array(open(os.path.join(self.case_path, self.topological_volumes_file)),
                                 delimiter=' ', n_header_lines=0)
        sorted_ids = np.argsort(data[:, 0])
        self.eids_topological = list(data[sorted_ids, 0].astype(int))
        self.topological_tissue_volumes_list = data[sorted_ids, 1]

