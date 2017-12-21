"""
Classes for sampled sets produced by OpenFOAM
"""

import numpy as np
import os

from HbO2.parse.readSampleFiles import loadSampleFile, createSampleDirPath
from parse.case import last_time_dir


class SampledSet(object):
    """Class to read sampled sets from an OpenFOAM case."""

    def __init__(self, case_path, dir_name):
        self.case_path = case_path
        self.dir_name = dir_name

    def values(self, timeDir, set_name, field_name):
        time = float(timeDir)
        return loadSampleFile(os.path.join(self.case_path, 'domain'),
                              self.dir_name, time, set_name, [field_name])

    def last_time_values(self, set_name, field_name):
        try:
            sample_dir_name = createSampleDirPath(os.path.join(self.case_path, 'domain'),
                                                  self.dir_name)
            dir_name = last_time_dir(sample_dir_name)
        except IndexError:
            print 'Error while reading case {} (set name = {}, field name = {})' \
                .format(self.case_path, set_name, field_name)
            raise
        return self.values(dir_name, set_name, field_name)


def average_sampled_sets(path_to_sets, time, set_names, field_name, weights):
    """
    Compute a weighted average of field values on sampled sets computed by OpenFOAM.
    """
    x = []
    average = []
    for (setName, weight) in zip(set_names, weights):
        sampled_data = loadSampleFile('.', path_to_sets, time, setName, [field_name])
        x      = np.asarray(sampled_data[:,0])
        values = np.asarray(sampled_data[:,1])
        if average == []:
            average = np.zeros(values.shape)
        average += weight*values
    average /= sum(weights)
    return x, average
