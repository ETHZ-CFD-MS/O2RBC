"""
Classes for OpenFOAM probes
"""

import numpy as np
import os
import re

from parse.case import time_dirs
from utilities.containers import check_equal, concatenate_array_list, strictly_increasing

float_min = -1e300
float_max = 1e300


class ProbeValue(object):

    def __init__(self, case_name, dir_name, field_name):
        self.case = case_name
        if os.path.isdir(os.path.join(case_name, 'domain', 'postProcessing')):
            self.dir_path = os.path.join(case_name, 'domain', 'postProcessing', dir_name)
        else:
            self.dir_path = os.path.join(case_name, 'domain', dir_name)
        self.time_dirs = time_dirs(self.dir_path)
        self.field_name = field_name
        self._times = None
        self._values = None
        self._probe_positions = None
        self.read()

    def read(self):
        times = [[]]*len(self.time_dirs)
        values = [[]]*len(self.time_dirs)
        probe_positions = [[]]*len(self.time_dirs)
        for i, time_dir in enumerate(self.time_dirs):
            times[i], values[i], probe_positions[i] = self.load_one_file(time_dir)
        if not check_equal(probe_positions):
            raise RuntimeError('The probe positions in {:s} are not equal'.format(self.dir_path))
        self._probe_positions = probe_positions[0]
        merged_times = concatenate_array_list(times)
        self._check_times(merged_times)
        self._times = merged_times
        self._values = concatenate_array_list(values)

    def times(self, min_time=float_min, max_time=float_max):
        """
        Return the probe times.
            min_time (float): minimal time to return
            max_time (float): maximal time to return

        Returns:
            np.ndarray with time values
        """
        return self._times[self.indices_in_range(min_time, max_time)]

    def values(self, position=None, min_time=float_min, max_time=float_max):
        """
        Return the field values at each time

        Args:
            position (int): index of the probe, starting with zero.
                If none, uses all positions
            min_time (float): minimal time to return
            max_time (float): maximal time to return

        Returns:
            np.ndarray with one column if the argument 'position' is used, or
            as many columns as probe positions.

        """
        indices = self.indices_in_range(min_time, max_time)
        if position is not None:
            return np.squeeze(self._values[indices, position])
        else:
            return self._values[indices, :]

    def positions(self, position=None):
        """
        Return the probe positions

        Args:
            position (int): index of the probe, starting with zero.
                If none, returns all positions.

        Returns:
        """
        if position is not None:
            return np.squeeze(self._probe_positions[position,:])
        else:
            return self._probe_positions

    def position_indices(self):
        """
        Return an index list for position indexing.
        """
        return range(self._probe_positions.shape[0])

    def time_range(self):
        """
        Return the range of times.

        Returns:
            tuple with minimum and maximum
        """
        return min(self.times()), max(self.times())

    def load_one_file(self, time_dir):
        """
        Load the content of one file with probe data with the following format:

        # Probe 0 (x_0 y_0 z_0)
        # Probe 1 (x_1 y_1 z_1)
        # ...
        # Probe n (x_n y_n z_n)
        #         Probe               0               1             ...               n
        #          Time
                 0.0005       79.172597       78.085073             ...       77.758809
                  0.001       78.826741       78.103028             ...         77.7598
                    ...

        Args:
            time_dir (str): time directory name

        Returns:
            times (np.ndarray): time values (n_times x 1)
            values (np.ndarray): field values at probes (n_times x n_probes)
            probe_positions (np.ndarray): probe coordinates (n_probes x 3)

        """
        path_to_probe = os.path.join(self.dir_path, time_dir, self.field_name)
        try:
            n_times = sum(not line.startswith('#') for line in open(path_to_probe))
            n_probes = sum(line.startswith('#') for line in open(path_to_probe)) - 2
        except IOError:
            print 'Could not open file {:s}'.format(path_to_probe)
            raise
        times = np.zeros((n_times,))
        probe_positions = np.zeros((n_probes, 3))
        with open(path_to_probe, 'r') as f:
            values = np.zeros((n_times, n_probes))
            for i, line in enumerate(f):
                if i < n_probes:
                    # extract coordinates enclosed in parentheses
                    match = re.search('\((.*)\)', line)
                    coord_string = match.group(1)
                    coords = [float(r) for r in coord_string.split()]
                    probe_positions[i, :] = coords
                elif i >= n_probes + 2:
                    t = i - n_probes - 2
                    split_line = [float(x) for x in line.split()]
                    times[t] = split_line[0]
                    values[t,:] = split_line[1:]
        return times, values, probe_positions

    def indices_in_range(self, min_time=float_min, max_time=float_max):

        indices = (self._times >= min_time) & (self._times <= max_time)
        if not np.any(indices):
            raise RuntimeError('No indices with time values in the range {:g}, {:g}'
                               .format(min_time, max_time))
        return np.where(indices)

    def _check_times(self, times):
        """
        Check that the times array
        Args:
            times (np.ndarray): time values

        Raises:
            ValueError if times are not strictly increasing
        """
        if not strictly_increasing(times):
            ValueError('The time values in {:s} are not strictly increasing')

