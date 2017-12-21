"""
Functions to read data files and return them in given data structures.
"""

import csv
import re

import numpy as np

# Regular expression to extract a key and a value from a line such as
# "key  value; // whatever"
key_space_value_semicolon_pattern = """
\s*     # allowed initial white spaces
(\w+)   # matches keyword (only alphanumeric)
\s+     # at least one white space
([^;]+) # matches value (only semi-colon not allowed)
;       # final semi-colon
"""


def readKeysFromFile(fileObject, pattern, filter_keys=None):
    """Return a dictionary with key-value pairs extracted from a given file object.

    The value is converted to a float whenever the conversion is possible.

    Args:
        fileObject: file object
        pattern (str): regular expression pattern, can be verbose
        filter_keys (list): iterable that specifies selected keys to extract.
            If it is None, all keys are extracted.

    Returns:
        dictionary with extracted key-value pairs.
    """
    matches = {}
    for line in fileObject:
        try:
            (key, value) = readKeyValueFromLine(line, pattern)
            try:
                value = float(value)
            except ValueError:
                pass
            if filter_keys is None or key in filter_keys:
                matches[key] = value
        except TypeError:
            pass
    return matches


def readKeyValueFromLine(line, pattern):
    """
    Extract a key and a value from a line using a regular expression.

    Args:
        line (str): line to parse
        pattern (str): verbose regular expression pattern with two groups, the first
                       for the key and the second for the value
    Returns:
        If found, a tuple (key, value)
    """
    m = re.match(pattern, line, re.VERBOSE)
    if m:
        return m.group(1, 2)


def load_numeric_to_dictionary(fo):
    """Read a file with one header line and corresponding numeric data.

    Args:
        fo: file object to read
    Returns:
        Dictionary with keys corresponding to the headers, and values which are
        numpy arrays with the columns values.
    """
    results = {}
    line = fo.readline()
    headers = line.split()
    for h in headers:
        results[h] = []
    for line in fo:
        values = [float(x) for x in line.split()]
        for h, x in zip(headers, values):
            results[h].append(x)
    for h in headers:
        results[h] = np.asarray(results[h])
    return results


def load_csv_to_dictionary(fo, delimiter='\t'):
    """
    Read a CSV file with a header line and return the columns in a dictionary
    with corresponding keys

    Args:
        fo: file object to read
        delimiter: delimiter in csv file

    Returns:
        dictionary with keys corresponding to column headers
    """
    num_lines = sum(1 for line in fo)
    fo.seek(0, 0)
    header_lines = 1
    n_data = num_lines - header_lines
    keys = []
    ncol = 0
    data = np.array(())
    reader = csv.reader(fo, delimiter=delimiter, quotechar='"')
    i = 0
    for line in reader:
        if i == 0:
            keys = line
            if keys[-1] == '':  # remove trailing separator
                keys.pop()
            ncol = len(keys)
            data = np.zeros((n_data, ncol))
        else:
            data[i-1, :] = [float(x) for x in line[0:ncol]]
        i += 1
    results = {}
    for i, key in enumerate(keys):
        results[key] = data[:,i]
    return results


def load_csv_to_array(fo, delimiter='\t', n_header_lines=0):
    """
    Read a CSV file and return the values in a numpy array.

    Args:
        fo: file object to read
        delimiter: delimiter in csv file
        n_header_lines (int): number of header lines

    Returns:
        np.ndarray
    """
    num_lines = sum(1 for line in fo)
    fo.seek(0, 0)
    n_data = num_lines - n_header_lines
    ncol = 0
    data = np.array(())
    reader = csv.reader(fo, delimiter=delimiter, quotechar='"')
    i = 0
    for line in reader:
        if i == 0:
            keys = line
            if line[-1] == '':  # remove trailing separator
                line.pop()
            ncol = len(line)
            data = np.zeros((n_data, ncol))
        if i >= n_header_lines:
            data[i - n_header_lines, :] = [float(x) for x in line[0:ncol]]
        i += 1
    return data
