"""
Utility functions for container manipulation.
"""
from itertools import izip

import numpy as np


def merge_list_from_dict_list(key, dicts):
    """
    Return a list which is the concatenated version of lists that correspond to a given key
    in a dictionary list.

    Args:
        key (str): key
        dicts (list): list of dictionaries that have a array-like value for the given key.

    Returns:
        A list with concatenated values of dicts[0][key], dicts[1][key]...
    """
    merge_list = dicts[0][key]
    for dictionary in dicts[1:]:
        merge_list = np.append(merge_list, dictionary[key], axis=0)
    return merge_list


def concatenate_array_list(array_list, axis=0):
    """
    Concatenate a list of numpy array
    Args:
        array_list (list): list of numpy arrays
        axis (int): axis along which values are appended

    Returns:
        np.ndarray, concatenated array
    """
    merged_array = array_list[0]
    if len(array_list) > 1:
        for vals in array_list[1:]:
            merged_array = np.append(merged_array, vals, axis=axis)
    return merged_array


def unique_everseen(seq):
    """
    Return a list of unique elements in seq while preserving the original order of the
    first occurrences.

    Args:
        seq (sequence):

    Returns:
        list of unique ordered elements
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def strictly_increasing(seq):
    """
    Check if the provided sequence has strictly increasing values.

    Args:
        seq (sequence):

    Returns:
        bool
    """
    return all(x < y for x, y in zip(seq, seq[1:]))


def check_equal(iterator):
    """
    Check whether the elements of the iterator are equal.

    Two additional nesting levels of iterators are supported. However this is not well tested
    yet...

    Args:
        iterator (seq):

    Returns:
        bool
    """
    try:
        iterator = iter(iterator)
        first = next(iterator)
        try:
            return all(first == rest for rest in iterator)
        except ValueError:
            try:
                return all(all(x1 == x2 for x1, x2 in izip(first, rest)) for rest in iterator)
            except ValueError:
                return all(all(all(x1 == x2) for x1, x2 in izip(first, rest)) for rest in iterator)
    except StopIteration:
        return True
