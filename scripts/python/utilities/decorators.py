"""
Module with generale use decorators
"""

from functools import wraps

import numpy as np


def lazy_property(func):
    """
    Decorator that makes a property lazily evaluated.
    """
    attr_name = '_lazy_' + func.__name__

    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, func(self))
        return getattr(self, attr_name)
    return _lazy_property


def lazy_function(func):
    """
    Decorator that makes a function with arguments lazily evaluated.

    Currently, only functions with hashable arguments are supported.
    """
    dict_name = '_lazy_dict_' + func.__name__

    @wraps(func)
    def _lazy_function(self, *args):
        if not hasattr(self, dict_name):
            setattr(self, dict_name, {})
        if args not in getattr(self, dict_name):
            getattr(self, dict_name)[args] = func(self, *args)
        return getattr(self, dict_name)[args]
    return _lazy_function


def remove_lazy_properties(obj):
    """
    Removes the lazy properties stored in an object.

    Args:
        func (obj):
    """
    for attr_name in dir(obj):
        if attr_name.startswith('_lazy_'):
            delattr(obj, attr_name)


def accept_scalar_arg(f):
    """
    Decorator that returns a function that accepts scalar and array
    arguments, from a function that initially only accepts array arguments.
    """
    @wraps(f)
    def f_scalar_or_array_arg(self, x):
        x, scalar_input = to_array(x)
        y = f(self, x)
        if scalar_input:
            np.squeeze(y)
            y = y[0]
        return y
    return f_scalar_or_array_arg


def to_array(x):
    """
    Return a numpy array constructed from x and a flag is_scalar.

    If the number of dimensions of the resulting array is 0, the flag is_scalar
    is set to True.
    """
    x = np.asarray(x)
    is_scalar = False
    if x.ndim == 0:
        is_scalar = True
        x = x[None]
    return x, is_scalar

