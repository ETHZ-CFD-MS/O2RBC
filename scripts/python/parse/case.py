import os

def time_dirs(path):
    """Return the sorted list of time folders in the given path.

    Args:
        path: relative path with respect to the case name

    Returns:
        Sorted list of strings with time directories
    """
    try:
        dirs = next(os.walk(path))[1]
    except StopIteration:
        print 'Empty directory {}'.format(path)
        raise
    return sorted(filter(is_float, dirs))


def last_time_dir(path):
    """Return the last time folder in the given path

    Args:
        path: relative path with respect to the case name

    Returns:
        string, name of last time directory
    """
    dirs = time_dirs(path)
    try:
        return dirs[-1]
    except IndexError:
        print 'Error while reading time directories in {:s}'\
            .format(path)
        raise


def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False