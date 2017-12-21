"""
Utility functions for plotting
"""

import numpy as np
from matplotlib.axis import Axis
from matplotlib.ticker import MultipleLocator, Base


def add_multi_panel_parser_options(parser):
    """
    Add parser options for multi panel figures

    Args:
        parser (argparse.ArgumentParser): argument parser
    """
    group = parser.add_argument_group('multipanel', 'Options for figures with multiple panels')
    group.add_argument('--multipanel', action='store_true',
                       help='Whether to create multiple panels in one figure')
    group.add_argument('--panelLayout', nargs=2, type=int,
                       help='Number of panels in y and x-direction')


def add_panel_layout_parser_options(parser):
    parser.add_argument('--panelLayout', nargs=2, type=int,
                        help='Number of panels in y and x-direction')


def scale_rgb_color(color):
    """
    Scale a color given by a tuple with values from 0 to 255 to a tuple with values
    between 0 and 1.

    Args:
        color (tuple): tuple with three values between 0 and 255

    Returns:
        3-tuple with values between 0 and 1
    """
    return (color[0]/255., color[1]/255., color[2]/255.)


def annotate_axis_corner(ax, text, xy=(0.92, 0.91)):
    """
    Annotate an axis corner with a string.

    Args:
        ax (Axis): axis to annotate
        text (str): annotation text
        xy (tuple): location of the annotation, in axis fraction.
    """
    ax.annotate(r'$\mathbf{{{:s}}}$'.format(text), xy=xy, xycoords='axes fraction')


def rounded_bounds(x, roundingBase):
    """
    Rounds the function x to base, where base may be any positive float number.

    Args:
        x (array like): array for which to find bounds
        roundingBase: base for rounding
    Returns:
        2-tuple with rounded minimum and maximum.
    """
    minValue = round_to_base(np.min(x), roundingBase, func=np.floor)
    maxValue = round_to_base(np.max(x), roundingBase, func=np.ceil)
    return minValue, maxValue


def set_rounded_axis_limits(ax, base, axis='x'):
    """
    Set the axis limits by rounding the data interval to the nearest lower/upper multiple of base,

    Args:
        ax (Axis):
        base (float): rounding base
        axis (str): 'x' or 'y'
    """
    coord_axis = getattr(ax, '{:s}axis'.format(axis))
    interval = coord_axis.get_data_interval()
    rounded_interval = rounded_bounds(interval, base)
    getattr(ax, 'set_{:s}lim'.format(axis))(rounded_interval)


def include_zero_in_axis_range(ax, axis='y'):
    """
    Change the axis limits to include zero in the bounds

    Args:
        ax (Axis): axis object to adapt
        axis (str): 'x' or 'y'
    """
    lim = getattr(ax, 'get_{:s}lim'.format(axis))()
    if lim[0] < 0 and lim[1] < 0:
        lim = (lim[0], 0)
    elif lim[0] > 0 and lim[1] > 0:
        lim = (0, lim[1])
    getattr(ax, 'set_{:s}lim'.format(axis))(lim)


def contourLevels(Z, step, roundingBase):
    """Computes the contour levels from min(Z) to max(Z).

    Args:
        Z: array with variable values
        step: step between the contours
        roundingBase: base for the rounding of the min/max value of Z

    Returns:
        Array with the contour levels
    """
    minValue = round_to_base(np.min(Z), roundingBase, func=np.floor)
    maxValue = round_to_base(np.max(Z), roundingBase, func=np.ceil)
    return np.arange(minValue, maxValue+1e-6, step)


def round_to_base(x, base, func=round):
    """Rounds the function x to base, where base may be any positive float number.

    Args:
        x: number to round
        base: base for rounding
        func: rounding function (round, np.floor, np.ceil...)
    Returns:
        x rounded to a multiple of base
    """
    return base*func(float(x)/base)


def scientific_format(x, pos):
    """Returns a LaTeX string with argument x in scientific format.

    Based on http://stackoverflow.com/questions/25983218/scientific-notation-colorbar-in-matplotlib
    """
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


class MultipleWithMaxNLocator(MultipleLocator):
    """
    Locator that uses a multiple and also enforces a maximum number of ticks.
    If the maximum number of ticks is exceeded, the multiple is increment by itself.
    The determination of the number of ticks is approximative.
    """

    def __init__(self, base=1.0, maxn=10):
        super(MultipleWithMaxNLocator, self).__init__(base=base)
        self._maxn = maxn

    def tick_values(self, vmin, vmax):
        if vmax < vmin:
            vmin, vmax = vmax, vmin
        n_tick = 1e10
        original_base = self._base.get_base()
        base = 0
        locs = np.zeros(0,)
        while n_tick > self._maxn:
            base += original_base
            base_obj = Base(base)
            vmin = base_obj.ge(vmin)
            n = (vmax - vmin + 0.001 * base) // base
            locs = vmin - base + np.arange(n + 3) * base
            n_tick = len(locs) - 2
        return self.raise_if_exceeds(locs)

