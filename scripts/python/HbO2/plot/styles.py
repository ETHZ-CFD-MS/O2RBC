
from matplotlib.font_manager import FontProperties
import warnings

from plot.utils import scale_rgb_color


class StyleSchemeCOSH(dict):

    def __init__(self):
        super(StyleSchemeCOSH, self).__init__()

        # line widths
        mean_width = 0.35
        mean_pm_std_width = 1
        std_width = 1
        cap_width = 1
        diff_width = 1
        individual_rbc_width = 0.3
        individual_rbc_whole_path_width = 1.0

        # colors

        # brown/orange/dark red palette
        rbc_inter_color = scale_rgb_color((150, 72, 0))
        nonlinear_color = scale_rgb_color((205, 106, 30))
        explicit_color = scale_rgb_color((156, 34, 90))
        linear_color = scale_rgb_color((255, 156, 0))
        no_model_color = scale_rgb_color((238, 189, 168))
        exp_fit_color = scale_rgb_color((200, 100, 100))

        # green/blue palette
        # nonlinear_color = scale_rgb_color((50, 216, 212))
        # explicit_color = scale_rgb_color((0, 100, 0))
        # linear_color = scale_rgb_color((50, 255, 72))
        # no_model_color = scale_rgb_color((176, 208, 159))

        individual_rbc_color = scale_rgb_color((178, 222, 166))
        # individual_rbc_whole_path_color = scale_rgb_color((55, 220, 0))
        individual_rbc_whole_path_color = scale_rgb_color((100, 100, 100))
        # individual_rbc_symbol_color = scale_rgb_color((20, 100, 0))
        individual_rbc_symbol_color = scale_rgb_color((0, 0, 0))
        twin_axes_color = scale_rgb_color((0, 0, 176))

        simul_radii_face_color = scale_rgb_color((0, 0, 220))
        simul_radii_edge_color = 'k'
        int_func_radii_face_color = scale_rgb_color((160, 0, 0))
        # int_func_radii_random_inlet_face_color = scale_rgb_color((240, 80, 80))
        int_func_radii_random_inlet_face_color = simul_radii_face_color
        int_func_radii_edge_color = 'k'
        int_topol_radii_face_color = 'k'
        int_topol_radii_edge_color = 'k'
        int_topol_radii_marker_face_color = (1, 1, 1, 0)
        int_topol_radii_marker_edge_color = (0.2, 0.2, 0.2)

        mean_dashes = ()
        rbc_inter_dashes = (2.4, 1.2, 0.6, 1.2)
        nonlinear_dashes = rbc_inter_dashes
        explicit_dashes = (3.6, 3.6)
        linear_dashes = (0.75, 1.125)
        no_model_dashes = (0.6, 2.4)
        exp_fit_dashes = ()

        labels = {'sim': 'simulation',
                  'RBCInter': 'RBC inter.',
                  'linearized': 'linearized',
                  'nonlinear': 'nonlinear',
                  'explicit': 'explicit',
                  'linearizedODE': 'linearized',
                  'equal_outfluxes': 'equal flux',
                  'expfit': 'exp. fit'}

        self['meanSim'] = {'dashes': mean_dashes, 'color': 'k', 'linewidth': mean_width,
                           'label': labels['sim']}
        self['meanPmStdSim'] = {'linestyle': '-', 'color': 'k', 'linewidth': mean_pm_std_width,
                                'label': labels['sim']}
        self['meanRBCInter'] = {'dashes': mean_dashes, 'color': rbc_inter_color,
                                'linewidth': mean_width, 'label': labels['RBCInter']}
        self['meanPmStdRBCInter'] = {'dashes': rbc_inter_dashes, 'color': rbc_inter_color,
                                     'linewidth': mean_pm_std_width, 'label': labels['RBCInter']}
        self['meanPmStdNoModel'] = {'dashes': no_model_dashes, 'color': no_model_color,
                                    'linewidth': mean_pm_std_width, 'label': labels['equal_outfluxes']}
        self['stdSim'] = {'linestyle': '-', 'color': 'k', 'linewidth': std_width,
                          'label': labels['sim']}
        self['stdRBCInter'] = {'dashes': rbc_inter_dashes, 'color': rbc_inter_color,
                               'linewidth': std_width, 'label': labels['RBCInter']}
        self['stdLinear'] = {'dashes': linear_dashes, 'color': linear_color,
                             'linewidth': std_width, 'label': labels['linearized']}
        self['stdExpFit'] = {'dashes': exp_fit_dashes, 'color': exp_fit_color,
                             'linewidth': std_width, 'label': labels['expfit']}
        self['stdNoModel'] = {'dashes': no_model_dashes, 'color': no_model_color,
                              'linewidth': std_width, 'label': labels['equal_outfluxes']}
        self['individualRBC'] = {'linestyle': '-', 'color': individual_rbc_color,
                                 'linewidth': individual_rbc_width, 'alpha': 1.0}

        self['meanNonlinear'] = {'dashes': nonlinear_dashes, 'color': nonlinear_color,
                                 'linewidth': mean_width, 'label': labels['nonlinear']}
        self['meanExplicit'] = {'linestyle': '-', 'dashes': explicit_dashes, 'color': explicit_color,
                                'linewidth': mean_width, 'label': labels['explicit']}
        self['capSim'] = {'linestyle': '-', 'color': 'k', 'linewidth': cap_width, 'label': labels['sim']}
        self['capNonlinear'] = {'dashes': nonlinear_dashes, 'color': nonlinear_color,
                                'linewidth': cap_width, 'label': labels['nonlinear']}
        self['capExplicit'] = {'dashes': explicit_dashes, 'color': explicit_color,
                               'linewidth': cap_width, 'label': labels['explicit']}
        self['capNoModel'] = {'dashes': no_model_dashes, 'color': no_model_color,
                              'linewidth': cap_width, 'label': labels['equal_outfluxes']}
        self['diffSim'] = {'linestyle': '-', 'color': 'k', 'linewidth': diff_width,
                           'label': labels['sim']}
        self['diffNonlinear'] = {'dashes': nonlinear_dashes, 'color': nonlinear_color,
                                 'linewidth': diff_width, 'label': labels['nonlinear']}
        self['diffExplicit'] = {'dashes': explicit_dashes, 'color': explicit_color,
                                'linewidth': diff_width, 'label': labels['explicit']}
        self['diffLinear'] = {'dashes': linear_dashes, 'color': linear_color,
                              'linewidth': diff_width, 'label': labels['linearized']}
        self['diffExpFit'] = {'linestyle': ':', 'color': 'y', 'linewidth': diff_width,
                              'label': labels['expfit']}
        self['diffNoModel'] = self['capNoModel']

        self['twinAxis'] = {'linestyle': '-', 'color': twin_axes_color}
        self['MVN1'] = {'markeredgecolor': scale_rgb_color((0, 0, 160))}
        self['MVN2'] = {'markeredgecolor': scale_rgb_color((140, 0, 70))}

        self['simul_radii_face_color'] = simul_radii_face_color
        self['simul_radii_edge_color'] = simul_radii_edge_color
        self['int_func_radii_face_color'] = int_func_radii_face_color
        self['int_func_radii_random_inlet_face_color'] = int_func_radii_random_inlet_face_color
        self['int_func_radii_edge_color'] = int_func_radii_edge_color
        self['int_topol_radii_face_color'] = int_topol_radii_face_color
        self['int_topol_radii_edge_color'] = int_topol_radii_edge_color
        self['int_topol_radii_marker_face_color'] = int_topol_radii_marker_face_color
        self['int_topol_radii_marker_edge_color'] = int_topol_radii_marker_edge_color

        self['individualRBCWholePath'] = {'linestyle': '-', 'color': individual_rbc_whole_path_color,
                                          'linewidth': individual_rbc_whole_path_width, 'alpha': 0.01}
        self['distalHbMarker'] = {'marker': '.',
                                  'markeredgecolor': individual_rbc_symbol_color + (0.20,),
                                  'markerfacecolor': (1, 1, 1, 0),
                                  'markersize': 2.5,
                                  'markeredgewidth': 0.5}


def create_COSH_legend(ax, **kwargs):
    """
    Create a legend for the axis in argument.

    Args:
        ax (pyplot.Axis): axis object

    Returns:
        created legend object

    """
    fontP = FontProperties()
    fontP.set_size('xx-small')
    props = dict({'prop': fontP,
                  'fancybox': True,
                  'borderpad': 0.2,
                  'labelspacing': 0.3,
                  'handlelength': 3.6},
                 **kwargs)
    legend = ax.legend(**props)
    try:
        legend.get_frame().set_linewidth(0.5)
    except AttributeError:
        warnings.warn("No legends found")
    return legend

def set_COSH_rc_params():
    import matplotlib as mpl
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'

