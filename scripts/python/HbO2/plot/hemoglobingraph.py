"""
Plotter objects for hemoglobin saturation in a capillary network.
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.transforms as transforms
import numpy as np
import scipy.stats
import warnings

from HbO2.plot import styles
from HbO2.plot.labels import setXLabel, setYLabel
from HbO2.plot.simulated import GraphCasePlotter
from HbO2.postprocess.hemoglobingraph import HemoglobinOnSegmentsPostProcessor, HemoglobinOnWholePathsPostProcessor
from plot.utils import rounded_bounds, round_to_base, annotate_axis_corner

style_scheme = styles.StyleSchemeCOSH()


class HemoglobinOnSegmentsPlotter(GraphCasePlotter):
    """
    Plotter object for hemoglobin saturation on segments in a capillary network.

    Attributes:
        n_profile_min (int): minimal number of required profiles for some plots
    """

    edge_type_symbols = {'inlet': {'marker': 'd',
                                   'markeredgecolor': 'k',
                                   'markerfacecolor': (1, 1, 1, 0)},
                         'convBif': {'marker': '>',
                                     'markeredgecolor': 'k',
                                     'markerfacecolor': (1, 1, 1, 0)},
                         'noType': {'marker': 'x',
                                    'markeredgecolor': 'k',
                                    'markerfacecolor': (1, 1, 1, 0)},
                         'outlet': {'marker': '*',
                                    'markeredgecolor': 'k',
                                    'markerfacecolor': (1, 1, 1, 0)},
                         'divBif': {'marker': '<',
                                    'markeredgecolor': 'k',
                                    'markerfacecolor': (1, 1, 1, 0)}}

    edge_type_method = {'inlet': 'inlet_edges',
                        'convBif': 'post_converging_bifurcation_edges'}
    # 'outlet': 'outlet_edges',
    # 'divBif': 'post_diverging_bifurcation_edges'}

    def __init__(self, postprocessor, **kwargs):
        """
        Constructor

        Args:
            postprocessor (HemoglobinOnSegmentsPostProcessor): postprocessor
            **n_profile_min (int): minimal number of profile required for plotting
        """
        super(HemoglobinOnSegmentsPlotter, self).__init__(postprocessor, enforce_positive_flow=True)
        self.n_profile_min = kwargs.get('n_profile_min', 2)

    def edge_ids(self):
        eids = self.postprocessor.rbcDataPostProcessor.edgeIdsWithPaths(nPathMin=self.n_profile_min)
        return [ei for ei in eids if self.n_profiles(ei) >= self.n_profile_min]

    def n_profiles(self, si):
        return self.postprocessor.n_rbc_average(si)

    def plot_all_profiles(self):
        for ei in self.edge_ids():
            self.plot_hb_profiles(ei)
            self.plot_hb_mean_pm_std(ei)

    def plot_hb_mean_pm_std(self, si):
        """
        Plot the hemoglobin saturation mean plus/minus its standard deviation on segment si.
        """
        self.plotFieldAverageWithStdProfile('Hb_mean', si,
                                            nAverage=self.n_profiles(si),
                                            meanStyle=style_scheme['meanSim'],
                                            stdStyle=style_scheme['meanPmStdSim'])

    def plot_hb_profiles(self, si):
        """
        Plot individual Hb profiles on segment si.
        """
        self.plotFieldProfiles('Hb_mean', si, self.n_profiles(si),
                               style=style_scheme['individualRBC'])

    def plot_hb_profile_integrated(self, si):
        """
        Plot individual Hb profiles on segment si
        """
        x_edge = self.xPlot(si)
        xmin = min(x_edge)
        xmax = max(x_edge)
        graph = self.postprocessor.integrator.graph
        hb_up = graph.es(edge_index=si)['upstream_mean_hb']
        hb_down = graph.es(edge_index=si)['downstream_mean_hb']
        if hb_up and hb_down:
            plt.plot([xmin, xmax], [hb_up, hb_down])
        else:
            warnings.warn('Missing integrated hemoglobin values on edge {:d}'.format(si))

    def restrict_x_limits_to_defined_values(self, si):
        """
        Restrict the x-axis limits to values that are defined for all plotted RBC paths
        """
        s = self.sCoords(si)
        defined_s_values = self.postprocessor.rbcDataPostProcessor.\
                           sCoordIntervalWithDefinedValues('Hb_mean', s, si,
                                                           nPath=self.n_profiles(si))
        defined_x_values = self.sToXPlot(np.array(defined_s_values), si)
        plt.gca().set_xlim(min(defined_x_values), max(defined_x_values))

    def plot_std_difference_vs_upstream_std(self, **kwargs):
        self.plot(self.postprocessor.upstream_std, self.postprocessor.std_difference, **kwargs)
        setXLabel(r'\sigma_{S,a}', '')
        setYLabel(r'\sigma_{S,v} - \sigma_{S,a}', '')

    def plot_std_slope_vs_upstream_std(self, **kwargs):
        self.plot(self.postprocessor.upstream_std, self.postprocessor.std_slope,
                  y_factor=1e-3, **kwargs)
        plt.gca().yaxis.set_major_locator(MultipleLocator(0.1))
        setXLabel(r'\sigma_{S,a}', '')
        setYLabel(r'(\sigma_{S,v} - \sigma_{S,a})/l', r'\mathrm{mm}^{-1}')

    def plot_std_slope_vs_upstream_mean(self, **kwargs):
        self.plot(self.postprocessor.upstream_mean, self.postprocessor.std_slope,
                  y_factor=1e-3, **kwargs)
        setXLabel(r'\mu_{S,a}', '')
        setYLabel(r'(\sigma_{S,v} - \sigma_{S,a})/l', r'\mathrm{mm}^{-1}')

    def plot_std_slope_vs_mean_slope(self, **kwargs):
        self.plot(self.postprocessor.mean_slope, self.postprocessor.std_slope,
                  x_factor=1e-4, y_factor=1e-4, **kwargs)
        setXLabel(r'(\mu_{S,a} - \mu_{S,v})/l', r'(100 \; \mathrm{\muup m)^{-1}}')
        setYLabel(r'(\sigma_{S,v} - \sigma_{S,a})/l', r'(100 \; \mathrm{\muup m)^{-1}}')

    def plot_std_slope_vs_length(self, **kwargs):
        self.plot(self.postprocessor.scoord_interval_length, self.postprocessor.std_slope,
                  x_factor=1e6, y_factor=1e-3, **kwargs)
        setXLabel(r'\mathrm{edge\;length}', 'um')
        setYLabel(r'(\sigma_{S,v} - \sigma_{S,a})/l', r'\mathrm{mm}^{-1}')

    def plot_hb_drop_vs_time_to_previous_rbc(self, si, **kwargs):
        threshold = kwargs.get('threshold', np.inf)
        style = kwargs.get('style', {'color': 'b', 'markersize': 3})
        time_difference = self.postprocessor.rbcDataPostProcessor.\
                          timeToPreviousRBC(si, self.n_profiles(si))
        hb_drop = self.postprocessor.hb_difference(si)[1:]
        plt.plot(1e3*time_difference, hb_drop, '.', **style)
        # plot linear regression and annotate with r^2 value
        res = self.postprocessor.linregress_hb_drop_with_time_to_previous_rbc(si, threshold=threshold)
        slope = res[0]
        intercept = res[1]
        r_value = res[2]
        x = np.linspace(min(time_difference), max(time_difference[time_difference < threshold]), 100)
        plt.plot(1e3*x, slope*x + intercept, 'k-')
        plt.gca().annotate(r'$r$ = {:5.3g}'.format(r_value), xy=(0.70, 0.15), xycoords='axes fraction')
        setXLabel(r'\mathrm{time\;to\;previous\;RBC}', r'ms')
        setYLabel(r'S_a - S_v', '')

    def plot(self, x_func, y_func, x_factor=1., y_factor=1., **kwargs):
        """
        Plot the result of the function y_func agains values of the function x_func.

        Args:
            x_func (func): function of the edge index for the x-axis
            y_func (func): function of the edge index for the y-axis
            x_factor (float, optional): multiplication factor for the x-values
            y_factor (float, optional): multiplication factor for the y-values
            **kwargs: optional plotting arguments
        """
        eids_no_type = self.edge_ids()
        for edge_type, edge_method_name in self.edge_type_method.iteritems():
            eids = getattr(self.postprocessor.rbcDataPostProcessor.rbc_path_analyzer, edge_method_name)()
            x_values = x_factor*self.filter_with_n_path_min(np.array([x_func(ei) for ei in eids]), eids)
            y_values = y_factor*self.filter_with_n_path_min(np.array([y_func(ei) for ei in eids]), eids)
            plt.plot(x_values, y_values, linestyle='None',
                     **dict(self.edge_type_symbols[edge_type], **kwargs))
            eids_no_type = [ei for ei in eids_no_type if ei not in eids]

        x_values = x_factor*self.filter_with_n_path_min(np.array([x_func(ei) for ei in eids_no_type]),
                                                        eids_no_type)
        y_values = y_factor*self.filter_with_n_path_min(np.array([y_func(ei) for ei in eids_no_type]),
                                                        eids_no_type)
        plt.plot(x_values, y_values, linestyle='None',
                 **dict(self.edge_type_symbols['noType'], **kwargs))
        plt.gca().ticklabel_format(axis='x', style='sci', scilimits=(-3, 3))
        plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))


    def filter_with_n_path_min(self, array, eids):
        """
        Filter array values on edges based on the minimal required number of paths per edge.

        Args:
            array (np.ndarray): array to filter
            eids (list): edge indices

        Returns:
            np.ndarray, filtered array
        """
        n_sampled_path = np.array([self.postprocessor.n_sampled_paths(ei) for ei in eids])
        return array[n_sampled_path >= self.n_profile_min]


class HemoglobinOnWholePathsPlotter(GraphCasePlotter):
    """
    Plotter object for hemoglobin saturation on whole paths in a capillary network.
    """

    marker_dict = {'marker': 'o',
                   'markeredgecolor': (0, 0, 0, 0.4),
                   'markerfacecolor': (1, 1, 1, 0),
                   'markersize': 3.0}

    distal_hb_marker_dict = {'marker': '.',
                             'markeredgecolor': 'k',
                             'markerfacecolor': (1, 1, 1, 0),
                             'markersize': 1.5,
                             'markeredgewidth': 0.3}

    def __init__(self, postprocessor):
        """
        Constructor

        Args:
            postprocessor (HemoglobinOnWholePathsPostProcessor): postprocessor
        """
        super(HemoglobinOnWholePathsPlotter, self).__init__(postprocessor)
        if not isinstance(postprocessor, HemoglobinOnWholePathsPostProcessor):
            ValueError('The given postprocessor has type {:s}'.format(type(postprocessor)))

    def plot_distal_hb_vs_transit_time_with_paths(self, path_ids, **kwargs):
        ax = plt.gca()
        path_style = kwargs.get('path_style', {})
        symbol_style = kwargs.get('symbol_style', {})
        error_bars = kwargs.get('error_bars', True)
        for path_i in path_ids:
            self.plotFieldWholePathVsTransitTime('Hb_mean', path_i, style=path_style)
        self.plot_distal_hb_vs_transit_time(**symbol_style)
        setYLabel('S', '')
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        x0, x1 = ax.set_xlim(ax.xaxis.get_data_interval())
        ax.set_xlim(x0 - 0.015*(x1 - x0), x1 + 0.05*(x1 - x0))

        if error_bars:
            mean_hb_v = self.postprocessor.mean_distal_hb()
            std_hb_v = self.postprocessor.std_distal_hb()
            mean_tt = self.postprocessor.mean_transit_time()
            std_tt = self.postprocessor.std_transit_time()
            # (see https://matplotlib.org/users/transforms_tutorial.html)
            trans_x = transforms.blended_transform_factory(
                ax.transAxes, ax.transData)
            ax.errorbar(0.97, mean_hb_v, yerr=std_hb_v,
                        transform=trans_x,
                        fmt='o', markersize=2, linewidth=1,
                        elinewidth=0.5,
                        capsize=2, capthick=0.5,
                        color='k')
            # trans_y = transforms.blended_transform_factory(
            #     ax.transData, ax.transAxes)
            y_error = ax.get_ylim()[0] + 0.03*(ax.get_ylim()[1] - ax.get_ylim()[0])
            ax.errorbar(1e3*mean_tt, y_error, xerr=1e3*std_tt,
                        # transform=trans_y,
                        fmt='o', markersize=2, linewidth=1,
                        elinewidth=0.5,
                        capsize=2, capthick=0.5,
                        color='k')

    def plot_distal_hb_vs_transit_time(self, **kwargs):
        self.plot(self.postprocessor.transit_times, self.postprocessor.distal_hb,
                  x_factor=1e3, **dict(self.distal_hb_marker_dict, **kwargs))
        setXLabel('\mathrm{transit\;time}', 'ms')
        setYLabel('\mathrm{outflow\;saturation}\;S_v', '')

    def plot_distal_hb_vs_path_length_with_paths(self, path_ids, **kwargs):
        ax = plt.gca()
        path_style = kwargs.get('path_style', {})
        symbol_style = kwargs.get('symbol_style', {})
        error_bars = kwargs.get('error_bars', True)
        for path_i in path_ids:
            self.plotFieldWholePathVsPathLength('Hb_mean', path_i, style=path_style)
        self.plot_distal_hb_vs_path_length(**symbol_style)
        setYLabel('S', '')
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        x0, x1 = ax.set_xlim(ax.xaxis.get_data_interval())
        ax.set_xlim(x0 - 0.015*(x1 - x0), x1 + 0.05*(x1 - x0))

        if error_bars:
            mean_hb_v = self.postprocessor.mean_distal_hb()
            std_hb_v = self.postprocessor.std_distal_hb()
            mean_tp = self.postprocessor.mean_transit_path_length()
            std_tp = self.postprocessor.std_transit_path_length()
            # (see https://matplotlib.org/users/transforms_tutorial.html)
            trans_x = transforms.blended_transform_factory(
                ax.transAxes, ax.transData)
            ax.errorbar(0.97, mean_hb_v, yerr=std_hb_v,
                        transform=trans_x,
                        fmt='o', markersize=2, linewidth=1,
                        elinewidth=0.5,
                        capsize=2, capthick=0.5,
                        color='k')
            # trans_y = transforms.blended_transform_factory(
            #     ax.transData, ax.transAxes)
            y_error = ax.get_ylim()[0] + 0.03*(ax.get_ylim()[1] - ax.get_ylim()[0])
            ax.errorbar(1e6*mean_tp, y_error, xerr=1e6*std_tp,
                        # transform=trans_y,
                        fmt='o', markersize=2, linewidth=1,
                        elinewidth=0.5,
                        capsize=2, capthick=0.5,
                        color='k')

    def plot_distal_hb_vs_path_length(self, **kwargs):
        self.plot(self.postprocessor.transit_path_lengths, self.postprocessor.distal_hb,
                  x_factor=1e6, **dict(self.distal_hb_marker_dict, **kwargs))
        setXLabel('\mathrm{path\;length}\;s', 'um')
        setYLabel('\mathrm{outflow\;saturation}\;S_v', '')

    def plot_hb_drop_vs_transit_time(self, **kwargs):
        self.plot(self.postprocessor.transit_times, self.postprocessor.hb_drop,
                  x_factor=1e3, **kwargs)
        setXLabel('\mathrm{transit\;time}', 'ms')
        setYLabel('\mathrm{saturation\;drop}\;S_a - S_v', '')

    def plot_hb_drop_vs_proximal_hb(self):
        self.plot(self.postprocessor.proximal_hb, self.postprocessor.hb_drop)
        setXLabel('\mathrm{inflow\;saturation}\;S_a', '')
        setYLabel('\mathrm{saturation\;drop}\;S_a - S_v', '')

    def plot_mean_hb_slope_transit_time_vs_proximal_hb(self):
        self.plot(self.postprocessor.proximal_hb, self.postprocessor.mean_hb_slope_transit_time,
                  y_factor=1e0)
        setXLabel('\mathrm{inflow\;saturation}\;S_a', '')
        setYLabel('(S_a - S_v)/tt', r'\mathrm{s}^{-1}')

    def plot_mean_hb_slope_path_length_vs_proximal_hb(self):
        self.plot(self.postprocessor.proximal_hb, self.postprocessor.mean_hb_slope_path_length,
                  y_factor=1e-3)
        setXLabel('\mathrm{inflow\;saturation}\;S_a', '')
        setYLabel('(S_a - S_v)/ts', r'\mathrm{mm}^{-1}')

    def plot_simulated_hb_distal_distribution(self, bin_width=0.01, **kwargs):
        values = self.postprocessor.distal_hb()
        bounds = rounded_bounds(values, bin_width)
        bins = np.arange(bounds[0], bounds[1]+0.5*bin_width, bin_width)
        plt.hist(values, bins=bins, histtype='bar', **kwargs)
        setXLabel('\mathrm{outflow\;saturation}\;S_v', '')

    def plot_integrated_hb_distal_distribution(self, bin_width=0.01, **kwargs):
        hb, weights = self.postprocessor.integrator.distal_hb_distribution()
        bounds = rounded_bounds(hb, bin_width)
        bins = np.arange(bounds[0], bounds[1]+0.5*bin_width, bin_width)
        plt.hist(hb, bins=bins, histtype='bar', weights=weights, **kwargs)
        setXLabel('\mathrm{outflow\;saturation}\;S_v', '')

    def plot_compared_hb_distal_distribution(self, bin_width=0.02):
        """
        Plot a comparison between the simulated and the integrated hemoglobin distribution.

        Args:
            bin_width (float): bin width for histogram
        """
        ax = plt.gca()
        sim_hb = self.postprocessor.distal_hb()
        # compute integrated distal hemoglobin distributions
        self.postprocessor.integrator.use_topological_radii = True
        self.postprocessor.integrator.compute()
        int_topol_hb, topol_weights = self.postprocessor.integrator.distal_hb_distribution()
        int_topol_mean = self.postprocessor.integrator.average_distal_hb()
        int_topol_std = self.postprocessor.integrator.std_distal_hb()
        self.postprocessor.integrator.use_topological_radii = False
        self.postprocessor.integrator.compute()
        int_func_hb, func_weights = self.postprocessor.integrator.distal_hb_distribution()
        int_func_mean = self.postprocessor.integrator.average_distal_hb()
        int_func_std = self.postprocessor.integrator.std_distal_hb()

        # plot histograms
        bounds = rounded_bounds(np.hstack((sim_hb, int_topol_hb, int_func_hb)), bin_width)
        bins = np.arange(bounds[0], bounds[1]+0.5*bin_width, bin_width)
        bar_fraction = 0.7
        bar_width = 1./3.*bin_width*bar_fraction
        simul_count, _ = np.histogram(sim_hb, bins=bins, normed=True)
        int_topol_count, _ = np.histogram(int_topol_hb, weights=topol_weights,
                                          bins=bins, normed=True)
        int_func_count, _ = np.histogram(int_func_hb, weights=func_weights,
                                          bins=bins, normed=True)
        x1 = 0.5*(bins[:-1] + bins[1:]) - bar_width
        x2 = 0.5*(bins[:-1] + bins[1:])
        x3 = 0.5*(bins[:-1] + bins[1:]) + bar_width
        ax.bar(x1, simul_count, bar_width,
               edgecolor=style_scheme['simul_radii_edge_color'],
               facecolor=style_scheme['simul_radii_face_color'],
               linewidth=0.3,
               label='moving RBCs')
        ax.bar(x2, int_func_count, bar_width,
               edgecolor=style_scheme['int_func_radii_edge_color'],
               facecolor=style_scheme['int_func_radii_face_color'],
               linewidth=0.3,
               label='ODE (functional)')
        ax.bar(x3, int_topol_count, bar_width,
               edgecolor=style_scheme['int_topol_radii_edge_color'],
               facecolor=style_scheme['int_topol_radii_face_color'],
               linewidth=0.3,
               label='ODE (geometric)')
        xlim = max(0.0, round_to_base(ax.get_xlim()[0], 0.1, func=np.floor)), ax.get_xlim()[1]
        ax.set_xlim(xlim)

        # plot mean and standard deviation
        y_max = np.max(np.hstack((simul_count, int_topol_count, int_func_count)))
        sim_mean = np.mean(self.postprocessor.distal_hb())
        sim_std = np.std(self.postprocessor.distal_hb())
        ax.errorbar(sim_mean, 1.1*y_max, xerr=sim_std,
                    fmt='o', markersize=2, linewidth=1,
                    elinewidth=0.5,
                    capsize=2, capthick=0.5,
                    color=style_scheme['simul_radii_face_color'])
        ax.errorbar(int_func_mean, 1.06*y_max, xerr=int_func_std,
                    fmt='o', markersize=2, linewidth=1,
                    elinewidth=0.5,
                    capsize=2, capthick=0.5,
                    color=style_scheme['int_func_radii_face_color'])
        ax.errorbar(int_topol_mean, 1.02*y_max, xerr=int_topol_std,
                    fmt='o', markersize=2, linewidth=1,
                    elinewidth=0.5,
                    capsize=2, capthick=0.5,
                    color=style_scheme['int_topol_radii_edge_color'])
        styles.create_COSH_legend(ax, bbox_to_anchor=(0.45, 0.95))
        setXLabel('\mathrm{outflow\;saturation}\;S_v', '')
        setYLabel('\mathrm{probability\;density}', '')
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(5))

    def plot_compared_hb_distal_distribution_multipanel(self, bin_width=0.02,
                                                        panel_annotation=('A', 'B', 'C')):
        """
        Plot a comparison between the simulated and the integrated hemoglobin distribution.

        Args:
            bin_width (float): bin width for histogram
            panel_annotation (sequence): panel annotations
        """
        f, axarr = plt.subplots(4, sharex=True)
        hist_axes = axarr[0:3]
        bp_ax = axarr[-1]
        sim_hb = self.postprocessor.distal_hb()
        # compute integrated distal hemoglobin distributions
        self.postprocessor.integrator.use_topological_radii = True
        self.postprocessor.integrator.compute()
        int_topol_hb, topol_weights = self.postprocessor.integrator.distal_hb_distribution()
        int_topol_hb_repeat = self.postprocessor.integrator.distal_hb_distribution_with_repeat()
        int_topol_mean = self.postprocessor.integrator.average_distal_hb()
        int_topol_std = self.postprocessor.integrator.std_distal_hb()
        self.postprocessor.integrator.use_topological_radii = False
        self.postprocessor.integrator.compute()
        int_func_hb, func_weights = self.postprocessor.integrator.distal_hb_distribution()
        int_func_hb_repeat = self.postprocessor.integrator.distal_hb_distribution_with_repeat()
        int_func_mean = self.postprocessor.integrator.average_distal_hb()
        int_func_std = self.postprocessor.integrator.std_distal_hb()

        # plot histograms
        bounds = rounded_bounds(np.hstack((sim_hb, int_topol_hb, int_func_hb)), bin_width)
        bins = np.arange(bounds[0], bounds[1]+0.5*bin_width, bin_width)
        simul_count, _ = np.histogram(sim_hb, bins=bins, normed=True)
        int_topol_count, _ = np.histogram(int_topol_hb, weights=topol_weights,
                                          bins=bins, normed=True)
        int_func_count, _ = np.histogram(int_func_hb, weights=func_weights,
                                         bins=bins, normed=True)
        bar_fraction = 0.8
        bar_width = bin_width*bar_fraction
        x = 0.5*(bins[:-1] + bins[1:])
        hist_axes[0].bar(x, simul_count, bar_width,
                         edgecolor=style_scheme['simul_radii_edge_color'],
                         facecolor=style_scheme['simul_radii_face_color'],
                         linewidth=0.3,
                         label='moving RBCs')
        hist_axes[1].bar(x, int_func_count, bar_width,
                         edgecolor=style_scheme['int_func_radii_edge_color'],
                         facecolor=style_scheme['int_func_radii_face_color'],
                         linewidth=0.3,
                         label='ODE (functional)')
        hist_axes[2].bar(x, int_topol_count, bar_width,
                         edgecolor=style_scheme['int_topol_radii_edge_color'],
                         facecolor=style_scheme['int_topol_radii_face_color'],
                         linewidth=0.3,
                         label='ODE (geometric)')
        xlim = max(0.0, round_to_base(axarr[0].get_xlim()[0], 0.1, func=np.floor)), \
               axarr[0].get_xlim()[1]
        axarr[0].set_xlim(xlim)

        # plot mean and standard deviation
        y_max = np.max(np.hstack((simul_count, int_topol_count, int_func_count)))
        sim_mean = np.mean(self.postprocessor.distal_hb())
        sim_std = np.std(self.postprocessor.distal_hb())
        y_error_bar = 1.05*y_max
        hist_axes[0].errorbar(sim_mean, y_error_bar, xerr=sim_std,
                              fmt='o', markersize=2, linewidth=1,
                              elinewidth=0.5,
                              capsize=2, capthick=0.5,
                              color=style_scheme['simul_radii_face_color'])
        hist_axes[1].errorbar(int_func_mean, y_error_bar, xerr=int_func_std,
                              fmt='o', markersize=2, linewidth=1,
                              elinewidth=0.5,
                              capsize=2, capthick=0.5,
                              color=style_scheme['int_func_radii_face_color'])
        hist_axes[2].errorbar(int_topol_mean, y_error_bar, xerr=int_topol_std,
                              fmt='o', markersize=2, linewidth=1,
                              elinewidth=0.5,
                              capsize=2, capthick=0.5,
                              color=style_scheme['int_topol_radii_edge_color'])
        setXLabel('\mathrm{outflow\;saturation}\;S_v', '')
        axarr[0].xaxis.set_major_locator(MultipleLocator(0.1))
        for ax, annotation in zip(hist_axes, panel_annotation):
            styles.create_COSH_legend(ax, loc='upper left')
            ax.yaxis.set_major_locator(MultipleLocator(5))
            # ax.yaxis.set_major_locator(MultipleLocator(10))
            setYLabel('f_{S_v}', '', ax=ax)
        for ax, annotation in zip(axarr, panel_annotation):
            annotate_axis_corner(ax, annotation, xy=(0.956, 0.79))
            # if annotation != 'C':
            #     annotate_axis_corner(ax, annotation, xy=(0.956, 0.79))
            # else:
            #     annotate_axis_corner(ax, annotation, xy=(0.956, 0.72))

        # boxplot
        data = [int_topol_hb_repeat, int_func_hb_repeat, sim_hb]
        offset = len(panel_annotation) - len(data)
        labels = panel_annotation[-1-offset::-1]
        medianprops = {'color': 'w'}
        flierprops = {'marker': '.', 'markersize': 2}
        bp = bp_ax.boxplot(data, whis=[5, 95], vert=False, labels=labels,
                           patch_artist=True, medianprops=medianprops, flierprops=flierprops,
                           showfliers=False)
        fill_colors = [style_scheme['int_topol_radii_edge_color'],
                       style_scheme['int_func_radii_face_color'],
                       style_scheme['simul_radii_face_color']]
        for patch, color in zip(bp['boxes'], fill_colors):
            patch.set_facecolor(color)
        print "5th percentile of topological:  ", np.percentile(int_topol_hb_repeat, 5)
        print "First quartile of topological:  ", np.percentile(int_topol_hb_repeat, 25)
        print "Median of topological:          ", np.percentile(int_topol_hb_repeat, 50)
        print "Third quartile of topological:  ", np.percentile(int_topol_hb_repeat, 75)
        print "95th percentile of topological: ", np.percentile(int_topol_hb_repeat, 95)

        # statistical testing
        F = int_topol_std**2/sim_std**2
        pvalue = scipy.stats.f.sf(F, np.size(int_topol_hb) - 1, np.size(sim_hb) - 1)
        print "p-value for ODE model (topological radii) vs. moving RBC: {:g}".format(pvalue)
        F = int_func_std**2/sim_std**2
        pvalue = scipy.stats.f.sf(F, np.size(int_func_hb) - 1, np.size(sim_hb) - 1)
        print "p-value for ODE model (functional radii) vs. moving RBC: {:g}".format(pvalue)
        F = int_topol_std**2/int_func_std**2
        pvalue = scipy.stats.f.sf(F, np.size(int_topol_hb) - 1, np.size(int_func_hb) - 1)
        print "p-value for ODE model, topological vs. functional radii: {:g}".format(pvalue)

    def plot_integrated_hb_profile_vs_path_length(self, **kwargs):
        self.plot_integrated_hb_profile('path_coord', x_factor=1e6, **kwargs)

    def plot_integrated_hb_profile_vs_transit_time(self, **kwargs):
        self.plot_integrated_hb_profile('transit_time', x_factor=1e3, **kwargs)

    def plot_integrated_hb_profile(self, attribute_suffix, x_factor=1., **kwargs):
        style = kwargs.get('style', {'color': 'k'})
        mean_linewidth = kwargs.get('mean_linewidth', 2)
        graph = self.postprocessor.integrator.graph
        mean_flow = np.mean(graph.es()['rbc_flow'])
        for edge in graph.es():
            for int_hb in edge['hb_distribution']:
                linewidth = int_hb.flow_part*mean_linewidth/mean_flow
                plt.plot([x_factor*getattr(int_hb, 'upstream_' + attribute_suffix),
                          x_factor*getattr(int_hb, 'downstream_' + attribute_suffix)],
                         [int_hb.upstream_hb, int_hb.downstream_hb], '-',
                         linewidth=linewidth, **style)

    def plot(self, x_func, y_func, x_factor=1., y_factor=1., **kwargs):
        """
        Plot the result of the function y_func agains values of the function x_func.

        Args:
            x_func (func): function for the x-axis
            y_func (func): function for the y-axis
            x_factor (float, optional): multiplication factor for the x-values
            y_factor (float, optional): multiplication factor for the y-values
            **kwargs: optional plotting arguments
        """
        x_values = x_factor*x_func()
        y_values = y_factor*y_func()
        plt.plot(x_values, y_values, linestyle='None', **dict(self.marker_dict, **kwargs))
