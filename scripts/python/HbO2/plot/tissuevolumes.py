"""
Plotter objects for hemoglobin saturation in a capillary network.
"""

from __future__ import division
import copy
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np

from HbO2.plot import styles
from HbO2.plot.labels import setXLabel, setYLabel
from plot.utils import annotate_axis_corner

style_scheme = styles.StyleSchemeCOSH()


class TissueVolumePlotter(object):
    """
    Plotter object for tissue volumes in a capillary network.
    """

    volume_type_symbols = {'functional': {'marker': '^',
                                          'markeredgecolor': style_scheme['int_func_radii_face_color'],
                                          'markerfacecolor': (1, 1, 1, 0)},
                           'topological': {'marker': 'o',
                                           'markersize': 2.5,
                                           'markeredgecolor': style_scheme['int_topol_radii_marker_edge_color'],
                                           'markerfacecolor': style_scheme['int_topol_radii_marker_face_color'],
                                           'alpha': 0.5},
                           'difference': {'marker': 's',
                                          'markeredgecolor': 'k',
                                          'markerfacecolor': (1, 1, 1, 0)}}

    def __init__(self, results):
        """
        Constructor

        Args:
            results: simulation results
        """
        super(TissueVolumePlotter, self).__init__()
        self.results = results

    def plot_histogram_tissue_radii(self, **kwargs):
        style_scheme = styles.StyleSchemeCOSH()
        x_range = (0, 40)
        bin_width = 2
        bar_fraction = 0.7
        bar_width = 0.5*bin_width*bar_fraction
        bins = np.array(np.linspace(x_range[0], x_range[1], (x_range[1] - x_range[0])/bin_width + 1))
        ax = plt.gca()
        func_radii = 1e6*np.array(self.results.functional_tissue_radii())
        topol_radii = 1e6*np.array(self.results.topological_tissue_radii())
        functional_count, _ = np.histogram(func_radii, bins=bins)
        topological_count, _ = np.histogram(topol_radii, bins=bins)
        x1 = 0.5*(bins[:-1] + bins[1:]) - 0.5*bar_width
        x2 = 0.5*(bins[:-1] + bins[1:]) + 0.5*bar_width
        plt.bar(x1, functional_count, bar_width,
                edgecolor=style_scheme['int_func_radii_edge_color'],
                facecolor=style_scheme['int_func_radii_face_color'],
                linewidth=0.4,
                label='functional')
        plt.bar(x2, topological_count, bar_width,
                edgecolor=style_scheme['int_topol_radii_edge_color'],
                facecolor=style_scheme['int_topol_radii_face_color'],
                linewidth=0.4,
                label='geometric')

        # plot mean and standard deviation
        y_max = np.max(np.hstack((functional_count, topological_count)))
        ax.errorbar(np.mean(func_radii), 1.1*y_max, xerr=np.std(func_radii),
                    fmt='o', markersize=2, linewidth=1,
                    elinewidth=0.5,
                    capsize=2, capthick=0.5,
                    color=style_scheme['int_func_radii_face_color'])
        ax.errorbar(np.mean(topol_radii), 1.05*y_max, xerr=np.std(topol_radii),
                    fmt='o', markersize=2, linewidth=1,
                    elinewidth=0.5,
                    capsize=2, capthick=0.5,
                    color=style_scheme['int_topol_radii_edge_color'])

        ax.set_xlim(x_range)
        majorFormatter = FormatStrFormatter('%d')
        ax.yaxis.set_major_formatter(majorFormatter)
        styles.create_COSH_legend(ax, bbox_to_anchor=(1, 0.90))
        setXLabel('\mathrm{tissue\;radius}\;r_t', 'um')
        setYLabel('\mathrm{frequency}', '')

    def plot_functional_vs_topological_tissue_radii(self, **kwargs):
        self.plot(self.results.topological_tissue_radii,
                  self.results.functional_tissue_radii,
                  x_factor=1e6, y_factor=1e6, **dict(self.volume_type_symbols['difference'], **kwargs))
        setXLabel('r_{t,\mathrm{geom}}', 'um')
        setYLabel('r_{t,\mathrm{func}}', 'um')

    def plot_tissue_radii_vs_hb_upstream(self, **kwargs):
        self.plot(self.results.hb_mean_upstream, self.results.topological_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['topological'], **kwargs))
        self.plot(self.results.hb_mean_upstream, self.results.functional_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['functional'], **kwargs))
        setXLabel('S_a', '')
        setYLabel('\mathrm{tissue\;radius}\;r_t', 'um')

    def plot_tissue_radii_vs_hb_downstream(self, **kwargs):
        self.plot(self.results.hb_mean_downstream, self.results.topological_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['topological'], **kwargs))
        self.plot(self.results.hb_mean_downstream, self.results.functional_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['functional'], **kwargs))
        setXLabel('S_v', '')
        setYLabel('\mathrm{tissue\;radius}\;r_t', 'um')

    def plot_tissue_radii_vs_hb_mean(self, **kwargs):
        kwargs_topol, kwargs_func = self.append_label_with_geom_and_func(**kwargs)
        self.plot(self.results.hb_mean_on_edge, self.results.topological_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['topological'], **kwargs_topol))
        self.plot(self.results.hb_mean_on_edge, self.results.functional_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['functional'], **kwargs_func))
        setXLabel(r'\mathrm{average\;saturation}\;S', '')
        setYLabel('\mathrm{tissue\;radius}\;r_t', 'um')

    def append_label_with_geom_and_func(self, **kwargs):
        kwargs_topol = copy.deepcopy(kwargs)
        kwargs_func = copy.deepcopy(kwargs)
        if 'label' in kwargs:
            label_text = kwargs['label']
            label_topological = label_text + ' (geometric)'
            label_functional = label_text + ' (functional)'
            kwargs_topol.update(label=label_topological)
            kwargs_func.update(label=label_functional)
        return kwargs_topol, kwargs_func

    def plot_tissue_radii_vs_hb_drop(self, **kwargs):
        self.plot(self.results.hb_mean_drop, self.results.topological_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['topological'], **kwargs))
        self.plot(self.results.hb_mean_drop, self.results.functional_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['functional'], **kwargs))
        setXLabel(r'S_a - S_v', '')
        setYLabel('\mathrm{tissue\;radius}\;r_t', 'um')

    def plot_tissue_radii_vs_rbc_flow(self, **kwargs):
        kwargs_topol, kwargs_func = self.append_label_with_geom_and_func(**kwargs)
        self.plot(self.results.rbc_flow, self.results.topological_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['topological'], **kwargs_topol))
        self.plot(self.results.rbc_flow, self.results.functional_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['functional'], **kwargs_func))
        setXLabel('\mathrm{RBC\;flow}', '\mathrm{s}^{-1}')
        setYLabel('\mathrm{tissue\;radius}\;r_t', 'um')

    def plot_tissue_radii_vs_rbc_velocity(self, **kwargs):
        self.plot(self.results.rbc_velocity, self.results.topological_tissue_radii,
                  x_factor=1e3, y_factor=1e6,
                  **dict(self.volume_type_symbols['topological'], **kwargs))
        self.plot(self.results.rbc_velocity, self.results.functional_tissue_radii,
                  x_factor=1e3, y_factor=1e6,
                  **dict(self.volume_type_symbols['functional'], **kwargs))
        setXLabel('\mathrm{RBC\;velocity}\;v_{\mathrm{rbc}}', '\mathrm{mm/s}')
        setYLabel('\mathrm{tissue\;radius}\;r_t', 'um')

    def plot_tissue_radii_vs_hematocrit(self, **kwargs):
        self.plot(self.results.hematocrit, self.results.topological_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['topological'], **kwargs))
        self.plot(self.results.hematocrit, self.results.functional_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['functional'], **kwargs))
        setXLabel('\mathrm{tube\;hematocrit}\;H_T', '')
        setYLabel('\mathrm{tissue\;radius}\;r_t', 'um')

    def plot_tissue_radii_vs_arrival_transit_time(self, **kwargs):
        self.plot(self.results.mean_arrival_transit_time, self.results.topological_tissue_radii,
                  x_factor=1e3, y_factor=1e6,
                  **dict(self.volume_type_symbols['topological'], **kwargs))
        self.plot(self.results.mean_arrival_transit_time, self.results.functional_tissue_radii,
                  x_factor=1e3, y_factor=1e6,
                  **dict(self.volume_type_symbols['functional'], **kwargs))
        setXLabel('\mathrm{arrival\;transit\;time}\;tt', 'ms')
        setYLabel('\mathrm{tissue\;radius}\;r_t', 'um')

    def plot_tissue_radii_vs_arrival_path_length(self, **kwargs):
        self.plot(self.results.mean_arrival_path_length, self.results.topological_tissue_radii,
                  x_factor=1e6, y_factor=1e6,
                  **dict(self.volume_type_symbols['topological'], **kwargs))
        self.plot(self.results.mean_arrival_path_length, self.results.functional_tissue_radii,
                  x_factor=1e6, y_factor=1e6,
                  **dict(self.volume_type_symbols['functional'], **kwargs))
        setXLabel('\mathrm{arrival\;path\;length}\;ts', 'um')
        setYLabel('\mathrm{tissue\;radius}\;r_t', 'um')

    def plot_tissue_radii_vs_edge_index(self, **kwargs):
        self.plot(self.results.edge_ids, self.results.topological_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['topological'], **kwargs))
        self.plot(self.results.edge_ids, self.results.functional_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['functional'], **kwargs))
        setXLabel('e_i', '')
        setYLabel('\mathrm{tissue\;radius}\;r_t', 'um')

    def plot_tissue_radius_difference_vs_hb_mean(self, **kwargs):
        self.plot(self.results.hb_mean_on_edge,
                  self.results.functional_minus_topological_tissue_radii,
                  y_factor=1e6, **dict(self.volume_type_symbols['difference'], **kwargs))
        setXLabel(r'\mathrm{average\;saturation}\;S', '')
        setYLabel('r_{t,\mathrm{func}} - r_{t,\mathrm{topol}}', 'um')

    def plot(self, x_func, y_func, x_factor=1., y_factor=1., **kwargs):
        """
        Plot the result of the function y_func against values of the function x_func.

        Args:
            x_func (func): function of the edge index for the x-axis
            y_func (func): function of the edge index for the y-axis
            x_factor (float, optional): multiplication factor for the x-values
            y_factor (float, optional): multiplication factor for the y-values
            **kwargs: optional plotting arguments
        """
        x_values = x_factor*np.asarray(x_func())
        y_values = y_factor*np.asarray(y_func())
        plt.plot(x_values, y_values, linestyle='None', **kwargs)


class TissueVolumeCombinedPlotter(object):
    """
    Plotter object for different non-pooled simulations for tissue volumes in a capillary
    network.
    """

    def __init__(self, results_list):
        """
        Constructor

        Args:
            results (list): simulation results
        """
        super(TissueVolumeCombinedPlotter, self).__init__()
        self.results_list = results_list

    def plot_histogram_tissue_radii(self, functional_colors, functional_labels, **kwargs):
        ax = plt.gca()
        style_scheme = styles.StyleSchemeCOSH()
        n_bar_types = 1 + len(self.results_list)
        x_range = (0, 40)
        bin_width = 2
        bar_fraction = 0.7
        bar_width = bin_width*bar_fraction/n_bar_types
        bins = np.array(np.linspace(x_range[0], x_range[1], (x_range[1] - x_range[0])/bin_width + 1))
        bin_center = 0.5*(bins[:-1] + bins[1:])
        y_max = 0.0

        for i, results in enumerate(self.results_list):
            func_radii = 1e6*np.array(results.functional_tissue_radii())
            functional_count, _ = np.histogram(func_radii, bins=bins)
            y_max = max(y_max, np.max(functional_count))
            offset = (-0.5 + (i + 0.5)/n_bar_types)*bin_width*bar_fraction
            x_func = bin_center + offset
            plt.bar(x_func, functional_count, bar_width,
                    edgecolor='k',
                    facecolor=functional_colors[i],
                    linewidth=0.4,
                    label=functional_labels[i])

        offset = (0.5 - 0.5/n_bar_types)*bin_width*bar_fraction
        x_topol = bin_center + offset
        topol_radii = 1e6*np.array(self.results_list[0].topological_tissue_radii())
        topological_count, _ = np.histogram(topol_radii, bins=bins)
        plt.bar(x_topol, topological_count, bar_width,
                edgecolor=style_scheme['int_topol_radii_edge_color'],
                facecolor=style_scheme['int_topol_radii_face_color'],
                linewidth=0.4,
                label='geometric')

        # plot mean and standard deviation
        y_max = max(y_max, np.max(topological_count))
        for i, results in enumerate(self.results_list):
            func_radii = 1e6*np.array(results.functional_tissue_radii())
            y_error_bar = (1.1 - i/(len(self.results_list))*0.07)*y_max
            ax.errorbar(np.mean(func_radii), y_error_bar, xerr=np.std(func_radii),
                        fmt='o', markersize=2, linewidth=1,
                        elinewidth=0.5,
                        capsize=2, capthick=0.5,
                        color=functional_colors[i])

        y_error_bar = 1.03*y_max
        ax.errorbar(np.mean(topol_radii), y_error_bar, xerr=np.std(topol_radii),
                    fmt='o', markersize=2, linewidth=1,
                    elinewidth=0.5,
                    capsize=2, capthick=0.5,
                    color=style_scheme['int_topol_radii_edge_color'])

        ax.set_xlim(x_range)
        styles.create_COSH_legend(ax, bbox_to_anchor=(0.48, 0.85),
                                  handlelength=2.3, handleheight=0.4)
        setXLabel('\mathrm{tissue\;radius}\;r_t', 'um')
        setYLabel('\mathrm{frequency}', '')

    def plot_histogram_tissue_radii_multipanel(self, functional_colors, functional_labels, **kwargs):
        f, axarr = plt.subplots(len(self.results_list)+2, sharex=True)
        panel_annotation = 'ABCD'
        style_scheme = styles.StyleSchemeCOSH()
        x_range = (0, 40)
        bin_width = 2
        bar_fraction = 0.7
        bar_width = bin_width*bar_fraction
        bins = np.array(np.linspace(x_range[0], x_range[1], (x_range[1] - x_range[0])/bin_width + 1))
        bin_center = 0.5*(bins[:-1] + bins[1:])
        y_max = 0.0
        hist_axes = axarr[0:-1]

        # boxplot
        data = [1e6*np.array(self.results_list[0].topological_tissue_radii()),
                1e6*np.array(self.results_list[1].functional_tissue_radii()),
                1e6*np.array(self.results_list[0].functional_tissue_radii())]
        labels = 'CBA'
        medianprops = {'color': 'w'}
        flierprops = {'marker': '.', 'markersize': 2}
        bp = axarr[-1].boxplot(data, whis=[5, 95], vert=False, labels=labels,
                               patch_artist=True, medianprops=medianprops,
                               showfliers=False, flierprops=flierprops)
        fill_colors = [style_scheme['int_topol_radii_face_color'],
                       functional_colors[1],
                       functional_colors[0]]
        for patch, color in zip(bp['boxes'], fill_colors):
            patch.set_facecolor(color)

        # histograms
        for i, (results, ax) in enumerate(zip(self.results_list, hist_axes)):
            func_radii = 1e6*np.array(results.functional_tissue_radii())
            functional_count, _ = np.histogram(func_radii, bins=bins)
            y_max = max(y_max, np.max(functional_count))
            x_func = bin_center
            ax.bar(x_func, functional_count, bar_width,
                   edgecolor='k',
                   facecolor=functional_colors[i],
                   linewidth=0.4,
                   label=functional_labels[i])

        x_topol = bin_center
        topol_radii = 1e6*np.array(self.results_list[0].topological_tissue_radii())
        topological_count, _ = np.histogram(topol_radii, bins=bins)
        axarr[-2].bar(x_topol, topological_count, bar_width,
                      edgecolor=style_scheme['int_topol_radii_edge_color'],
                      facecolor=style_scheme['int_topol_radii_face_color'],
                      linewidth=0.4,
                      label='geometric')
        y_max = max(y_max, np.max(topological_count))

        # plot mean and standard deviation
        y_error_bar = 1.04*y_max
        for i, (results, ax) in enumerate(zip(self.results_list, hist_axes)):
            func_radii = 1e6*np.array(results.functional_tissue_radii())
            ax.errorbar(np.mean(func_radii), y_error_bar, xerr=np.std(func_radii),
                        fmt='o', markersize=2, linewidth=1,
                        elinewidth=0.5,
                        capsize=2, capthick=0.5,
                        color=functional_colors[i])

        axarr[-2].errorbar(np.mean(topol_radii), y_error_bar, xerr=np.std(topol_radii),
                           fmt='o', markersize=2, linewidth=1,
                           elinewidth=0.5,
                           capsize=2, capthick=0.5,
                           color=style_scheme['int_topol_radii_edge_color'])

        axarr[0].set_xlim(x_range[0] - 1, x_range[1] + 1)
        y_lim_max = 1.1*y_max
        setXLabel('\mathrm{tissue\;radius}\;r_t', 'um')
        for ax in hist_axes:
            majorLocator = MultipleLocator(10)
            ax.yaxis.set_major_locator(majorLocator)
            ax.set_ylim((0, y_lim_max))
            styles.create_COSH_legend(ax, loc='upper left', bbox_to_anchor=(0.01, 0.92))
            setYLabel('\mathrm{frequency}', '', ax=ax)
        for ax, annotation in zip(axarr, panel_annotation):
            annotate_axis_corner(ax, annotation, xy=(0.94, 0.82))
