"""
Module for plotting simulated results produced by OpenFOAM.
"""

import warnings

import matplotlib.pyplot as plt
import numpy as np
from HbO2.postprocess.rbcdata.postprocess import FlowReversalError

from HbO2.plot import labels
from HbO2.postprocess.case import CasePostProcessor


class AxisymmetricCasePlotter(object):
    """
    Plot simulation results of an axisymmetric OpenFOAM case.

    Attributes:
        postprocessor (CasePostProcessor): case postprocessor for graph cases
        sim_params (HbO2SimulationParameters): simulation parameters
        n_points_on_edge (int): number of points for plotting on each edge
    """

    def __init__(self, postprocessor):
        """
        Constructor

        Args:
            postprocessor (CasePostProcessor): postprocessor
        """
        self.postprocessor = postprocessor
        self.sim_params = postprocessor.simParams
        self.n_points_on_edge = 101

    def x_coords(self):
        """
        Return the x-coordinates that are used for computing field values

        Returns:
            np.ndarray, positions from 0 to the domain length
        """
        return np.linspace(0, self.sim_params['domainLength'], self.n_points_on_edge)

    def x_plot(self):
        """
        Return the x-coordinates that are used for plotting

        Returns:
            np.ndarray, positions for plotting, from 0 to the domain length
        """
        return 1e6*self.x_coords()

    def plot_field_profiles(self, field_name, n_profiles, **kwargs):
        """
        Plot a profile of a field average.

        Args:
            field_name (str): name of the field to plot
            n_profiles (int): number of profiles to plot
            **style (dict): line style for plotting
        """
        rbcDataPostProcessor = self.postprocessor.rbcDataPostProcessor
        yValues = rbcDataPostProcessor.fieldOnLastPaths(field_name, self.x_coords(), n_profiles)
        self.plot_spatial_profile(yValues, field_name, **kwargs)

    def plot_field_average_profile(self, field_name, n_average, **kwargs):
        """
        Plot a profile of a field average.

        Args:
            field_name (str): name of the field to plot
            n_average (int): number of RBCs to use for averaging
            **style (dict): line style for plotting
        """
        rbcDataPostProcessor = self.postprocessor.rbcDataPostProcessor
        yValues = rbcDataPostProcessor.fieldAverage(field_name, self.x_coords(), nAverage=n_average)
        self.plot_spatial_profile(yValues, field_name, **kwargs)

    def plot_field_std_profile(self, field_name, n_average, **kwargs):
        """
        Plot a profile of a field average.

        Args:
            field_name (str): name of the field to plot
            n_average (int): number of RBCs to use for averaging
            **style (dict): line style for plotting
        """
        rbcDataPostProcessor = self.postprocessor.rbcDataPostProcessor
        yValues = rbcDataPostProcessor.fieldStd(field_name, self.x_coords(), nAverage=n_average)
        self.plot_spatial_profile(yValues, field_name, **kwargs)

    def plot_field_average_with_std_profile(self, field_name, n_average, **kwargs):
        """
        Plot a profile of a field average plus/minus the standard deviation.

        Args:
            field_name (str): name of the field to plot
            n_average (int): number of RBCs to use for averaging
            **meanStyle (dict): line style for plotting the mean
            **stdStyle (dict): line style for plotting the standard deviation
        """
        meanStyle = kwargs.get('meanStyle', {'style': '-'})
        stdStyle = kwargs.get('stdStyle', {'style': '--'})
        rbcDataPostProcessor = self.postprocessor.rbcDataPostProcessor
        mean = rbcDataPostProcessor.fieldAverage(field_name, self.x_coords(), nAverage=n_average)
        std = rbcDataPostProcessor.fieldStd(field_name, self.x_coords(), nAverage=n_average)
        self.plot_spatial_profile(mean, field_name, style=meanStyle)
        self.plot_spatial_profile(mean + std, field_name, style=stdStyle)
        self.plot_spatial_profile(mean - std, field_name, style=stdStyle)

    def plot_spatial_profile(self, y, ylabel_name, **kwargs):
        """
        Plot a spatial profile of a field with a given label.

        Args:
            y (np.ndarray): values to plot
            ylabel_name (str): string for y-axis (without unit)
            **style (dict): line style for plotting
        """
        style = kwargs.get('style', {'linestyle': '-'})
        if len(y.shape) == 0 or y.shape[0] == 0:
            warnings.warn('No field value returned for variable {:s}, skipping...'\
                          .format(ylabel_name))
            return
        plt.plot(self.x_plot(), y, **style)
        labels.setXLabel('x', 'um')
        paramUnit = labels.openFOAMVarNameToUnit(ylabel_name)
        plottingParamName = labels.openFOAMVarNameToLatex(ylabel_name)
        labels.setYLabel(plottingParamName, paramUnit)


class GraphCasePlotter(object):
    """
    Plot simulation results of OpenFOAM cases with graphs

    Attributes:
        postprocessor (GraphCasePostProcessor): case postprocessor for graph cases
        simParams (HbO2SimulationParameters): simulation parameters
        nPointsOnEdge (int): number of points for plotting on each edge
        enforce_positive_flow (bool): whether to reverse plotting coordinates if flow
                                          direction is negative
    """

    def __init__(self, postprocessor, enforce_positive_flow=False):
        """
        Constructor

        Args:
            postprocessor (GraphCasePostProcessor): postprocessor
        """
        self.postprocessor = postprocessor
        self.simParams = postprocessor.simParams
        self.nPointsOnEdge = 101
        self.enforce_positive_flow = enforce_positive_flow

    def xPlot(self, ei):
        """
        Return the x-coordinates that are used for plotting on edge ei.

        Args:
            ei (int): edge index

        Returns:
            np.ndarray, positions for plotting, from 0 to the edge length
        """
        return self.sToXPlot(self.sCoords(ei), ei)

    def sToXPlot(self, s, ei):
        if self.enforce_positive_flow:
            try:
                positive_flow = self.postprocessor.rbcDataPostProcessor. \
                    rbc_path_analyzer.positive_flow(ei)
            except FlowReversalError:
                print "Flow reversal in edge {:d}, setting positive flow for plotting.".format(ei)
                positive_flow = True
            if not positive_flow:
                s = s[-1] + s[0] - s
        return 1e6*(s - self.simParams.sCoordOffset())

    def sCoords(self, ei):
        """
        Return the s-coordinates that are used for plotting on edge ei.

        Args:
            ei (int): edge index

        Returns:
            np.ndarray, positions for plotting, from 0 to the edge length
        """
        interval = self.postprocessor.coordinate_interval_in_box(ei)
        return np.asarray(np.linspace(interval[0], interval[1], self.nPointsOnEdge))

    def plotFieldProfiles(self, fieldName, edge, nProfiles, **kwargs):
        """
        Plot a profile of a field along a given edge.

        Args:
            fieldName (str): name of the field to plot
            edge (int): edge index in RBC paths
            nProfiles (int): number of profiles to plot
            **style (dict): line style for plotting
        """
        rbcDataPostProcessor = self.postprocessor.rbcDataPostProcessor
        sCoords = self.sCoords(edge)
        yValues = rbcDataPostProcessor.fieldOnEdge(fieldName, sCoords,
                                                   edge, nProfiles)
        self.plotSpatialProfile(yValues, fieldName, edge, **kwargs)
        # RBCIndices = self.postprocessor.rbcDataPostProcessor.rbc_path_analyzer.last_complete_path_indices_on_edge(
        #                 edge, nProfiles)
        # for RBCI in RBCIndices:
        #     x = self.postprocessor.rbc_data.pathCentersOnEdge(RBCI, edge)
        #     y = self.postprocessor.rbc_data.pathFieldOnEdge(RBCI, edge, fieldName)
        #     plt.plot(1e6*x, y, 'b--', linewidth=0.3, alpha=0.5)

    def plotFieldWholePathVsPathLength(self, fieldName, pathI, **kwargs):
        """
        Plot a profile of a field along a given edge.

        Args:
            fieldName (str): name of the field to plot
            pathI (int): index of the RBC path to plot
            **style (dict): line style for plotting
        """
        style = kwargs.get('style', {'linestyle': '-',
                                     'linewidth': 0.2,
                                     'color': (0, 0, 0, 0.2)})
        x_list = self.postprocessor.path_lengths_on_path(pathI)
        y_list = self.postprocessor.field_on_whole_path(fieldName, pathI)
        for x, y in zip(x_list, y_list):
            plt.plot(1e6*x, y, **style)
        labels.setXLabel('x', 'um')
        paramUnit = labels.openFOAMVarNameToUnit(fieldName)
        plottingParamName = labels.openFOAMVarNameToLatex(fieldName)
        labels.setYLabel(plottingParamName, paramUnit)

    def plotFieldWholePathVsTransitTime(self, fieldName, pathI, **kwargs):
        """
        Plot a profile of a field along a given edge.

        Args:
            fieldName (str): name of the field to plot
            pathI (int): index of the RBC path to plot
            **style (dict): line style for plotting
        """
        style = kwargs.get('style', {'linestyle': '-',
                                     'linewidth': 0.2,
                                     'color': (0, 0, 0, 0.2)})
        x_list = self.postprocessor.transit_times_on_path(pathI)
        y_list = self.postprocessor.field_on_whole_path(fieldName, pathI)
        for x, y in zip(x_list, y_list):
            plt.plot(1e3*x, y, **style)
        labels.setXLabel('t', 'ms')
        paramUnit = labels.openFOAMVarNameToUnit(fieldName)
        plottingParamName = labels.openFOAMVarNameToLatex(fieldName)
        labels.setYLabel(plottingParamName, paramUnit)

    def plotFieldAverageProfile(self, fieldName, edge, nAverage, **kwargs):
        """
        Plot a profile of a field average along a given edge.

        Args:
            fieldName (str): name of the field to plot
            edge (int): edge index in RBC paths
            nAverage (int): number of RBCs to use for averaging
            **style (dict): line style for plotting
        """
        rbcDataPostProcessor = self.postprocessor.rbcDataPostProcessor
        sCoords = self.sCoords(edge)
        yValues = rbcDataPostProcessor.fieldAverageOnEdge(fieldName, sCoords,
                                                          edge, nAverage=nAverage)
        self.plotSpatialProfile(yValues, fieldName, edge, **kwargs)

    def plotFieldStdProfile(self, fieldName, edge, nAverage, **kwargs):
        """
        Plot a profile of a field average along a given edge.

        Args:
            fieldName (str): name of the field to plot
            edge (int): edge index in RBC paths
            nAverage (int): number of RBCs to use for averaging
            **style (dict): line style for plotting
        """
        rbcDataPostProcessor = self.postprocessor.rbcDataPostProcessor
        sCoords = self.sCoords(edge)
        yValues = rbcDataPostProcessor.fieldStdOnEdge(fieldName, sCoords,
                                                      edge, nAverage=nAverage)
        self.plotSpatialProfile(yValues, fieldName, edge, **kwargs)

    def plotFieldAverageWithStdProfile(self, fieldName, edge, nAverage, **kwargs):
        """
        Plot a profile of a field average plus/minus the standard deviation
        along a given edge.

        Args:
            fieldName (str): name of the field to plot
            edge (int): edge index in RBC paths
            nAverage (int): number of RBCs to use for averaging
            **meanStyle (dict): line style for plotting the mean
            **stdStyle (dict): line style for plotting the standard deviation
        """
        meanStyle = kwargs.get('meanStyle', {'style': '-'})
        stdStyle = kwargs.get('stdStyle', {'style': '--'})
        rbcDataPostProcessor = self.postprocessor.rbcDataPostProcessor
        sCoords = self.sCoords(edge)
        mean = rbcDataPostProcessor.fieldAverageOnEdge(fieldName, sCoords,
                                                       edge, nAverage=nAverage)
        std = rbcDataPostProcessor.fieldStdOnEdge(fieldName, sCoords,
                                                  edge, nAverage=nAverage)
        self.plotSpatialProfile(mean, fieldName, edge, style=meanStyle)
        self.plotSpatialProfile(mean + std, fieldName, edge, style=stdStyle)
        self.plotSpatialProfile(mean - std, fieldName, edge, style=stdStyle)

    def plotSpatialProfile(self, y, ylabel_name, edge, **kwargs):
        """
        Plot a profile of a field average along a given edge.

        Args:
            fieldName (str): name of the field to plot
            y (np.ndarray): values to plot
            ylabel_name (str): string for y-axis (without unit)
            edge (int): internal edge index in graph
            **style (dict): line style for plotting
        """
        style = kwargs.get('style', {'linestyle': '-'})
        if len(y.shape) == 0 or y.shape[0] == 0:
            warnings.warn('No field value returned for edge {:d}, skipping...'.format(edge))
            return
        x = self.xPlot(edge)
        plt.plot(x, y, **style)
        labels.setXLabel('x', 'um')
        paramUnit = labels.openFOAMVarNameToUnit(ylabel_name)
        plottingParamName = labels.openFOAMVarNameToLatex(ylabel_name)
        labels.setYLabel(plottingParamName, paramUnit)

    def plotFieldTemporalProfile(self, fieldName, x, edgeIndex, nPath, **kwargs):
        """
        Plot a field value at x as a function of time, for the last nPath RBC paths.

        Args:
            fieldName (str): name of the field to plot
            x (float): x-coordinate at which the field is plotted
            edgeIndex (int): internal edge index in graph
            nPath (int): number of RBC paths to use
            **style (dict): line style for plotting
        """
        style = kwargs.get('style', {'linestyle': '-'})
        rbcDataPostProcessor = self.postprocessor.rbcDataPostProcessor
        s = x + self.simParams.sCoordOffset()
        t = rbcDataPostProcessor.timesOnEdgeCoordinate(x, edgeIndex, nPath)
        Hb = rbcDataPostProcessor.fieldOnEdge(fieldName, s, edgeIndex, nPath)
        plt.plot(t, Hb, **style)
        labels.setXLabel('t', 's')
        paramUnit = labels.openFOAMVarNameToUnit(fieldName)
        plottingParamName = labels.openFOAMVarNameToLatex(fieldName)
        labels.setYLabel(plottingParamName, paramUnit)

