"""
Plot the results of the model for diffusive interaction between capillaries
and the corresponding simulations.
"""

from collections import OrderedDict
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

from HbO2.plot import labels, styles
from HbO2.plot.parameterstudy import ParameterStudyPlotter
from HbO2.plot.simulated import GraphCasePlotter
from HbO2.plot.styles import set_COSH_rc_params
from HbO2.postprocess.diffusiveinteraction import DiffusiveInteractionPostProcessor
from plot.utils import include_zero_in_axis_range

set_COSH_rc_params()


style_scheme = styles.StyleSchemeCOSH()
mean_line_styles = {'krogh': style_scheme['meanNonlinear'],
                    'simple': style_scheme['meanExplicit'],
                    'linearized': style_scheme['meanExplicit'],
                    'equal_outfluxes': style_scheme['capNoModel'],
                    }
cap_line_styles = {'krogh': style_scheme['capNonlinear'],
                   'simple': style_scheme['capExplicit'],
                   'linearized': style_scheme['capExplicit'],
                   'equal_outfluxes': style_scheme['capNoModel'],
                   }
diff_line_styles = {'krogh': style_scheme['diffNonlinear'],
                    'simple': style_scheme['diffExplicit'],
                    'linearized': style_scheme['diffExplicit'],
                    'linearizedODE': style_scheme['diffLinear'],
                    'equal_outfluxes': style_scheme['diffNoModel'],
                    }


class DiffusiveInteractionModelPlotter(object):
    """
    Plot model computations of diffusive interaction between capillaries.

    Attributes:
        integrator: instance of COSHDiffusiveInteractionIntegrator
    """

    labels = {'krogh': 'nonlinear',
              'simple': 'explicit',
              'linearizedODE': 'linearized',
              'equal_outfluxes': 'equal fluxes',
              'expfit': 'exp. fit'}

    def __init__(self, integrator):
        self.integrator = integrator
        self.simParams = integrator.simParams
        self.xValues = np.linspace(0, self.simParams['domainLength'], 100)
        self._Hb = []

    @property
    def Hb(self):
        return self._Hb

    @Hb.getter
    def Hb(self):
        if self._Hb == []:
            self._Hb = self.integrator.saturationAtX(self.xValues)
        return self._Hb

    def xPlot(self):
        return 1e6*self.xValues

    def setHbLabels(self):
        labels.setXLabel('x', 'um')
        labels.setYLabel('S', '')

    def setHbDiffLabels(self):
        labels.setXLabel('x', 'um')
        labels.setYLabel('\Delta S', '')

    def plotAll(self):
        self.plotHb()
        self.plotHbMean()
        self.plotHbDifference()
        if self.simParams.cocurrentFlow():
            self.plotHbDifferenceLinearized()
            self.plotHbNoInteractionModel()

    def plotHb(self, **kwargs):
        """
        Plot a profile of the modeled hemoglobin saturation.

        Args:
            **style (dict): style for plotting
        """
        style = kwargs.get('style', cap_line_styles[self.integrator.model_name()])
        plt.plot(self.xPlot(), self.Hb, **style)
        self.setHbLabels()

    def plotHbMean(self, **kwargs):
        """
        Plot a profile of the mean modeled hemoglobin saturation.

        Args:
            **style (dict): style for plotting
        """
        style = kwargs.get('style', mean_line_styles[self.integrator.model_name()])
        plt.plot(self.xPlot(), np.mean(self.Hb, axis=1), 'b', **style)
        self.setHbLabels()

    def plotHbDifference(self, **kwargs):
        """
        Plot a profile of the modeled difference hemoglobin saturation

        Args:
            **style (dict): style for plotting
        """
        model_name = self.integrator.model_name()
        style = kwargs.get('style', diff_line_styles[model_name])
        if self.simParams.cocurrentFlow():
            hbDifference = np.abs(self.Hb[:, 1] - self.Hb[:, 0])
        else:
            hbDifference = np.abs(self.Hb[::-1, 1] - self.Hb[:, 0])
        plt.plot(self.xPlot(), hbDifference, **style)
        self.setHbDiffLabels()

    def plotHbDifferenceLinearized(self, **kwargs):
        """
        Plot a profile of the difference in hemoglobin saturation based on the linearized model

        Args:
            **style (dict): style for plotting
        """
        style = kwargs.get('style', diff_line_styles['linearizedODE'])
        hbDifference = self.integrator.linearizedSaturationDifference(self.xValues)
        plt.plot(self.xPlot(), hbDifference, **style)
        self.setHbDiffLabels()

    def plotHbNoInteractionModel(self, **kwargs):
        """
        Plot a profile of the hemoglobin saturation from the analytical model without
        diffusive interaction modeling.

        Args:
            **style (dict): style for plotting
        """
        style = kwargs.get('style', cap_line_styles['equal_outfluxes'])
        for Hb in self.integrator.inletHb:
            self.integrator.kroghSol.PO2Inlet = self.integrator.kroghSol.chem.hillP(Hb)
            S = self.integrator.kroghSol.saturationAtX(self.xValues)
            plt.plot(self.xPlot(), S, **style)
        self.setHbLabels()

    def plotOutfluxes(self, **kwargs):
        """
        Plot a profile of the modeled oxygen flux out of the capillary.

        Args:
            **style (dict): style for plotting
        """
        style = kwargs.get('style', cap_line_styles[self.integrator.model_name()])
        outfluxes = np.zeros(self.Hb.shape)
        for i in range(outfluxes.shape[0]):
            outfluxes[i,:] = self.integrator.interactionModel.outflux_per_length(self.Hb[i, :],
                                                                                 x=self.xValues[i])
        plt.plot(self.xPlot(), 1e6*outfluxes, **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('j_t', r'10^{-6} \mathrm{mlO_2\,m^{-1}\, s^{-1}}')

    def plotTissueRadii(self, **kwargs):
        """
        Plot a profile of the modeled tissue radii.

        Args:
            **style (dict): style for plotting
        """
        style = kwargs.get('style', cap_line_styles[self.integrator.model_name()])
        outfluxes = np.zeros(self.Hb.shape)
        tissue_radii = np.zeros(self.Hb.shape)
        for i in range(outfluxes.shape[0]):
            outfluxes[i,:] = self.integrator.interactionModel.outflux_per_length(self.Hb[i, :],
                                                                                 x=self.xValues[i])
            tissue_radii[i,:] = self.integrator.interactionModel.\
                _tissue_radius_from_outflux(outfluxes[i, :])
        plt.plot(self.xPlot(), 1e6*tissue_radii, **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('r_t', 'um')


class DiffusiveInteractionSimulationPlotter(GraphCasePlotter):
    """
    Plot simulations of diffusive interactions between capillaries
    """

    def __init__(self, postprocessor):
        super(DiffusiveInteractionSimulationPlotter, self).__init__(postprocessor,
                                                                    enforce_positive_flow=False)
        self.n_path = 50

    def plotHb(self, **kwargs):
        """
        Plot simulated profiles of hemoglobin saturation.

        Args:
            **style (dict): style for plotting
        """
        self.plotFieldAverageProfile('Hb_mean', 0, self.n_path, **kwargs)
        self.plotFieldAverageProfile('Hb_mean', 1, self.n_path, **kwargs)
        self.plotFieldAverageProfile('Hb_mean', 2, self.n_path, **kwargs)
        self.plotFieldAverageProfile('Hb_mean', 3, self.n_path, **kwargs)

    def plotHbMean(self, **kwargs):
        """
        Plot simulated profiles of mean hemoglobin saturation.

        Args:
            **style (dict): style for plotting
        """
        style = kwargs.get('style', {'linestyle': '--'})
        rbcDataPostProcessor = self.postprocessor.rbcDataPostProcessor
        hb0 = rbcDataPostProcessor.fieldAverageOnEdge('Hb_mean', self.sCoords(0), 0,
                                                      nAverage=self.n_path)
        hb1 = rbcDataPostProcessor.fieldAverageOnEdge('Hb_mean', self.sCoords(1), 1,
                                                      nAverage=self.n_path)
        if not self.simParams.cocurrentFlow():
            hb1 = hb1[::-1]
        plt.plot(self.xPlot(0), 0.5*(hb0 + hb1), 'b', **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('S', '')

    def plotHbDifference(self, **kwargs):
        """
        Plot difference between simulated profiles of hemoglobin saturation.

        Args:
            **style (dict): style for plotting
        """
        style = kwargs.get('style', {'linestyle': '--'})
        rbcDataPostProcessor = self.postprocessor.rbcDataPostProcessor
        hb0 = rbcDataPostProcessor.fieldAverageOnEdge('Hb_mean', self.sCoords(0), 0,
                                                      nAverage=self.n_path)
        hb1 = rbcDataPostProcessor.fieldAverageOnEdge('Hb_mean', self.sCoords(1), 1,
                                                      nAverage=self.n_path)
        if self.simParams.cocurrentFlow():
            hbDifference = hb0 - hb1
        else:
            hbDifference = hb0 - hb1[::-1]
        plt.plot(self.xPlot(0), hbDifference, **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel(r'\Delta S', '')

    def plotExponentialFitToHbDifference(self, **kwargs):
        """
        Plot the exponential fit to the simulated hemoglobin saturation difference.
        """
        style = kwargs.get('style', {'dashes': (1,2), 'color': 'r'})
        fittedHbDiff = self.postprocessor.exponentialFitToHbDifference()
        plt.plot(1e6*self.postprocessor.xValues, fittedHbDiff, **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('\Delta S', '')

    def plotDistalHb(self, **kwargs):
        """
        Plot hemoglobin saturation at the distal end of the domain as a function of time.

        Args:
            **style (dict): style for plotting
        """
        L = self.simParams.geometry()['domainLength']
        self.plotFieldTemporalProfile('Hb_mean', L, 0, self.n_path, **kwargs)
        self.plotFieldTemporalProfile('Hb_mean', L, 1, self.n_path, **kwargs)
        self.plotFieldTemporalProfile('Hb_mean', L, 2, self.n_path, **kwargs)
        self.plotFieldTemporalProfile('Hb_mean', L, 3, self.n_path, **kwargs)


class DiffusiveInteractionComparisonPlotter(object):
    """
    Plot simulations of diffusive interactions between capillaries
    """

    def __init__(self, postprocessor, model_names):
        """
        Constructor

        Args:
            postprocessor (DiffusiveInteractionPostProcessor):
            model_names (sequence of str): names of the integration model
        """
        self.modelPlotters = OrderedDict()  # use an OrderedDict to conserve plotting order
        for name in model_names:
            self.modelPlotters[name] = DiffusiveInteractionModelPlotter(postprocessor.integrator(name))
        self.simulPlotter = DiffusiveInteractionSimulationPlotter(postprocessor)

    def plotHbProfiles(self, plot_mean=True):
        self.plotHbSimulation(plot_mean)
        self.plotHbModel(plot_mean)

    def plotHbDiffProfiles(self, exponential_fit=True, linearized_ode=True):
        self.plotHbDiffSimulation()
        self.plotHbDiffModel(linearized_ode=linearized_ode)
        if exponential_fit:
            self.plotHbDiffExpFit()
        ax = plt.gca()
        ax.set_ylim([0.0, ax.get_ylim()[1]])

    def plotDistalHbTemporalProfile(self):
        self.simulPlotter.plotDistalHb()

    def plotHbModel(self, plot_mean):
        for model, plotter in self.modelPlotters.iteritems():
            plotter.plotHb(style=cap_line_styles[model])
            # only plot the mean for the Krogh model since the results are indistinguishable
            if model == 'krogh' and plot_mean:
                plotter.plotHbMean(style=mean_line_styles[model])

    def plotHbDiffModel(self, linearized_ode=True):
        for model, plotter in self.modelPlotters.iteritems():
            plotter.plotHbDifference(style=diff_line_styles[model])
        if linearized_ode:
            plotter = self.modelPlotters[self.modelPlotters.keys()[0]]
            plotter.plotHbDifferenceLinearized(style=diff_line_styles['linearizedODE'])

    def plotHbSimulation(self, plot_mean):
        self.simulPlotter.plotHb(style=style_scheme['capSim'])
        if plot_mean:
            self.simulPlotter.plotHbMean(style=style_scheme['meanSim'])
        modelPlotter = self.modelPlotters[self.modelPlotters.keys()[0]]
        x_min = min(modelPlotter.xPlot())
        x_max = max(modelPlotter.xPlot())
        plt.gca().set_xlim([x_min, x_max])

    def plotHbDiffSimulation(self):
        self.simulPlotter.plotHbDifference(style=style_scheme['diffSim'])
        modelPlotter = self.modelPlotters[self.modelPlotters.keys()[0]]
        x_min = min(modelPlotter.xPlot())
        x_max = max(modelPlotter.xPlot())
        plt.gca().set_xlim([x_min, x_max])

    def plotHbDiffExpFit(self):
        self.simulPlotter.plotExponentialFitToHbDifference(style=style_scheme['diffExpFit'])
        modelPlotter = self.modelPlotters[self.modelPlotters.keys()[0]]
        x_min = min(modelPlotter.xPlot())
        x_max = max(modelPlotter.xPlot())
        plt.gca().set_xlim([x_min, x_max])

    def plotHbDiffExpFit(self):
        expFitHbDiffStyle = {'style': {'linestyle': ':', 'color': 'y'}}
        self.simulPlotter.plotExponentialFitToHbDifference(**expFitHbDiffStyle)
        modelPlotter = self.modelPlotters[self.modelPlotters.keys()[0]]
        x_min = min(modelPlotter.xPlot())
        x_max = max(modelPlotter.xPlot())
        plt.gca().set_xlim([x_min, x_max])

class DiffusiveInteractionParameterStudyPlotter(ParameterStudyPlotter):
    """
    Plotter for a parameter study on diffusive interaction.

    The attribute post_processor should have type DiffusiveInteractionParameterStudyPostProcessor.
    """

    error_plots_styles = {'krogh': style_scheme['diffNonlinear'],
                          'simple': style_scheme['diffExplicit'],
                          'linearized': style_scheme['diffExplicit'],
                          'linearizedODE': style_scheme['diffLinear'],
                          'equal_outfluxes': style_scheme['diffNoModel']}

    def __init__(self, post_processor, fig_options, model_names=None, plot_mean=False,
                 plot_exponential_fit=False):
        super(DiffusiveInteractionParameterStudyPlotter, self).__init__(
            post_processor,
            fig_options
        )
        if model_names is None:
            self.model_names = ['krogh', 'simple']
        else:
            self.model_names = model_names
        self.plot_mean = plot_mean
        self.plot_exponential_fit = plot_exponential_fit
        self.param_name_to_locator = {'RBCVelocity': MultipleLocator(0.2)}

    def plot_all(self):
        self.plotHbProfiles()
        methods = ['plotAbsModelsErrorInFinalHbDifference',
                   'plotRelModelsErrorInFinalHbDifference',
                   'plotSaturationDifferenceDecayLength',
                   'plotComparisonSaturationDifferenceDecayLength',
                   'plotSaturationDifferenceDecayTime',
                   'plotComparisonSaturationDifferenceDecayTime']
        for method_name in methods:
            getattr(self, method_name)()
            self.fig_options.saveFig(method_name)
            plt.cla()

    def plotHbProfiles(self):
        for pp, caseName in zip(self.post_processor.post_processors,
                                self.post_processor.param_study.caseNames()):
            print "Plotting profiles for %s" % caseName
            plotter = DiffusiveInteractionComparisonPlotter(pp, self.model_names)
            plotter.plotHbProfiles(self.plot_mean)
            self.fig_options.saveFig('plotHbMeanProfiles_%s' % caseName)
            plt.clf()
            plotter.plotHbDiffProfiles(exponential_fit=self.plot_exponential_fit)
            self.fig_options.saveFig('plotHbDiffProfiles_%s' % caseName)
            plt.clf()

    def plotFinalHbDifferenceFromSimul(self, **kwargs):
        hb_diff = np.asarray([pp.finalHbDifferenceFromSimul()
                              for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(hb_diff, **kwargs)
        labels.setYLabel('\Delta S_v', '')

    def plotNormalizedFinalHbDifferenceFromSimul(self, **kwargs):
        hb_diff = 1e2*np.asarray([pp.finalHbDifferenceFromSimul() / pp.initialHbDifference()
                                  for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(hb_diff, **kwargs)
        labels.setYLabel('\Delta S_v/\Delta S_a', '\%')

    def plotAbsModelsErrorInFinalHbDifference(self):
        for model in self.model_names :
            if model != 'equal_outfluxes':
                self.plotAbsModelErrorInFinalHbDifference(model, style=self.error_plots_styles[model])
        self.plotAbsLinearizedODEErrorInFinalHbDifference(style=self.error_plots_styles['linearizedODE'])

    def plotRelModelsErrorInFinalHbDifference(self):
        for model in self.model_names:
            if model != 'equal_outfluxes':
                self.plotRelModelErrorInHbDifferenceDrop(model, style=self.error_plots_styles[model])
        self.plotRelLinearizedODEErrorInHbDifferenceDrop(style=self.error_plots_styles['linearizedODE'])

    def plotAbsModelErrorInFinalHbMean(self, model_name, **kwargs):
        abs_error = np.asarray([pp.absModelErrorInFinalHbMean(model_name)
                                for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(abs_error, **kwargs)
        labels.setYLabel('\mathrm{abs.\;error\;in\;mean}\;S_v', '')

    def plotRelModelErrorInFinalHbMean(self, model_name, **kwargs):
        rel_error = 1e2*np.asarray([pp.relModelErrorInFinalHbMean(model_name)
                                    for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(rel_error, **kwargs)
        labels.setYLabel('\mathrm{error\;in\;mean}\;S_v', '\%')

    def plotAbsModelErrorInFinalHbDifference(self, model_name, **kwargs):
        abs_error = np.asarray([pp.absModelErrorInFinalHbDifference(model_name)
                                for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(abs_error, **kwargs)
        labels.setYLabel('\mathrm{abs.\;error\;in}\;\Delta S_v', '')

    def plotRelModelErrorInHbDifferenceDrop(self, model_name, **kwargs):
        rel_error = 1e2*np.asarray([pp.relModelErrorInHbDifferenceDrop(model_name)
                                    for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(rel_error, **kwargs)
        labels.setYLabel('\mathrm{error\;in}\;\Delta S_a - \Delta S_v', '\%')

    def plotAbsLinearizedODEErrorInFinalHbDifference(self, **kwargs):
        abs_error = np.asarray([pp.absLinearizedODEErrorInFinalHbDifference()
                                for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(abs_error, **kwargs)
        labels.setYLabel('\mathrm{abs.\;error\;in}\;\Delta S_v', '')

    def plotRelLinearizedODEErrorInHbDifferenceDrop(self, **kwargs):
        rel_error = 1e2*np.asarray([pp.relLinearizedODEErrorInHbDifferenceDrop()
                                    for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(rel_error, **kwargs)
        labels.setYLabel('\mathrm{error\;in}\;\Delta S_a - \Delta S_v', '\%')

    def plotSaturationDifferenceDecayLength(self, **kwargs):
        decayLength = 1e6*np.asarray([pp.saturationDifferenceDecayLength()
                                      for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(decayLength, **kwargs)
        labels.setYLabel(r'L_{CI}', 'um')

    def plotTheoreticalSaturationDifferenceDecayLength(self, **kwargs):
        decayLength = 1e6*np.asarray([pp.theoreticalSaturationDifferenceDecayLength()
                                      for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(decayLength, **kwargs)
        labels.setYLabel(r'L_{CI}', 'um')

    def plotComparisonSaturationDifferenceDecayLength(self, **kwargs):
        self.plotSaturationDifferenceDecayLength(style=dict({'linestyle': ''}, **kwargs))
        self.plotTheoreticalSaturationDifferenceDecayLength(style=dict({'linestyle': '--'}, **kwargs))

    def plotSaturationDifferenceDecayTime(self, **kwargs):
        decayTime = 1e3*np.asarray([pp.saturationDifferenceDecayTime()
                                    for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(decayTime, **kwargs)
        labels.setYLabel(r'\tau_{CI}', 'ms')

    def plotTheoreticalSaturationDifferenceDecayTime(self, **kwargs):
        decayTime = 1e3*np.asarray([pp.theoreticalSaturationDifferenceDecayTime()
                                    for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(decayTime, **kwargs)
        labels.setYLabel(r'\tau_{CI}', 'ms')

    def plotComparisonSaturationDifferenceDecayTime(self, **kwargs):
        self.plotSaturationDifferenceDecayTime(style=dict({'linestyle': '-'}, **kwargs))
        self.plotTheoreticalSaturationDifferenceDecayTime(style=dict({'dashes': (3.6, 3.6)}, **kwargs))
