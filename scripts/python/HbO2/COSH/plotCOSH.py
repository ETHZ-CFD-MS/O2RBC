"""
Plot the results of the COSH ODE model and the corresponding simulations.
"""

import matplotlib.pyplot as plt
import numpy as np

from HbO2.COSH.postprocess import COSHPostProcessor
from HbO2.plot.plasmaoxygenation import MultipleLocator

from HbO2.plot import labels, styles
from HbO2.plot.parameterstudy import ParameterStudyPlotter
from HbO2.postprocess.factory.case import make_post_processor
from plot.utils import annotate_axis_corner, contourLevels, MultipleWithMaxNLocator, \
                       set_rounded_axis_limits, include_zero_in_axis_range

style_scheme = styles.StyleSchemeCOSH()


class COSHPlotter(object):
    """Plotting functions for a single OpenFOAM case for the study of COSH"""

    def __init__(self, postprocessor, **kwargs):
        self.postprocessor = postprocessor
        n_points = kwargs.get('nPoints', 101)
        self.xValues = np.linspace(0, self.postprocessor.simParams['domainLength'], n_points)
        self.postprocessor.integrator.KRI = self.postprocessor.KRI
        print "KRI = {:g}".format(self.postprocessor.KRI)

    def plotODEMeanHbProfile(self, **kwargs):
        """Plot the mean hemoglobin saturation profiles from the ODE model.
        """
        lineStyle = kwargs.get('style', style_scheme['meanRBCInter'])
        SMean = self.postprocessor.meanSaturationODE()
        plt.plot(1e6*self.xValues, SMean, **lineStyle)
        labels.setXLabel('x', 'um')
        labels.setYLabel('S', '')

    def plotODEMeanPlusMinusStdHbProfile(self, **kwargs):
        """Plot the mean +/- std of the hemoglobin saturation from the ODE model.
        """
        lineStyle = kwargs.get('style', style_scheme['meanPmStdRBCInter'])
        SMean = self.postprocessor.meanSaturationODE()
        SStd  = self.postprocessor.stdSaturationODE()
        plt.plot(1e6*self.xValues, SMean + SStd, **lineStyle)
        plt.plot(1e6*self.xValues, SMean - SStd, **lineStyle)
        labels.setXLabel('x', 'um')
        labels.setYLabel('S', '')

    def plotODEStdHbProfile(self, **kwargs):
        """Plots the standard deviation of the hemoglobin saturation from the
        ODE model.
        """
        style = kwargs.get('style', style_scheme['stdRBCInter'])
        SStd  = self.postprocessor.stdSaturationODE()
        plt.plot(1e6 * self.xValues, SStd, **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('\sigma_S', '')

    def plotODELinearizedStdHbProfile(self, **kwargs):
        """Plots the standard deviation of the hemoglobin saturation from the
        linearized ODE.
        """
        style = kwargs.get('style', style_scheme['stdLinear'])
        SLinStd  = self.postprocessor.linearizedStdSaturationODE()
        plt.plot(1e6 * self.xValues, SLinStd, **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('\sigma_S', '')

    def plotNoModelHbProfile(self, **kwargs):
        style = kwargs.get('style', style_scheme['meanPmStdNoModel'])
        S = self.postprocessor.saturationNoModel()
        plt.plot(1e6*self.xValues, S,  **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('S', '')

    def plotNoModelMeanPlusMinusStdHbProfile(self, **kwargs):
        """Plot the mean +/- std of the hemoglobin saturation with no RBC interaction modeling.
        """
        style = kwargs.get('style', style_scheme['meanPmStdNoModel'])
        SMean = self.postprocessor.meanSaturationNoModel()
        SStd  = self.postprocessor.stdSaturationNoModel()
        plt.plot(1e6*self.xValues, SMean + SStd, **style)
        plt.plot(1e6*self.xValues, SMean - SStd, **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('S', '')

    def plotNoModelStdHbProfile(self, **kwargs):
        style = kwargs.get('style', style_scheme['stdNoModel'])
        plt.plot(1e6*self.xValues,
                 self.postprocessor.initialStdSaturation()*np.ones(self.xValues.shape),
                 **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('\sigma_S', '')

    def plotSimMeanHbProfile(self, **kwargs):
        """Plot the mean simulated hemoblogin saturation"""
        style = kwargs.get('style', style_scheme['meanSim'])
        # S     = self.rbcDataPostProcessor.fieldAverageAlternating('Hb_mean',
        #                                     self.sValues(), 2, self.nAverage)
        SMean = self.postprocessor.rbcDataPostProcessor.fieldAverage(
            'Hb_mean',
            self.postprocessor.sValues(),
            self.postprocessor.nAverage
        )
        plt.plot(1e6*self.xValues, SMean, **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('S', '')

    def plotSimMeanPlusMinusHbProfile(self, **kwargs):
        """Plot the mean simulated hemoblogin saturation"""
        style = kwargs.get('style', style_scheme['meanPmStdSim'])
        SMean = self.postprocessor.rbcDataPostProcessor.fieldAverage(
            'Hb_mean',
            self.postprocessor.sValues(),
            self.postprocessor.nAverage
        )
        SStd = self.postprocessor.rbcDataPostProcessor.fieldStd(
            'Hb_mean',
            self.postprocessor.sValues(),
            self.postprocessor.nAverage
        )
        plt.plot(1e6*self.xValues, SMean + SStd, **style)
        plt.plot(1e6*self.xValues, SMean - SStd, **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('S', '')

    def plotSimHbProfile(self, nProfiles, **kwargs):
        """Plot individual simulated hemoblogin saturation profiles."""
        style = kwargs.get('style', style_scheme['individualRBC'])
        S = self.postprocessor.rbcDataPostProcessor.fieldOnLastPaths(
            'Hb_mean',
            self.postprocessor.sValues(),
            nProfiles
        )
        plt.plot(1e6*self.xValues, S, **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('S', '')

    def plotSimStdHbProfile(self, **kwargs):
        """Plot the standard deviation of the simulated hemoblogin saturation"""
        style = kwargs.get('style', style_scheme['stdSim'])
        SStd = self.postprocessor.rbcDataPostProcessor.fieldStd(
            'Hb_mean',
            self.postprocessor.sValues(),
            self.postprocessor.nAverage
        )
        plt.plot(1e6 * self.xValues, SStd, **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('\sigma_S', '')

    def plotExponentialFitStdHb(self, **kwargs):
        """
        Plot the exponential fit to the standard deviation of the simulated hemoglobin
        saturation.
        """
        style = kwargs.get('style', style_scheme['stdExpFit'])
        fittedStd = self.postprocessor.exponentialFitToStdHb()
        plt.plot(1e6 * self.xValues, fittedStd, **style)
        labels.setXLabel('x', 'um')
        labels.setYLabel('\sigma_S', '')


class COSHParamStudyPlotter(ParameterStudyPlotter):
    """Plot results of a parameters study for the study of COSH"""

    def __init__(self, post_processor, fig_options, nAverage=20):
        super(COSHParamStudyPlotter, self).__init__(post_processor, fig_options)
        self.nAverage = nAverage

    def plot_all(self):
        self.plotKRI()  # plot this first as it also optimizes KRI
        self.plotKOS()
        # self.plotProfiles()
        methods = ['plotRelErrorInStdDrop',
                   'plotOscillationRadius',
                   'plotRelativeOscillationRadius',
                   'plotIntegralOscillationRadius',
                   'plotIntegralOscillationDistance',
                   'plotStdSaturationDecayLength',
                   'plotStdSaturationDecayTime',
                   'plotKOSVersusIntegralOscillationDistance']
        for method_name in methods:
            getattr(self, method_name)()
            self.fig_options.saveFig(method_name)
            plt.cla()

    def plotKRI(self):
        KRI = np.asarray([pp.KRI for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(1e-6 * KRI)
        labels.setYLabel('K_{RI}', 'IVR')
        self.fig_options.saveFig('plotKRI')
        plt.clf()

    def plotKOS(self):
        KOS = np.asarray([pp.KOS for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(1e-6 * KOS)
        labels.setYLabel('K_{OS}', 'IVR')
        self.fig_options.saveFig('plotKOS')
        plt.clf()

    def plotOscillationRadius(self):
        oscillRads = 1e6*np.asarray([pp.oscillationRadius()
                                     for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(oscillRads)
        labels.setYLabel('\mathrm{oscill.\;radius}', 'um')

    def plotRelativeOscillationRadius(self):
        oscillRads = 1e6*np.asarray([pp.relativeOscillationRadius()
                                     for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(oscillRads)
        labels.setYLabel('\mathrm{rel.\;oscill.\;radius}', 'um')

    def plotIntegralOscillationRadius(self):
        oscillRads = 1e6*np.asarray([pp.integralOscillationRadius()
                                     for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(oscillRads)
        labels.setYLabel('\mathrm{int.\;oscill.\;radius}', 'um')

    def plotIntegralOscillationDistance(self):
        oscillDistance = 1e6*np.asarray([pp.integralOscillationPenetrationDistance()
                                         for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(oscillDistance)
        labels.setYLabel('\mathrm{int.\;oscill.\;length}', 'um')

    def plotRelErrorInStdDrop(self):
        self.plotRelODEModelErrorInStdDrop(style=style_scheme['stdRBCInter'])
        self.plotRelLinearizedModelErrorInStdDrop(style=style_scheme['stdLinear'])
        self.plotRelExponentialFitErrorInStdDrop(style=style_scheme['stdExpFit'])
        ax = plt.gca()
        ax.yaxis.set_major_locator(MultipleLocator(1))
        set_rounded_axis_limits(ax, 1.0, 'y')
        include_zero_in_axis_range(ax, 'y')
        labels.setYLabel('\mathrm{error\;in}\;\sigma_{S,a} - \sigma_{S,v}', '\%')

    def plotRelODEModelErrorInStdDrop(self, **kwargs):
        relError = 1e2*np.asarray([pp.relODEModelErrorInStdDrop()
                                   for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(relError, **kwargs)

    def plotRelLinearizedModelErrorInStdDrop(self, **kwargs):
        relError = 1e2*np.asarray([pp.relLinearizedModelErrorInStdDrop()
                                   for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(relError, **kwargs)

    def plotRelExponentialFitErrorInStdDrop(self, **kwargs):
        relError = 1e2*np.asarray([pp.relExponentialFitErrorInStdDrop()
                                   for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(relError, **kwargs)

    def plotStdSaturationDecayLength(self, **kwargs):
        decayLength = 1e6*np.asarray([pp.stdSaturationDecayLength()
                                      for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(decayLength, **kwargs)
        labels.setYLabel(r'L_{RI}', 'um')

    def plotStdSaturationDecayTime(self, **kwargs):
        decayTime = 1e3*np.asarray([pp.stdSaturationDecayTime()
                                    for pp in self.post_processor.post_processors])
        self.plot_variable_single_parameter(decayTime, **kwargs)
        labels.setYLabel(r'\tau_{RI}', 'ms')

    def plotProfiles(self):
        for caseName in self.post_processor.param_study.caseNames():
            print "Plotting profiles for %s" % caseName
            postprocessor = make_post_processor(caseName)
            plotter = COSHPlotter(postprocessor)
            plotter.plotODEMeanHbProfile()
            plotter.plotSimMeanHbProfile()
            self.fig_options.saveFig('plotCOSHHbProfiles_%s' % caseName)
            plt.clf()

            plotter.plotODEStdHbProfile()
            plotter.plotODELinearizedStdHbProfile()
            plotter.plotSimStdHbProfile()
            plotter.plotExponentialFitStdHb()
            self.fig_options.saveFig('plotCOSHStdProfiles_%s' % caseName)
            plt.gca().set_yscale('log')
            self.fig_options.saveFig('plotCOSHStdProfilesLogScale_%s' % caseName)
            plt.clf()

    def plotKOSVersusIntegralOscillationDistance(self):
        """Plot KOS as a function of the integral oscillation length and a linear fit."""
        oscillLengths = 1e6*np.asarray([pp.integralOscillationPenetrationDistance()
                                        for pp in self.post_processor.post_processors])
        KOS = 1e-6*np.asarray([pp.KOS for pp in self.post_processor.post_processors])
        plt.plot(oscillLengths, KOS, '.')

        slope, intercept, rvalue, pvalue, stderr = \
            self.post_processor.KOSFitter.fitKOS(oscillLengths, KOS)
        rfit = np.linspace(np.min(oscillLengths), np.max(oscillLengths), 100)
        plt.plot(rfit, slope*rfit + intercept, 'k-')
        labels.setXLabel(r'\Delta r_{osc}', 'um')
        labels.setYLabel(r'K_{OS}', 'IVR')


class LDvRBCCOSHParamStudyPlotter(COSHParamStudyPlotter):
    """Plot results of a LD_U parameter study for the study of COSH"""

    def __init__(self, paramStudyPostProcessor, fig_options):
        super(LDvRBCCOSHParamStudyPlotter, self).__init__(paramStudyPostProcessor,
                                                          fig_options)
        self.cmap = plt.cm.jet
        # self.cmap = plt.cm.hot

    def plot_all(self):
        # self.plotProfiles()
        self.plotKRI()
        self.plotKOS()
        self.plotKOS_LD()
        self.plotKOS_U()
        self.plotMultiPanelKOSContourAndOscillationDistance()
        self.plotOscillationRadii()
        self.plotRelativeOscillationRadii()
        self.plotIntegralOscillationRadii()
        self.plotIntegralOscillationLengths()
        self.plotStdSaturationDecayLength()
        self.plotStdSaturationDecayTime()
        self.plotNormalizedTissueOscillationIntegral()
        self.plotKOSVersusIntegralOscillationDistance()

    def plotKRI(self):
        U, LD = self.post_processor.param_study.meshGrid()
        Z = 1e-6*self.reshapeListToArray(
            [pp.KRI for pp in self.post_processor.post_processors])
        # levels = contourLevels(Z, 0.5, 1)
        # plt.contourf(LD, 1e3*U, Z, levels)
        # plt.colorbar()
        levels = contourLevels(Z, 1, 1)
        cplot = plt.contour(LD, 1e3*U, Z, levels, cmap=self.cmap)
        plt.clabel(cplot, inline=1, fontsize=12, fmt='%g')
        self.setAxisLabels()
        self.fig_options.saveFig('plotKRI')
        plt.clf()

    def plotKOS(self):
        U, LD = self.post_processor.param_study.meshGrid()
        Z = 1e-6*self.reshapeListToArray(
            [pp.KOS for pp in self.post_processor.post_processors])
        # levels = contourLevels(Z, 0.25, 0.5)
        # plt.contourf(LD, 1e3*U, Z, levels)
        # plt.colorbar()
        levels = contourLevels(Z, 0.5, 0.5)
        cplot = plt.contour(LD, 1e3*U, Z, levels, cmap=self.cmap)
        plt.clabel(cplot, inline=1, fontsize=12, fmt='%g')
        self.setAxisLabels()
        self.fig_options.saveFig('plotKOS')
        plt.clf()

    def plotKOS_LD(self):
        LDValues = self.post_processor.param_study['LD']
        UValues  = 1e3*self.post_processor.param_study['U']
        Z = 1e-6*self.reshapeListToArray(
            [pp.KOS for pp in self.post_processor.post_processors])
        for i, U in enumerate(UValues):
            plt.plot(LDValues, Z[:,i], label='U=%g'%U, color='k')
        plt.ylim([0, np.max(Z)])
        ax = plt.gca()
        ax.annotate("U",
                    xy=(0.5*(LDValues[0] + LDValues[-1]), 0.7*np.max(Z)),
                    xycoords='data',
                    xytext=(0.35*(LDValues[0] + LDValues[-1]), 0.05*np.max(Z)),
                    textcoords='data',
                    arrowprops=dict(arrowstyle='<-', color='b',
                        connectionstyle="arc3"),
            )
        labels.setXLabel('LD', '')
        labels.setYLabel('K_{OS}', 'IVR')
        self.fig_options.saveFig('plotKOS_LD')
        plt.clf()

    def plotKOS_U(self):
        LDValues = self.post_processor.param_study['LD']
        UValues  = 1e3*self.post_processor.param_study['U']
        Z = 1e-6*self.reshapeListToArray(
            [pp.KOS for pp in self.post_processor.post_processors])
        for i, LD in enumerate(LDValues):
            plt.plot(UValues, Z[i,:], color='k', label='LD=%g' % LD)
        plt.ylim([0, np.max(Z)])
        ax = plt.gca()
        ax.annotate("LD",
                    xy=(0.5*(UValues[0] + UValues[-1]), 0.7*np.max(Z)),
                    xycoords='data',
                    xytext=(0.3*(UValues[0] + UValues[-1]), 0.05*np.max(Z)),
                    textcoords='data',
                    arrowprops=dict(arrowstyle='<-', color='b',
                        connectionstyle="arc3"),
            )
        labels.setXLabel('U', 'mm/s')
        labels.setYLabel('K_{OS}', 'IVR')
        self.fig_options.saveFig('plotKOS_U')
        plt.clf()

    def plotMultiPanelKOSContourAndOscillationDistance(self, manual_labels=False):
        fig, axs = plt.subplots(1, 2)
        U, LD = self.post_processor.param_study.meshGrid()
        Z = 1e-6*self.reshapeListToArray(
            [pp.KOS for pp in self.post_processor.post_processors])
        levels = contourLevels(Z, 0.5, 0.5)
        cplot = axs[0].contour(LD, 1e3*U, Z, levels, cmap=self.cmap)
        # cplot = axs[0].contour(LD, 1e3*U, Z, levels, colors='k')
        if manual_labels:
            plt.clabel(cplot, inline=True, inline_spacing=24,
                       fontsize=10, fmt='%g', manual=True)
        else:
            plt.clabel(cplot, inline=1, fontsize=12, fmt='%g')
        self.setAxisLabels(ax=axs[0])
        annotate_axis_corner(axs[0], 'A')

        oscillLengths = 1e6*np.asarray([pp.integralOscillationPenetrationDistance()
                                        for pp in self.post_processor.post_processors])
        KOS = 1e-6*np.asarray([pp.KOS for pp in self.post_processor.post_processors])
        axs[1].plot(oscillLengths, KOS, '.', color='b', markersize=3)

        slope, intercept, rvalue, pvalue, stderr = \
            self.post_processor.KOSFitter.fitKOS(oscillLengths, KOS)
        rfit = np.linspace(np.min(oscillLengths), np.max(oscillLengths), 100)
        axs[1].xaxis.set_major_locator(MultipleWithMaxNLocator(0.5, 10))
        axs[1].plot(rfit, slope*rfit + intercept, 'k-')
        labels.setXLabel(r'\Delta r_{\mathrm{osc}}', 'um', ax=axs[1])
        labels.setYLabel(r'K_{OS}', 'IVR', ax=axs[1])
        set_rounded_axis_limits(axs[1], base=0.5, axis='x')
        set_rounded_axis_limits(axs[1], base=1.0, axis='y')
        annotate_axis_corner(axs[1], 'B')
        self.fig_options.saveFig('plotKOSContourAndOscillationDistance')
        plt.clf()


    def plotOscillationRadii(self):
        amplitude = COSHPostProcessor.oscillationAmplitude
        U, LD = self.post_processor.param_study.meshGrid()
        Z = 1e6*self.reshapeListToArray(
            [pp.oscillationRadius() for pp in self.post_processor.post_processors])
        levels = contourLevels(Z, 0.5, 1)
        plt.contourf(LD, 1e3*U, Z, levels)
        plt.colorbar()
        self.setAxisLabels()
        self.fig_options.saveFig('plotOscillationRadii_{:g}mmHg'.format(amplitude))
        plt.clf()

    def plotRelativeOscillationRadii(self):
        amplitude = COSHPostProcessor.relOscillationAmplitude
        U, LD = self.post_processor.param_study.meshGrid()
        Z = 1e6*self.reshapeListToArray(
            [pp.relativeOscillationRadius()
             for pp in self.post_processor.post_processors])
        levels = contourLevels(Z, 0.5, 1)
        plt.contourf(LD, 1e3*U, Z, levels)
        plt.colorbar()
        self.setAxisLabels()
        self.fig_options.saveFig('plotRelOscillationRadii_{:g}'.format(amplitude))
        plt.clf()

    def plotIntegralOscillationRadii(self):
        U, LD = self.post_processor.param_study.meshGrid()
        Z = 1e6*self.reshapeListToArray(
            [pp.integralOscillationRadius() for pp in self.post_processor.post_processors])
        levels = contourLevels(Z, 0.1, 0.5)
        plt.contourf(LD, 1e3*U, Z, levels)
        plt.colorbar()
        self.setAxisLabels()
        self.fig_options.saveFig('plotIntegralOscillationRadii')
        plt.clf()

    def plotIntegralOscillationLengths(self):
        U, LD = self.post_processor.param_study.meshGrid()
        Z = 1e6*self.reshapeListToArray(
            [pp.integralOscillationPenetrationDistance()
             for pp in self.post_processor.post_processors])
        levels = contourLevels(Z, 0.1, 0.5)
        plt.contourf(LD, 1e3*U, Z, levels)
        plt.colorbar()
        self.setAxisLabels()
        self.fig_options.saveFig('plotIntegralOscillationLengths')
        plt.clf()

    def plotStdSaturationDecayLength(self):
        U, LD = self.post_processor.param_study.meshGrid()
        Z = 1e6*self.reshapeListToArray([pp.stdSaturationDecayLength()
                                         for pp in self.post_processor.post_processors])
        levels = contourLevels(Z, 10, 50)
        plt.contourf(LD, 1e3*U, Z, levels)
        plt.colorbar()
        self.setAxisLabels()
        self.fig_options.saveFig('plotStdSaturationDecayLength')
        plt.clf()

    def plotStdSaturationDecayTime(self):
        U, LD = self.post_processor.param_study.meshGrid()
        Z = 1e3*self.reshapeListToArray([pp.stdSaturationDecayTime()
                                         for pp in self.post_processor.post_processors])
        levels = contourLevels(Z, 10, 50)
        plt.contourf(LD, 1e3*U, Z, levels)
        plt.colorbar()
        self.setAxisLabels()
        self.fig_options.saveFig('plotStdSaturationDecayTime')
        plt.clf()

    def plotStdSaturationDecayLength(self):
        U, LD = self.postProcessor.paramStudy.meshGrid()
        Z = 1e6*self.reshapeListToArray([pp.stdSaturationDecayLength()
                                         for pp in self.postProcessor.postProcessors])
        levels = contourLevels(Z, 10, 50)
        plt.contourf(LD, 1e3*U, Z, levels)
        plt.colorbar()
        self.setAxisLabels()
        self.figOptions.saveFig('plotStdSaturationDecayLength')
        plt.clf()

    def plotStdSaturationDecayTime(self):
        U, LD = self.postProcessor.paramStudy.meshGrid()
        Z = 1e3*self.reshapeListToArray([pp.stdSaturationDecayTime()
                                         for pp in self.postProcessor.postProcessors])
        levels = contourLevels(Z, 10, 50)
        plt.contourf(LD, 1e3*U, Z, levels)
        plt.colorbar()
        self.setAxisLabels()
        self.figOptions.saveFig('plotStdSaturationDecayTime')
        plt.clf()

    def plotNormalizedTissueOscillationIntegral(self):
        U, LD = self.post_processor.param_study.meshGrid()
        Z = self.reshapeListToArray(
            [pp.normalizedTissueOscillationIntegral()
             for pp in self.post_processor.post_processors])
        levels = contourLevels(Z, 0.005, 0.01)
        plt.contourf(LD, 1e3*U, Z, levels)
        plt.colorbar()
        self.setAxisLabels()
        self.fig_options.saveFig('plotNormalizedTissueOscillationIntegral')
        plt.clf()

    def reshapeListToArray(self, l):
        """Reshape a list with values corresponding to the study parameters to 
        an array for plotting"""
        nLD = len(self.post_processor.param_study['LD'])
        nU  = len(self.post_processor.param_study['U'])
        return np.reshape(l, (nLD, nU))

    def setAxisLabels(self, **kwargs):
        labels.setXLabel('LD', '', **kwargs)
        labels.setYLabel('vrbc', 'mm/s', **kwargs)
