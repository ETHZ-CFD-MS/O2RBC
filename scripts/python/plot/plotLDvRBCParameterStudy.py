#!/usr/bin/env python
"""
Plot the results of a LD_U parameter study.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from HbO2.model.kroghSolution import KroghSolution2DCone
from HbO2.plot.labels import setXLabel, setYLabel
from HbO2.postprocess.factory.parameterstudy import make_param_study_post_processor
from HbO2.postprocess.parameterstudy import LDvRBCParameterStudyPostProcessor
from HbO2.setup.parameterStudy import LDvRBCParameterStudy
from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric
from plot.figureoptions import FigureOptions


class plotLDvRBCParameterStudy(object):
    flowLevels = np.arange(20, 241, 20)
    UFactor = 1./1e-3
    hypoxiaThreshold = 2.0 # [mmHg]

    def __init__(self, postProcessor):
        """
        Initialize an instance: constructs the attribute kroghSol and choose the radius
        for plotting PO2 as a function of the geometry.

        Args:
            postProcessor (LDvRBCParameterStudyPostProcessor): postprocessor for plotting
        """
        self.postProcessor = postProcessor
        simParams = IOHbO2ParametersAxisymmetric(self.postProcessor.param_study['baseCasePath'])
        self.kroghSol = KroghSolution2DCone(simParams)
        self.kroghSol.convO2Transport = self.postProcessor.kroghSols[0].convO2Transport
        if simParams.geometry().isGlomerulus():
            self.R = 12.6e-6
            self.PO2Levels = np.arange(0., 80., 8)
            self.PO2Levels[0] = 2.0
        else:
            self.R = 17.6e-6
            self.PO2Levels = np.arange(0., 80., 8)
            self.PO2Levels[0] = 2.0

        self.cmap = plt.cm.jet

    def contourPlot(self):
        data = np.loadtxt(self.postProcessor.PO2MeanOutputFilePath(), delimiter="\t", skiprows=1)
        n_LD = len(set(data[:, 0]))
        n_U  = len(set(data[:, 1]))
        LD = np.reshape(data[:, 0], (n_LD, n_U))
        U  = np.reshape(data[:, 1], (n_LD, n_U))
        Z  = np.reshape(data[:, self.postProcessor.probeIdx + 2], (n_LD, n_U))

        plt.contourf(LD, self.UFactor*U, Z, [0, self.hypoxiaThreshold], colors='0.9')
        C1 = plt.contour(LD, self.UFactor*U, Z, self.PO2Levels, cmap=self.cmap)
        plt.clabel(C1, inline=1, fontsize=12, fmt='%g')
        plt.xlabel(r'$\mathrm{LD} \; [-]$')
        plt.ylabel(r'$v_{\mathrm{rbc}} \; [\mathrm{mm\,s^{-1}}]$')

        L_RBC = self.kroghSol.RBCLength()
        flow = U*LD/L_RBC
        C2 = plt.contour(LD, self.UFactor*U, flow, self.flowLevels,
                         colors='k', linewidths=0.5, linestyles='dotted')
        plt.clabel(C2, inline=1, fontsize=7, fmt='%g')

    def contourPlotAnalyticalSolutionSecomb(self, **kwargs):
        K0 = kwargs.get('K0', self.kroghSol.intravascularResistanceLDHalf)
        xPosition = self.postProcessor.probePositionFromIndex(self.postProcessor.probeIdx)
        self.kroghSol.intravascularResistanceLDHalf = K0
        paramStudy = self.postProcessor.param_study
        LDValues = [v[0] for v in paramStudy.paramValues()]
        UValues  = [v[1] for v in paramStudy.paramValues()]
        LDAnalytical = np.linspace(min(LDValues), max(LDValues), 50)
        UAnalytical  = np.linspace(min(UValues), max(UValues), 50)

        tissuePO2 = self.kroghSol.PO2TissueArray(xPosition, self.R, LDAnalytical, UAnalytical)

        U, LD = np.meshgrid(UAnalytical, LDAnalytical)
        C2 = plt.contour(LD, self.UFactor*U, tissuePO2, self.PO2Levels, 
                linestyles='dashed', cmap=self.cmap)
        # plt.clabel(C2, inline=1, fontsize=8, fmt='%g')
        return self.kroghSol.intravascularResistanceLDHalf

    def contourPlotMTC(self):
        MTCFactor = 1e2  # convert from mlO2 m^-2 s^-1 to 1e-6 mlO2 cm^-2 s^-1
        data = np.loadtxt(self.postProcessor.MTCOutputFilePath(),
                          delimiter="\t", skiprows=1)
        n_LD = len(set(data[:,0]))
        n_U  = len(set(data[:,1]))
        LD = np.reshape(data[:,0], (n_LD, n_U))
        U  = np.reshape(data[:,1], (n_LD, n_U))
        MTC = np.reshape(data[:,2], (n_LD, n_U))
        Nu  = np.reshape(data[:,3], (n_LD, n_U))

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twiny()
        MTCContours = np.arange(0.4, MTCFactor*max(MTC.flatten()), 0.4)
        NuContours = np.arange(0.2, max(Nu.flatten()), 0.2)
        C1 = plt.contour(LD, self.UFactor*U, MTCFactor*MTC, MTCContours,
                         linestyles='--', colors='k')
        C2 = plt.contour(LD, self.UFactor*U, Nu, NuContours, colors='k')
        fig.canvas.draw()
        ax1.set_xlabel(r'$\mathrm{LD} \; [-]$')
        ax1.set_ylabel(r'$v_{\mathrm{rbc}} \; [\mathrm{mm\,s^{-1}}]$')
        # ax1.set_xlim([0.1, 0.6])
        # ax1.set_ylim([0.4, 2.4])
        xmin = min(LD.flatten())
        xmax = max(LD.flatten())
        ymin = self.UFactor*min(U.flatten())
        ymax = self.UFactor*max(U.flatten())
        ax1.set_xlim([xmin, xmax])
        ax1.set_ylim([ymin, ymax])

        simParams = IOHbO2ParametersAxisymmetric(self.postProcessor.param_study['baseCasePath'])
        radiusPlasma = simParams.geometry()['radiusPlasma']
        radiusRBC = simParams['radiusRBC']
        def inverse_tick_function(Ht):
            return (radiusPlasma/radiusRBC)**2*Ht
        Ht_tick_values = np.linspace(0.1, 0.5, 9)
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(inverse_tick_function(Ht_tick_values))
        ax2.set_xticklabels(['%.2f' % x for x in Ht_tick_values])
        ax2.set_xlabel(r'$H_t \; [-]$')

        fmt={}  # put empty labels to C1 to create space for C2's labels
        for l in C1.levels:
            fmt[l] = '  '
        # fmt='%g'
        plt.clabel(C1, inline=True, inline_spacing=48,
                   fontsize=8, fmt=fmt, manual=True)
        plt.clabel(C2, inline=True, inline_spacing=44,
                   fontsize=10, fmt='%g', manual=True)

    def linePlotMTCNu(self, plotField='Nu'):
        MTCFactor = 1e2  # convert from mlO2 m^-2 s^-1 to 1e-6 mlO2 cm^-2 s^-1
        UPlot = [0.4e-3, 2.4e-3]
        # LDPlot = [0.1, 0.3, 0.5, 0.7, 0.9]
        LDPlot = np.linspace(0.1, 0.9, 9)

        data = np.loadtxt(self.postProcessor.MTCOutputFilePath(),
                          delimiter="\t", skiprows=1)
        n_LD = len(set(data[:,0]))
        n_U  = len(set(data[:,1]))
        LD = np.reshape(data[:,0], (n_LD, n_U))
        U  = np.reshape(data[:,1], (n_LD, n_U))
        MTC = np.reshape(data[:,2], (n_LD, n_U))
        Nu  = np.reshape(data[:,3], (n_LD, n_U))
        if plotField == 'Nu':
            field = Nu
        elif plotField == 'MTC':
            field = MTCFactor*MTC
        LDvec = LD[:,0]
        Uvec = U[0,:]
        idxLDPlot = [np.argmin(np.abs(LDvec - l)) for l in LDPlot]
        idxUPlot = [np.argmin(np.abs(Uvec - v)) for v in UPlot]

        styles = [{'linestyle':'-'},
                  {'linestyle':'--'},
                  {'linestyle':'-.', 'dashes': (4, 2, 1, 2)},
                  {'linestyle': '-', 'dashes': (1, 2)}]

        fig, axs = plt.subplots(3, 1)
        for i, style in zip(idxUPlot, styles):
            axs[0].plot(LD, field[:, i], **style)
        axs[0].set_xlim(np.min(LD), np.max(LD))
        setXLabel('LD', '-', ax=axs[0])
        if plotField == 'Nu':
            setYLabel(r'\mathrm{Nu}', '-', ax=axs[0])
        elif plotField == 'MTC':
            setYLabel(r'k', '', ax=axs[0])

        axtop = axs[0].twiny()
        simParams = IOHbO2ParametersAxisymmetric(self.postProcessor.param_study['baseCasePath'])
        radiusPlasma = simParams.geometry()['radiusPlasma']
        radiusRBC = simParams['radiusRBC']
        def inverse_tick_function(Ht):
            return (radiusPlasma/radiusRBC)**2*Ht
        Ht_tick_values = np.linspace(0.1, 0.5, 9)
        axtop.set_xlim(axs[0].get_xlim())
        axtop.set_xticks(inverse_tick_function(Ht_tick_values))
        axtop.set_xticklabels(['%.2f' % x for x in Ht_tick_values])
        axtop.set_xlabel(r'$H_t \; [-]$')

        field_U_ratio = field[:,-1]/field[:,0]
        field_LD_ratio = field/np.transpose(np.tile(field[:,0], (n_U, 1)))
        axs[1].plot(LDvec, field_U_ratio, **styles[0])
        axs[1].set_xlim(np.min(LD), np.max(LD))
        setXLabel('LD', '-', ax=axs[1])
        if plotField == 'Nu':
            setYLabel(r'\mathrm{Nu}_{v_{\mathrm{rbc}}=2.4}/\mathrm{Nu}_{v_{\mathrm{rbc}}=0.2}', '', ax=axs[1])
        elif plotField == 'MTC':
            setYLabel(r'k_{v_{\mathrm{rbc}}=2.4}/k_{v_{\mathrm{rbc}}=0.2}', '', ax=axs[1])
    
        for i in idxLDPlot:
            axs[2].plot(1e3 * Uvec, np.transpose(field_LD_ratio[i, :]),
                        color='k', linewidth=0.5)
        axs[2].set_xlim(np.min(1e3*U), np.max(1e3*U))
        axs[2].set_ylim(1., axs[2].get_ylim()[1])
        axs[2].set_xticks(1e3*Uvec[1:])
        setXLabel('vrbc', 'mm/s', ax=axs[2])

        axs[2].annotate("LD",
                        xy = (2.1, 1.00), 
                        xycoords='data',
                        xytext = (1.78, 1.26),
                        textcoords='data',
                        arrowprops=dict(arrowstyle='->', color='k',
                            connectionstyle="arc3"),
            )

        if plotField == 'Nu':
            setYLabel(r'\mathrm{Nu}/\mathrm{Nu}_{v_{\mathrm{rbc}}=0.2}', '', ax=axs[2])
            axs[0].annotate(r'$\mathbf{A}$', xy=(0.93, 0.78), xycoords='axes fraction')
            axs[1].annotate(r'$\mathbf{B}$', xy=(0.93, 0.89), xycoords='axes fraction')
            axs[2].annotate(r'$\mathbf{C}$', xy=(0.93, 0.89), xycoords='axes fraction')
        elif plotField == 'MTC':
            setYLabel(r'k/k_{v_{\mathrm{rbc}}=0.2}', '', ax=axs[2])
            axs[0].annotate(r'$\mathbf{A}$', xy=(0.91, 0.89), xycoords='axes fraction')
            axs[1].annotate(r'$\mathbf{B}$', xy=(0.91, 0.89), xycoords='axes fraction')
            axs[2].annotate(r'$\mathbf{C}$', xy=(0.91, 0.89), xycoords='axes fraction')

    def quiverPlotPO2Gradient(self, x):
        LDValues = self.postProcessor.param_study.LDValues
        UValues  = self.postProcessor.param_study.UValues

        ddLD, ddvRBC = self.kroghSol.gradPO2(x, LDValues, UValues)
        ddLD   = np.transpose(ddLD)
        ddvRBC = np.transpose(ddvRBC)

        LD, U = np.meshgrid(LDValues, UValues)
        L_RBC = self.kroghSol.RBCLength
        flow = U*LD/L_RBC
        M = np.zeros(LD.shape, dtype='bool')
        M[flow < 20] = True
        M[LD < 0.15] = True
        
        ddLD   = np.ma.masked_array(ddLD,   mask=M)
        ddvRBC = np.ma.masked_array(ddvRBC, mask=M)

        Q = plt.quiver(LD, self.UFactor*U, ddLD, 1/self.UFactor*ddvRBC)
        plt.axis('equal')

    ## Plot the convective part and the IV resistance part of PO2 at the wall as a function
    #  oxygen consumption. Also plot the derivatives.
    def convectiveAndIVResistancePartsWRTM(self, x, LD, vRBC):
        MValues = np.linspace(0, 3000, 21)
        convectiveDrop   = np.zeros(MValues.shape)
        IVResistanceDrop = np.zeros(MValues.shape)
        ddLD             = np.zeros(MValues.shape)
        ddvRBC           = np.zeros(MValues.shape)
        for i, M in enumerate(MValues):
            self.kroghSol.M = M
            convectiveDrop[i]   = self.kroghSol.convectivePO2Drop(x)
            IVResistanceDrop[i] = self.kroghSol.intravascResistancePO2Drop(x)
            ddLD[i]             = self.kroghSol.ddLDPO2(x)
            ddvRBC[i]           = self.kroghSol.ddvRBCConvectivePart(x)

        plt.figure()
        plt.plot(MValues, convectiveDrop, label='Convective drop of PO2')
        plt.plot(MValues, IVResistanceDrop, label='IV resistance drop of PO2')
        plt.xlabel(r'M_0')
        plt.legend()
        plt.title('PO2 drops for x = %g, LD = %g, vRBC = %g' % (x, LD, vRBC))

        plt.figure()
        plt.plot(MValues, LD*ddLD, label='Norm. deriv. w.r.t. LD')
        plt.plot(MValues, vRBC*ddvRBC, label='Norm. deriv. w.r.t. vRBC')
        plt.xlabel(r'M_0')
        plt.legend()
        plt.title('Normalized derivatives for x = %g, LD = %g, vRBC = %g' % (x, LD, vRBC))

    ## Plot the convective part and the IV resistance part of PO2 at the wall as a function
    #  oxygen consumption. Also plot the derivatives.
    def convectiveAndIVResistancePartsWRTx(self, LD, vRBC):
        xValues = np.linspace(0, self.kroghSol.L, 20)
        convectiveDrop   = np.zeros(xValues.shape)
        IVResistanceDrop = np.zeros(xValues.shape)
        ddLD             = np.zeros(xValues.shape)
        ddvRBC           = np.zeros(xValues.shape)
        for i, x in enumerate(xValues):
            convectiveDrop[i]   = self.kroghSol.convectivePO2Drop(x)
            IVResistanceDrop[i] = self.kroghSol.intravascResistancePO2Drop(x)
            ddLD[i]             = self.kroghSol.ddLDPO2(x)
            ddvRBC[i]           = self.kroghSol.ddvRBCConvectivePart(x)

        plt.figure()
        plt.plot(xValues, convectiveDrop, label='Convective drop of PO2')
        plt.plot(xValues, IVResistanceDrop, label='IV resistance drop of PO2')
        plt.xlabel(r'x')
        plt.legend()
        plt.title('PO2 drops for LD = %g, vRBC = %g' % (LD, vRBC))

        plt.figure()
        plt.plot(xValues, LD*ddLD, label='Norm. deriv. w.r.t. LD')
        plt.plot(xValues, vRBC*ddvRBC, label='Norm. deriv. w.r.t. vRBC')
        plt.xlabel(r'x')
        plt.legend()
        plt.title('Normalized derivatives for LD = %g, vRBC = %g' % (LD, vRBC))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paramFile',
                        help='Path to parameter study file', default='params.json')
    parser.add_argument('--probe', '-p', type=int, help='index of the probe to plot (starts with index 0)')
    parser.add_argument('--field', '-f', default='PO2', help='Field to plot (PO2, MTC or Nu)')
    parser.add_argument('--convO2Transport', '-c', action='store_true',
                        help='Whether to include oxygen convective transport in the model')
    figOptions = FigureOptions(parser)
    args = parser.parse_args()
    fileName = args.paramFile
    probeIdx = args.probe
    field = args.field
    convO2 = args.convO2Transport
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    paramStudy = LDvRBCParameterStudy(parameterStudyFile=fileName)
    postProcessor = make_param_study_post_processor(paramStudy)
    postProcessor.probeIdx = probeIdx
    plotParamStudy = plotLDvRBCParameterStudy(postProcessor)

    K0 = postProcessor.fitIntravascularResistance()
    print K0

    if field == 'PO2':
        plotParamStudy.contourPlot()
        plotParamStudy.contourPlotAnalyticalSolutionSecomb(K0=K0)
        K0 = plotParamStudy.kroghSol.intravascularResistanceLDHalf
        plotName = 'LD_U_PO2_contourPlot_K0_%g_probe%i' % (K0, probeIdx)
    elif field == 'MTC':
        plotParamStudy.contourPlotMTC()
        plotName = 'LD_U_PO2_contourPlot_MTC'

    figOptions.adjustAxes()
    figOptions.saveFig(plotName)
