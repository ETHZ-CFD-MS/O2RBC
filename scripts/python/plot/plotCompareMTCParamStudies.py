#!/usr/bin/env python
#
# Plot on the same graph MTC from different parameter studies
#

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np

from HbO2.plot.labels import setXLabel, setYLabel
from HbO2.postprocess.parameterstudy import LDvRBCParameterStudyPostProcessor
from HbO2.setup.parameterStudy import LDvRBCParameterStudy
from plot.figureoptions import FigureOptions

pathToFolder = "/local/aluecker/OpenFOAM/aluecker-2.3.0/run/cbf/HbO2/eulerAxisymmetric/LD_U_study"
paramStudies = ['param_studyCortexFineNew', 'param_studyGlomerulusFineNew']

def plotMTC(postProcessor, U_plot, **kwargs):
    linestyle = kwargs.get('linestyle', '-')
    data = np.loadtxt(postProcessor.MTCOutputFilePath(),
                      delimiter="\t", skiprows=1)
    LDValues = postProcessor.param_study['LD']
    UValues = postProcessor.param_study['U']
    n_LD = len(postProcessor.param_study['LD'])
    n_U = len(postProcessor.param_study['U'])
    MTC = np.reshape(data[:,2], (n_LD, n_U))
    Nu  = np.reshape(data[:,3], (n_LD, n_U))

    i_U = np.argmin((postProcessor.param_study['U'] - U_plot)**2)
    plt.plot(LDValues, 1e2*MTC[:, i_U], linestyle=linestyle)
    print 1e2*MTC[:, i_U]
    setXLabel('LD', '-')
    setYLabel('k', r'10^6 \mathrm{ml\,O_2\,cm^{-2}\,s^{-1}\,mm\,Hg^{-1}}')
    plt.xlim(min(LDValues), max(LDValues))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-U', type=float, help='RBC velocity', default=1e-3)
    figOptions = FigureOptions(parser)
    args = parser.parse_args()
    U_plot = args.U
    figOptions.parseOptions(args)

    plotName = 'plotMTCCortexGlomerulus_U_{:g}'.format(U_plot)
    styles = [{'linestyle': '-'},
              {'linestyle': '--'}]

    for param_study_name, style in zip(paramStudies, styles):
        param_study_path = os.path.join(pathToFolder, param_study_name)
        param_study = LDvRBCParameterStudy(os.path.join(param_study_path, 'params.json'))
        postProcessor = LDvRBCParameterStudyPostProcessor(param_study)
        plotMTC(postProcessor, U_plot, **style)

    figOptions.adjustAxes()
    figOptions.saveFig(plotName)
