#!/usr/bin/env python
#
# Plot RBC fields as a function of time.
#

import argparse

import matplotlib.pyplot as plt

from HbO2.postprocess.case import CasePostProcessor


def plotRBCFieldOverTime(fields_dict, fieldName, min_time, max_time):

    times = fields_dict['times']

    # time indices that should be plotted
    plotTimeIdx = [i for i,t in enumerate(times) if t >= min_time and
                                                    t <= max_time]

    if fieldName not in fields_dict:
        print 'Error, key %s is not in dictionary fields_dict' % fieldName
        return

    plt.plot(times[plotTimeIdx], fields_dict[fieldName][plotTimeIdx])
    plt.xlabel('time [s]')
    plt.ylabel(fieldName.replace('_', '\_'))
    plt.xlim([min_time, max_time])
    plt.grid(True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--field', '-f', nargs='+', 
                        help='Field to plot', default=['PO2'])
    parser.add_argument('--RBC', type=int, nargs='+', 
                        help='RBCs to plot', default=[0])
    parser.add_argument('--min-time', type=float, help='Smallest time that should be used', \
                        default = 0.0)
    parser.add_argument('--max-time', type=float, help='Largest time that should be used', \
                        default = 4.0)
    parser.add_argument('--show', help='Show the plot', action='store_true')

    args            = parser.parse_args()
    
    fieldNames      = args.field
    RBCIdx          = args.RBC
    min_time        = args.min_time
    max_time        = args.max_time
    showPlot        = args.show

    # treat special field names
    if fieldNames[0] == 'Hb':
        fieldNames = ['Hb_max', 'Hb_mean', 'Hb_min']
        fieldLabel = 'Hb'
    elif fieldNames[0] == 'PO2':
        fieldNames = ['PO2_max', 'PO2_mean', 'PO2_min']
        fieldLabel = 'PO2'
    else:
        fieldLabel = fieldNames[0]

    postProcessor = CasePostProcessor('.')
    RBCFieldDicts = postProcessor.rbc_data.RBCFields

    fig = plt.figure()

    for fieldName in fieldNames:
        for i in RBCIdx:
            plotRBCFieldOverTime(RBCFieldDicts[i], fieldName, min_time, max_time)

    plt.ylabel(fieldLabel.replace('_', '\_'))
    fig.tight_layout(pad=0.4)

    plotName = 'plot%sOverTime.png' % fieldLabel
    plt.savefig(plotName)

    if showPlot:
        plt.show()

