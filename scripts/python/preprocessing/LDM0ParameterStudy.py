#!/usr/bin/env python
#
# Encapsulates data that define a parameter study for linear density
# and oxygen consumption.
#

import numpy as np
import argparse
import json


class LDM0ParameterStudy:

    def __init__(self, parameterStudyFile='params.json'):
        self.parameterStudyFile = parameterStudyFile
        jsonData=open(parameterStudyFile)
        data = json.load(jsonData)
        jsonData.close()

        self.LDValues = np.linspace(data['LDMin'], data['LDMax'], data['nLD'])
        self.M0Values = np.linspace(data['M0Min'], data['M0Max'], data['nM0'])
        self.U = data['U']
        self.baseCaseName = data['baseCase']
        self.timeAverageStart = data['averageStart']
        self.LDomain = data['LDomain']

    @classmethod
    def fromparser(cls, parser):
        parser.add_argument('--paramFile', 
                help='Path to parameter study file', default='params.json')
        args = parser.parse_args()
        return cls(args.paramFile)

    def nLD(self):
        return len(self.LDValues)

    def nM0(self):
        return len(self.M0Values)

    @staticmethod
    def buildCaseName(LD, M0):
        return 'case_LD_%g_M0_%g' % (LD, M0)



