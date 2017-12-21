#!/usr/bin/env python
"""
Translates the data produced by VGM to OpenFOAM dictionaries.

Usage: translateVGMData.py <pathToPickleFile>
"""

from HbO2.setup.VGMDataTranslator import VGMDataTranslator

translator = VGMDataTranslator()
translator.run()
