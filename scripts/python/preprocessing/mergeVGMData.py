#!/usr/bin/env python
#
# Merge the graph pickle and the RBC pickle produced by VGM.
#
# Usage: mergeVGMData.py --graph <pathToGraphPickle>
#                        --RBC   <pathToRBCPickle>

from HbO2.setup.VGMDictMerger import VGMDictMerger

merger = VGMDictMerger()
merger.run()


