#!/usr/bin/env python
"""
Generic postprocessor for an OpenFOAM case.
"""

import argparse

from HbO2.postprocess.case import PostProcessorWriter
from HbO2.postprocess.factory.case import make_post_processor
from utilities.arguments import add_case_argument


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_case_argument(parser)
    parser.add_argument('--settingsFile', help='Relative path to the file with postprocessing settings',
                        default='postprocess.json')
    args = parser.parse_args()
    case_name = args.caseName
    settings_file = args.settingsFile
    postprocessor = make_post_processor(case_name, param_file=settings_file)
    writer = PostProcessorWriter(postprocessor)
    writer.write_results()
    print 'Wrote results to {:s}'.format(postprocessor.output_file_name)

