"""
Postprocessing of topological and functional tissue volumes
"""


from HbO2.parse.readhemoglobin import TissueVolumesReader
from HbO2.postprocess.results import HemoglobinOnSegmentsResults


class TissueVolumeResults(HemoglobinOnSegmentsResults):

    default_tissue_volume_file = 'topologicalTissueVolumes.txt'

    def __init__(self, case_path, **kwargs):
        super(TissueVolumeResults, self).__init__(case_path)
        self.tissue_volume_reader = TissueVolumesReader(
            self.case_path,
            kwargs.get('tissueVolumesFile', self.default_tissue_volume_file)
        )

    def topological_tissue_volumes(self):
        return self['topolTissueVolume']

    def topological_tissue_radii(self):
        return self['topolTissueRadius']

    def functional_tissue_volumes(self):
        return self['funcTissueVolume']

    def functional_tissue_radii(self):
        return self['funcTissueRadius']*(self['funcTissueRadius'] >= 0)

    def functional_tissue_radii_with_negative_values(self):
        return self['funcTissueRadius']

    def functional_minus_topological_tissue_radii(self):
        return self.functional_tissue_radii() - self.topological_tissue_radii()
