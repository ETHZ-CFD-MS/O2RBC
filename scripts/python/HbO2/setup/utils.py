import json
import os

graphDictPath = 'domain/constant/graphDict'
RBCPathParamsPath = 'RBCPathParams.json'


def isAxisymmetricCase(casePath):
    """
    Returns true if the current case is an axisymmetric case.
    """
    return not isGraphCase(casePath)


def isGraphCase(casePath):
    """
    Returns true if the current case is a graph case.
    """
    return os.path.exists(os.path.join(casePath, graphDictPath))


def isGraphStraightCapillariesCase(casePath):
    """
    Returns true if the current case is a graph case with an array of straight capillaries
    """
    return isGraphCase(casePath) and os.path.exists(os.path.join(casePath, RBCPathParamsPath))


class GraphInletValue(object):

    value_types = ['uniform', 'edge']

    def __init__(self, value_type, uniform_value=None, value_dict=None):
        super(GraphInletValue, self).__init__()
        if value_type not in self.value_types:
            raise ValueError('Invalid value_type {:s}'.format(value_type))
        self.value_type = value_type
        self.uniform_value = uniform_value
        self.value_dict = value_dict

    @classmethod
    def from_json(cls, path_to_json):
        fp = open(path_to_json)
        json_data = json.load(fp)
        fp.close()
        value_type = json_data['type']
        uniform_value = None
        value_dict = None
        if value_type == 'uniform':
            uniform_value = json_data['value']
        elif value_type == 'edge':
            value_list = json_data['values']
            value_dict = {ei: val for ei, val in value_list}
        return cls(value_type, uniform_value=uniform_value, value_dict=value_dict)

    def inlet_value(self, edge_index):
        if self.value_type == 'uniform':
            return self.uniform_value
        elif self.value_type == 'edge':
            return self.value_dict[edge_index]
