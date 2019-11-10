import numpy as np
from refnx.reflect import Component, SLD, ReflectModel, Structure
from refnx.analysis import possibly_create_parameter, Parameters, Parameter

class NLayerModel(Component):
    def __init__(self, name, bs, thicks, roughs, totalThickness):
        super(NLayerModel, self).__init__()
        self.Nlayers = len(bs)
        self.bs = []
        self.thicks = []
        self.roughs = []
        self.name = name
        for i in range(self.Nlayers):
            self.bs.append( possibly_create_parameter(bs[i],
                name='%s - layer %d b_real' % (name,i)) )
            self.thicks.append( possibly_create_parameter(thicks[i],
                name='%s - thickness %d thickness' % (name,i)) )
            self.roughs.append( possibly_create_parameter(roughs[i],
                name='%s - layer %d roughness' % (name,i)) )
        self.totalThickness = possibly_create_parameter(totalThickness,
            name='%s - total thickness' % (name))

    def slabs(self, structure=None):
        layers = np.zeros((self.Nlayers,5))
        for i in range(self.Nlayers):
            layers[i, 0] = float(self.thicks[i])
            layers[i, 1] = float(self.bs[i])
            layers[i, 2] = float(0)
            layers[i, 3] = float(self.roughs[i])
            layers[i, 4] = float(0)
        return layers

    @property
    def parameters(self):
        p = Parameters(name=self.name)
        p.extend([self.bs,
                 self.thicks,
                 self.roughs,
                 self.totalThickness])
        return p

    def logp(self):
        returns = 0
        calcThickness = 0
        for thick in self.thick:
            calcThickness += float(thick)
        if self.totalThickness>calcThickness:
            returns = -np.inf
        return returns


