import numpy as np
from refnx.reflect import Component, SLD, ReflectModel, Structure
from refnx.analysis import possibly_create_parameter, Parameters, Parameter

class nLayerClass(Component):
    def __init__(self, layerNumber=None,bs=[],thicks=[],roughs=[],name=[]):
        super(nLayerClass,self).__init__()
        self.nLayers=layerNumber

        if layerNumber>=1:
            self.b1 = possibly_create_parameter(bs[i],
                name='%s - layer %d b_real' % (name,i))
            self.thick1=possibly_create_parameter(thicks[i],
                name='%s - thickness %d thickness' % (name,i))
            self.rough1=possibly_create_parameter(roughs[i],
                name='%s - layer %d roughness' % (name,i))

        if layerNumber>=2:
            self.b2=possibly_create_parameter(bs[i],
                name='%s - layer %d b_real' % (name,i))
            self.thick2=possibly_create_parameter(thicks[i],
                name='%s - thickness %d thickness' % (name,i))
            self.rough2=possibly_create_parameter(roughs[i],
                name='%s - layer %d roughness' % (name,i))

        if layerNumber>=3:
            self.b3=possibly_create_parameter(bs[i],
                name='%s - layer %d b_real' % (name,i))
            self.thick3=possibly_create_parameter(thicks[i],
                name='%s - thickness %d thickness' % (name,i))
            self.rough3=possibly_create_parameter(roughs[i],
                name='%s - layer %d roughness' % (name,i))

        if layerNumber>=4:
            self.b4=possibly_create_parameter(bs[i],
                name='%s - layer %d b_real' % (name,i))
            self.thick4=possibly_create_parameter(thicks[i],
                name='%s - thickness %d thickness' % (name,i))
            self.rough4=possibly_create_parameter(roughs[i],
                name='%s - layer %d roughness' % (name,i))

        if layerNumber>=5:
            self.b5=possibly_create_parameter(bs[i],
                name='%s - layer %d b_real' % (name,i))
            self.thick5=possibly_create_parameter(thicks[i],
                name='%s - thickness %d thickness' % (name,i))
            self.rough5=possibly_create_parameter(roughs[i],
                name='%s - layer %d roughness' % (name,i))


        self.totalThickness = possibly_create_parameter(totalThickness,
            name='%s - total thickness' % (name))

    def slabs(self, structure=None):
        layerNumber = self.nLayers
        layers = np.zeros((self.nLayers,5))
        i=1
        if layerNumber>=i:
            layers[i, 0] = float(self.thick1)
            layers[i, 1] = float(self.b1)
            layers[i, 2] = float(0)
            layers[i, 3] = float(self.rough1)
            layers[i, 4] = float(0)
        i=2
        if layerNumber>=i:
            layers[i, 0] = float(self.thick2)
            layers[i, 1] = float(self.b2)
            layers[i, 2] = float(0)
            layers[i, 3] = float(self.rough2)
            layers[i, 4] = float(0)
        i=3
        if layerNumber>=i:
            layers[i, 0] = float(self.thick3)
            layers[i, 1] = float(self.b3)
            layers[i, 2] = float(0)
            layers[i, 3] = float(self.rough3)
            layers[i, 4] = float(0)
        i=4
        if layerNumber>=i:
            layers[i, 0] = float(self.thick4)
            layers[i, 1] = float(self.b4)
            layers[i, 2] = float(0)
            layers[i, 3] = float(self.rough4)
            layers[i, 4] = float(0)
        i=5
        if layerNumber>=i:
            layers[i, 0] = float(self.thick5)
            layers[i, 1] = float(self.b5)
            layers[i, 2] = float(0)
            layers[i, 3] = float(self.rough5)
            layers[i, 4] = float(0)
        return layers

    @property
    def parameters(self):
        layerNumber = self.nLayers
        p = Parameters(name=self.name)
        i=1
        if layerNumber>=i:
            p.extend([self.b1,
                     self.thick1,
                     self.rough1])
        i=2
        if layerNumber>=i:
            p.extend([self.b2,
                     self.thick2,
                     self.rough2])
        i=3
        if layerNumber>=i:
            p.extend([self.b3,
                     self.thick3,
                     self.rough3])
        i=4
        if layerNumber>=i:
            p.extend([self.b4,
                     self.thick4,
                     self.rough4])
        i=5
        if layerNumber>=i:
            p.extend([self.b5,
                     self.thick5,
                     self.rough5])

        p.extend([self.totalThickness])
        return p

    def logp(self):
        layerNumber = self.nLayers
        returns = 0
        calcThickness = 0
        i=1
        if layerNumber>=i:
            calcThickness += float(thick1)
        i=2
        if layerNumber>=i:
            calcThickness += float(thick2)
        i=3
        if layerNumber>=i:
            calcThickness += float(thick3)
        i=4
        if layerNumber>=i:
            calcThickness += float(thick4)
        i=5
        if layerNumber>=i:
            calcThickness += float(thick5)
#         for thick in self.thick:
#             calcThickness += float(thick)
        if self.totalThickness>calcThickness:
            returns = -np.inf
        return returns

