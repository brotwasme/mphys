import numpy as np
from refnx.reflect import Component, SLD, ReflectModel, Structure, LipidLeaflet
from refnx.analysis import possibly_create_parameter, Parameters, Parameter

class LipidLeafletTest(LipidLeaflet):
    def __init__(self, apm, b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono,

                #  water_vm, waters_per_head, waters_per_tail,
                #  vm_mscl,
                 protein_sld,
                 PLRatio,
                
                 head_solvent=None, tail_solvent=None,
                 reverse_monolayer=False, name=''):

        super(LipidLeafletTest, self).__init__( apm,
                 b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono,

                 #water_vm, waters_per_head, waters_per_tail,
                
                 head_solvent, tail_solvent,
                 reverse_monolayer, name)
        self.PLRatio = possibly_create_parameter(
                PLRatio,
                name='%s - PLRatio, fraction lipid to protein' % name)
        
        self.protein_sld = possibly_create_parameter(
                protein_sld,
                name='%s - protein sld' % name)
        self.protein_sld = SLD(protein_sld, name="protein sld")
    
    def slabs(self, structure=None):
        layers = super(LipidLeafletTest, self).slabs(structure=structure)

        if self.reverse_monolayer:
            layers = np.flipud(layers)
            layers[:, 3] = layers[::-1, 3]


        # tail region
        volfrac = self.vm_tails.value / (self.apm.value *
                                         self.thickness_tails.value)
        layers[0, 4] = 1 - volfrac
        if self.head_solvent is not None:
            # we do the solvation here, not in Structure.slabs
            layers[0] = Structure.overall_sld(layers[0], self.head_solvent)
            layers[0, 4] = 0

        # tail region
        volfrac = self.vm_tails.value / (self.apm.value *
                                         self.thickness_tails.value)

        layers[1, 4] = 1 - volfrac
        if self.tail_solvent is not None:
            # we do the solvation here, not in Structure.slabs
            layers[1] = Structure.overall_sld(layers[1], self.tail_solvent)
            layers[1, 4] = 0

        if self.reverse_monolayer:
            layers = np.flipud(layers)
            layers[:, 3] = layers[::-1, 3]

