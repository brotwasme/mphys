import numpy as np
from refnx.reflect import Component, SLD, ReflectModel, Structure, LipidLeaflet
from refnx.analysis import possibly_create_parameter, Parameters, Parameter
from protein.LipidLeafletWithProtein_builtOn2 import LipidLeafletWithProtien

class LipidLeafletProteinProtrusion(LipidLeafletWithProtien):
    def __init__(self, apm, b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono,

                 water_vm, waters_per_head, waters_per_tail,
                 vm_mscl, PLRatio,
                
                 head_solvent=None, tail_solvent=None,
                 reverse_monolayer=False, name=''):

        super(LipidLeafletProteinProtrusion, self).__init__(apm, b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono,

                 water_vm, waters_per_head, waters_per_tail,
                 vm_mscl, PLRatio,
                
                 head_solvent, tail_solvent,
                 reverse_monolayer, name)

        