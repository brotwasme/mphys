import numpy as np
from refnx.reflect import Component, SLD, ReflectModel, Structure, LipidLeaflet
from refnx.analysis import possibly_create_parameter, Parameters, Parameter
from lipidBilayerAsGiven.LipidLeaflet_WaterPer_builtOn import LipidLeafletWithWaterPer

class LipidLeafletWithProtien(LipidLeafletWithWaterPer):
    def __init__(self, apm, b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono,

                 water_vm, waters_per_head, waters_per_tail, b_mscl, vm_mscl, PLRatio,
                
                 head_solvent=None, tail_solvent=None,
                 reverse_monolayer=False, name=''):

        super(LipidLeafletWithProtien, self).__init__( apm, b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono,

                 water_vm, waters_per_head, waters_per_tail,
                
                 head_solvent, tail_solvent,
                 reverse_monolayer, name)


        if isinstance(b_mscl, complex):
            self.b_mscl_real = possibly_create_parameter(
                b_mscl.real,
                name='%s - b_mscl_real' % name)
            self.b_mscl_imag = possibly_create_parameter(
                b_mscl.imag,
                name='%s - b_mscl_imag' % name)
        elif isinstance(b_mscl, SLD):
            self.b_mscl_real = b_mscl.real
            self.b_mscl_imag = b_mscl.imag
        else:
            self.b_mscl_real = possibly_create_parameter(
                b_mscl,
                name='%s - b_mscl_real' % name)
            self.b_mscl_imag = possibly_create_parameter(
                0,
                name='%s - b_mscl_imag' % name)
        
        self.vm_mscl = possibly_create_parameter(
                vm_mscl,
                name='%s - vm_mscl, protien volume' % name)
        
        self.PLRatio = possibly_create_parameter(
                PLRatio,
                name='%s - PLRatio, protien lipid ratio' % name)


    def sld_r(self):
        lipid_head_sld = float(self.b_heads_real) / self.vm_head() * 1.e6
        lipid_tail_sld = float(self.b_tails_real) / self.vm_tail() * 1.e6

        mscl_sld = float(self.b_mscl_real)/float(self.vm_mscl.value) * 1.e6

        mscl_head_sld = mscl_sld*(float(self.thickness_heads)/
                                  (2*(float(self.thickness_heads)+float(self.thickness_tails))))
        mscl_tail_sld = mscl_sld*(float(self.thickness_heads)/
                                  (2*(float(self.thickness_heads)+float(self.thickness_tails))))

        head_sld = (self.PLRatio*lipid_head_sld) + (1-self.PLRatio)*mscl_head_sld
        tail_sld = (self.PLRatio*lipid_tail_sld) + (1-self.PLRatio)*mscl_tail_sld
        return head_sld, tail_sld

    def sld_i(self):
        lipid_head_sld = float(self.b_heads_imag) / self.vm_head() * 1.e6
        lipid_tail_sld = float(self.b_tails_imag) / self.vm_tail() * 1.e6

        mscl_sld = float(self.b_mscl_imag)/float(self.vm_mscl.value) * 1.e6

        mscl_head_sld = mscl_sld*(float(self.thickness_heads)/
                                  (2*(float(self.thickness_heads)+float(self.thickness_tails))))
        mscl_tail_sld = mscl_sld*(float(self.thickness_heads)/
                                  (2*(float(self.thickness_heads)+float(self.thickness_tails))))

        head_sld = (self.PLRatio*lipid_head_sld) + (1-self.PLRatio)*mscl_head_sld
        tail_sld = (self.PLRatio*lipid_tail_sld) + (1-self.PLRatio)*mscl_tail_sld
        return head_sld, tail_sld


    @property
    def parameters(self):
        p = super(LipidLeafletWithProtien, self).parameters
        p.extend([self.b_mscl_real, self.b_mscl_imag, self.vm_mscl,
                    self.PLRatio])
        return p
