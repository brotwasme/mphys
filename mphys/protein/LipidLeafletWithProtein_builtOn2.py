import numpy as np
from refnx.reflect import Component, SLD, ReflectModel, Structure, LipidLeaflet
from refnx.analysis import possibly_create_parameter, Parameters, Parameter
from lipidBilayerAsGiven.LipidLeaflet_WaterPer_builtOn3 import LipidLeafletWithWaterPer

class LipidLeafletWithProtien(LipidLeafletWithWaterPer):
    def __init__(self, apm, b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono,

                 water_vm, waters_per_head, waters_per_tail,
                 vm_mscl, PLRatio,
                
                 head_solvent=None, tail_solvent=None,
                 reverse_monolayer=False, name=''):

        super(LipidLeafletWithProtien, self).__init__( apm,
                 b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono,

                 water_vm, waters_per_head, waters_per_tail,
                
                 head_solvent, tail_solvent,
                 reverse_monolayer, name)

#         if isinstance(b_mscl, complex):
#             self.b_mscl_real = possibly_create_parameter(
#                 b_mscl.real,
#                 name='%s - b_mscl_real' % name)
#             self.b_mscl_imag = possibly_create_parameter(
#                 b_mscl.imag,
#                 name='%s - b_mscl_imag' % name)
#         elif isinstance(b_mscl, SLD):
#             self.b_mscl_real = b_mscl.real
#             self.b_mscl_imag = b_mscl.imag
#         else:
#             self.b_mscl_real = possibly_create_parameter(
#                 b_mscl,
#                 name='%s - b_mscl_real' % name)
#             self.b_mscl_imag = possibly_create_parameter(
#                 0,
#                 name='%s - b_mscl_imag' % name)
        
        self.vm_mscl = possibly_create_parameter(
                vm_mscl,
                name='%s - vm_mscl, protien volume' % name)
        
        self.PLRatio = possibly_create_parameter(
                PLRatio,
                name='%s - PLRatio, protien lipid ratio' % name)

    def b_mscl_head(self):
        bc = 0.6646e-4  #Carbon
        bo = 0.5804e-4  #Oxygen
        bh = -0.3739e-4 #Hydrogen
        #bp = 0.513e-4   #Phosphorus
        bn = 0.936e-4   #Nitrogen
        bd = 0.6671e-4  #Deuterium
        bs = 2.847e-4   #Sulphur
        mol_frac,i = self.d2o_mol_fraction_head()
        return float(2745*bc + 675*bn + 641*bo + 25*bs + 
                    (4374.5-822.5*0.9)*bh +
                    mol_frac*822.5*0.9*bd + 
                    (1-mol_frac)*822.5*0.9*bh)

    def b_mscl_tail(self):
        bc = 0.6646e-4  #Carbon
        bo = 0.5804e-4  #Oxygen
        bh = -0.3739e-4 #Hydrogen
        #bp = 0.513e-4   #Phosphorus
        bn = 0.936e-4   #Nitrogen
        bd = 0.6671e-4  #Deuterium
        bs = 2.847e-4   #Sulphur
        mol_frac,i = self.d2o_mol_fraction_tail()
#         print("mol_frac",mol_frac)
        return float(2745*bc + 675*bn + 641*bo + 25*bs + 
                    (4374.5-822.5*0.9)*bh +
                    mol_frac*822.5*0.9*bd + 
                    (1-mol_frac)*822.5*0.9*bh)


    def sld_r(self):
        lipid_head_sld = float(self.b_heads_real) / float(self.vm_head()) * 1.e6
        lipid_tail_sld = float(self.b_tails_real) / float(self.vm_tail()) * 1.e6
        
        mscl_h_sld_r = self.b_mscl_head()/float(self.vm_mscl.value) * 1.e6
        mscl_t_sld_r = self.b_mscl_tail()/float(self.vm_mscl.value) * 1.e6

#         mscl_head_sld = mscl_h_sld_r*(float(self.thickness_heads)/
#                                   (2*(float(self.thickness_heads)+float(self.thickness_tails))))
#         mscl_tail_sld = mscl_t_sld_r*(float(self.thickness_heads)/
#                                   (2*(float(self.thickness_heads)+float(self.thickness_tails))))
#         print("protein:",mscl_h_sld_r,mscl_t_sld_r)
        head_sld = (self.PLRatio*lipid_head_sld) + (1-self.PLRatio)*mscl_h_sld_r
        tail_sld = (self.PLRatio*lipid_tail_sld) + (1-self.PLRatio)*mscl_t_sld_r
#         print("slds: ",head_sld,tail_sld)
        return float(head_sld), float(tail_sld)

    def sld_i(self):
        lipid_head_sld = float(self.b_heads_imag) / float(self.vm_head()) * 1.e6
        lipid_tail_sld = float(self.b_tails_imag) / float(self.vm_tail()) * 1.e6

# #         b_mscl_head_r, b_mscl_head_i = self.b_mscl_head()
# #         b_mscl_tail_r, b_mscl_tail_i = self.b_mscl_tail()
        
#         mscl_h_sld_i = float(b_mscl_head_i)/float(self.vm_mscl.value)# * 1.e6
#         mscl_t_sld_i = float(b_mscl_tail_i)/float(self.vm_mscl.value)# * 1.e6

#         mscl_head_sld = mscl_h_sld_i*(float(self.thickness_heads)/
#                                   (2*(float(self.thickness_heads)+float(self.thickness_tails))))
#         mscl_tail_sld = mscl_t_sld_i*(float(self.thickness_heads)/
#                                   (2*(float(self.thickness_heads)+float(self.thickness_tails))))

        #head_sld = lipid_head_sld#(self.PLRatio*lipid_head_sld) + (1-self.PLRatio)*0#mscl_head_sld
        #tail_sld = lipid_tail_sld#`(self.PLRatio*lipid_tail_sld) + (1-self.PLRatio)*0#mscl_tail_sld
        return lipid_head_sld, lipid_tail_sld


    @property
    def parameters(self):
        p = super(LipidLeafletWithProtien, self).parameters
        p.extend([self.vm_mscl,
                    self.PLRatio])
        return p
