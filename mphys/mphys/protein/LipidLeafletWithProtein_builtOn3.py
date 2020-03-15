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
                name='%s - PLRatio, fraction lipid to protein' % name)
        
        bc = 0.6646e-4  #Carbon
        bo = 0.5804e-4  #Oxygen
        bh = -0.3739e-4 #Hydrogen
        #bp = 0.513e-4   #Phosphorus
        bn = 0.936e-4   #Nitrogen
        #bd = 0.6671e-4  #Deuterium
        bs = 2.847e-4   #Sulphur
        self.b_mscl_base = float(2745*bc + 675*bn + 641*bo + 25*bs + 
                    (4374.5-822.5*0.9)*bh)

#     def vm_head(self):
        

    def b_mscl_head(self):
        #bc = 0.6646e-4  #Carbon
        #bo = 0.5804e-4  #Oxygen
        bh = -0.3739e-4 #Hydrogen
        #bp = 0.513e-4   #Phosphorus
        #bn = 0.936e-4   #Nitrogen
        bd = 0.6671e-4  #Deuterium
        #bs = 2.847e-4   #Sulphur
        mol_frac = self.d2o_mol_fraction_head() #,i
        return float(self.b_mscl_base +
                    mol_frac*822.5*0.9*bd
                    +(1-mol_frac)*822.5*0.9*bh)

    def b_mscl_tail(self):
        #bc = 0.6646e-4  #Carbon
        #bo = 0.5804e-4  #Oxygen
        bh = -0.3739e-4 #Hydrogen
        #bp = 0.513e-4   #Phosphorus
       # bn = 0.936e-4   #Nitrogen
        bd = 0.6671e-4  #Deuterium
       # bs = 2.847e-4   #Sulphur
        mol_frac = self.d2o_mol_fraction_tail() #,i
#         print("mol_frac",mol_frac)
        return float(self.b_mscl_base +
                    mol_frac*822.5*0.9*bd
                    +(1-mol_frac)*822.5*0.9*bh)

#     def vm_head(self):
# #         lipidAndWaterHeadVolume = super(LipidLeafletWithProtien, self).vm_head()
# #         lipidAndWaterTailVolume = super(LipidLeafletWithProtien, self).vm_tail()
#         lipidAndWaterHeadVolume = self.PLRatioi.value*self.vm_head.value + (1-self.PLRatioi.value)*self.vm_mscl.value
#         lipidAndWaterTailVolume = 
#         return 0

    def protien_frac_in_head(self):
        return self.thickness_heads.value/self.total_thickness()

    def protien_frac_in_tail(self):
        return self.thickness_tails.value/self.total_thickness()

    def vm_head(self):
        # print("head vm")
        return self.vm_heads.value*(self.PLRatio.value) + (self.vm_mscl.value*self.protien_frac_in_head())*(1-self.PLRatio.value)# + self.water_vm.value * self.waters_per_head.value

    def vm_tail(self):
        # print("tail vm")
        return self.vm_tails.value*(self.PLRatio.value) + (self.vm_mscl.value*self.protien_frac_in_tail())*(1-self.PLRatio.value)# + self.water_vm.value * self.waters_per_tail.value

    def total_vm(self):
        vm_lipid = 2*(self.vm_heads.value+self.vm_tails.value)
        vm_protienAndlipid = vm_lipid*(self.PLRatio.value) + self.vm_mscl.value*(1-self.PLRatio.value)
        # print("volums:", vm_protienAndlipid, 2*(self.vm_head()+self.vm_tail()))
        vm_water = 2*self.water_vm.value * (self.waters_per_head.value+self.waters_per_tail.value)
        # print("total vm")
        return vm_protienAndlipid + vm_water
        

    def sld_r(self):
        lipid_head_sld, lipid_tail_sld = super(LipidLeafletWithProtien, self).sld_r()
        
        mscl_h_sld_r = self.b_mscl_head()/float(self.vm_mscl.value) * 1.e6 #if solvation in head and tail are same then so are these, ignores water per 
        mscl_t_sld_r = self.b_mscl_tail()/float(self.vm_mscl.value) * 1.e6
        
        head_sld = (self.PLRatio*lipid_head_sld) + (1-self.PLRatio)*mscl_h_sld_r
        tail_sld = (self.PLRatio*lipid_tail_sld) + (1-self.PLRatio)*mscl_t_sld_r
#         print("slds: ",head_sld,tail_sld)
        return float(head_sld), float(tail_sld)

    def sld_i(self):
        lipid_head_sld, lipid_tail_sld = super(LipidLeafletWithProtien, self).sld_i()
        return lipid_head_sld, lipid_tail_sld


    @property
    def parameters(self):
        p = super(LipidLeafletWithProtien, self).parameters
        p.extend([self.vm_mscl,
                    self.PLRatio])
        return p
