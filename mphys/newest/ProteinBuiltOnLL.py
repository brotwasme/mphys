import numpy as np
from refnx.reflect import Component, SLD, ReflectModel, Structure, LipidLeaflet
from refnx.analysis import possibly_create_parameter, Parameters, Parameter
from lipidBilayerAsGiven.LipidLeaflet_WaterPer_builtOn4 import LipidLeafletWithWaterPer

class LipidLeafletWithProtien(LipidLeaflet):
    def __init__(self, apm, b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono,

                #  water_vm, waters_per_head, waters_per_tail,
                 vm_mscl, PLRatio,
                
                 head_solvent=None, tail_solvent=None,
                 reverse_monolayer=False, name=''):
        
        self.vm_mscl = possibly_create_parameter(
                vm_mscl,
                name='%s - vm_mscl, protien volume' % name)
        
        self.PLRatio = possibly_create_parameter(
                PLRatio,
                name='%s - PLRatio, fraction lipid to protein' % name)

        super(LipidLeafletWithProtien, self).__init__( apm,
                 b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono,

                #  water_vm, waters_per_head, waters_per_tail,
                
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
        
        bc = 0.6646e-4  #Carbon
        bo = 0.5804e-4  #Oxygen
        bh = -0.3739e-4 #Hydrogen
        #bp = 0.513e-4   #Phosphorus
        bn = 0.936e-4   #Nitrogen
        #bd = 0.6671e-4  #Deuterium
        bs = 2.847e-4   #Sulphur
        self.b_mscl_base = float(2745*bc + 675*bn + 641*bo + 25*bs + 
                    (4374.5-822.5*0.9)*bh)
        # bo = 0.5804e-4     #Oxygen
        # bh = -0.3741e-4    #Hydrogen
        bd = 0.6671e-4     #Deuterium
        D2O = (2*bd) + (1*bo)
        H2O = (2*bh) + (1*bo)
        self.D2O = float(D2O)
        self.H2O = float(H2O)

    def d2o_mol_fraction_calc(self, r): #, i):
        r = r.value*1e-6#, i.value*1e-6 ,i
#         print(self.name,"d2o_mol_fraction_calc", r, i)
        return (1/self.D2O-self.H2O)*((r*27.64)-self.H2O)#, (1/self.D2O-self.H2O)*((i*27.64)-self.H2O)
    
    def d2o_mol_fraction_head(self):
        return self.d2o_mol_fraction_calc(self.head_solvent.real)#, self.head_solvent.imag)

    def d2o_mol_fraction_tail(self):
        return self.d2o_mol_fraction_calc(self.tail_solvent.real)#, self.tail_solvent.imag)

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
        # print("mol_frac",mol_frac)
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

    def calc_protein_SLD(self, b):
        return complex(b*1.e6/self.vm_mscl.value, 0)

    def protein_head_SLD(self):
        return  self.calc_protein_SLD(self.b_mscl_head())

    def protein_tail_SLD(self):
        return  self.calc_protein_SLD(self.b_mscl_tail())

    def total_thickness(self):
        th_heads = self.thickness_heads.value
        th_tails = self.thickness_tails.value
        return 2*(th_heads+th_tails)

    def slabs(self, structure=None):
        """
        Slab representation of monolayer, as an array

        Parameters
        ----------
        structure : refnx.reflect.Structure
            The Structure hosting this Component
        """
        layers = super(LipidLeafletWithProtien, self).slabs(structure=structure)

        if self.reverse_monolayer:
            layers = np.flipud(layers)
            layers[:, 3] = layers[::-1, 3]
        
        # volume fractions
        # head region
        # volfrac = self.vm_heads.value / (self.apm.value *
        #                                  self.thickness_heads.value)
        layers[0, 4] = (1 - self.PLRatio.value)*self.thickness_heads.value/self.total_thickness()
        if self.protein_head_SLD() is not None:
            # we do the solvation here, not in Structure.slabs
            layers[0] = Structure.overall_sld(layers[0], self.protein_head_SLD())
            layers[0, 4] = 0

        # tail region
        # volfrac = self.vm_tails.value / (self.apm.value *
        #                                  self.thickness_tails.value)

        layers[1, 4] = (1 - self.PLRatio.value)*self.thickness_heads.value/self.total_thickness()
        if self.protein_tail_SLD() is not None:
            # we do the solvation here, not in Structure.slabs
            layers[1] = Structure.overall_sld(layers[1], self.protein_tail_SLD())
            layers[1, 4] = 0

        if self.reverse_monolayer:
            layers = np.flipud(layers)
            layers[:, 3] = layers[::-1, 3]
        
        return layers


    @property
    def parameters(self):
        p = super(LipidLeafletWithProtien, self).parameters
        p.extend([self.vm_mscl,
                    self.PLRatio])
        return p

    def logp(self):
        returns = 0
        returns += super(LipidLeafletWithProtien, self).logp()

        frac_h = self.d2o_mol_fraction_head()
        frac_t = self.d2o_mol_fraction_tail()

        if frac_h > 1 or  frac_t > 1 or frac_h < 0 or  frac_t < 0:
            returns += -np.inf
        # # penalise unphysical volume fractions.
        # volfrac_h = self.vm_head() / (self.apm.value *
        #                             self.thickness_heads.value)

        # # tail region
        # volfrac_t = self.vm_tail() / (self.apm.value *
        #                             self.thickness_tails.value)

        # if volfrac_h > 1 or volfrac_t > 1:
        #     returns += -np.inf
        if self.thickness_heads.value > self.thickness_tails.value:
            returns += -np.inf
        return returns