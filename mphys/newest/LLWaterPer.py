import numpy as np
from refnx.reflect import Component, SLD, ReflectModel, Structure, LipidLeaflet
from refnx.analysis import possibly_create_parameter, Parameters, Parameter

class LipidLeafletWithWaterPer(LipidLeaflet):
    def __init__(self, apm, b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono,

                 water_vm, waters_per_head, waters_per_tail,
                
                 head_solvent=None, tail_solvent=None,
                 reverse_monolayer=False, name=''):

        super(LipidLeafletWithWaterPer, self).__init__(apm, b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono, head_solvent, tail_solvent,
                 reverse_monolayer, name)

        self.water_vm = possibly_create_parameter(
                water_vm,
                name='%s - water_vm' % name)

        self.waters_per_head = possibly_create_parameter(
                waters_per_head,
                name='%s - waters_per_head' % name)
        
        self.waters_per_tail = possibly_create_parameter(
                waters_per_tail,
                name='%s - waters_per_tail' % name)
        
        bo = 0.5804e-4     #Oxygen
        bh = -0.3741e-4    #Hydrogen
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

    def wMol_calc(self, r):#, i):
#         print(self.name,"wMol_calc",r,i)
        return (r * self.D2O) + ((1-r)*self.H2O)#, (i * self.D2O) + ((1-i)*self.H2O)

    def wMol_head(self):
        d2o_mol_fraction_head_r = self.d2o_mol_fraction_head() #, d2o_mol_fraction_head_i 
        return self.wMol_calc(d2o_mol_fraction_head_r) #, d2o_mol_fraction_head_i)

    def wMol_tail(self):
        d2o_mol_fraction_tail_r = self.d2o_mol_fraction_tail() # , d2o_mol_fraction_tail_i 
        return self.wMol_calc(d2o_mol_fraction_tail_r) #, d2o_mol_fraction_tail_i)

    def calc_solvent_true(self, r):#, i):
#         print(self.name,"calc_solvent_true",r,i)
        return complex(r/self.water_vm.value,0)# i/self.water_vm.value)

    def head_solvent_true(self):
        wMol_r = self.wMol_head() #, wMol_i 
#         print("wMol_r, wMol_i:", wMol_r, wMol_i)
        return self.calc_solvent_true(wMol_r) #, wMol_i)

    def tail_solvent_true(self):
        wMol_r = self.wMol_tail() #, wMol_i
#         print("wMol_tail_r, wMol_tail_i:", wMol_tail_r, wMol_tail_i)
        return self.calc_solvent_true(wMol_r)#, wMol_i)

    def vm_head(self):
        return self.vm_heads.value

    def vm_tail(self):
        return self.vm_tails.value

    def total_vm(self):
        vm_head_lipidAndsolvent = self.vm_head() + self.water_vm.value * self.waters_per_head.value
        vm_tail_lipidAndsolvent = self.vm_tail() + self.water_vm.value * self.waters_per_tail.value
        return 2*(vm_head_lipidAndsolvent+vm_tail_lipidAndsolvent)
    
    def total_thickness(self):
        th_heads = self.thickness_heads.value
        th_tails = self.thickness_tails.value
        return 2*(th_heads+th_tails)
    
    def sld_r(self):
        head_sld = float(self.b_heads_real) / float(self.vm_heads.value) * 1.e6
        tail_sld = float(self.b_tails_real) / float(self.vm_tails.value) * 1.e6
        return head_sld, tail_sld

    def sld_i(self):
        head_sld = float(self.b_heads_imag) / float(self.vm_heads.value) * 1.e6
        tail_sld = float(self.b_tails_imag) / float(self.vm_tails.value) * 1.e6
        return head_sld, tail_sld

    def slabs(self, structure=None):
        """
        Slab representation of monolayer, as an array
        Parameters
        ----------
        structure : refnx.reflect.Structure
            The Structure hosting this Component
        """
        layers = np.zeros((2, 5))

        # thicknesses
        layers[0, 0] = float(self.thickness_heads)
        layers[1, 0] = float(self.thickness_tails)

        # real and imag SLD's
        head_sld_r, tail_sld_r = self.sld_r()
        head_sld_i, tail_sld_i = self.sld_i()
        layers[0, 1] = head_sld_r
        layers[0, 2] = head_sld_i

        layers[1, 1] = tail_sld_r
        layers[1, 2] = tail_sld_i

        # roughnesses
        layers[0, 3] = float(self.rough_preceding_mono)
        layers[1, 3] = float(self.rough_head_tail)

        # volume fractions
        # head region
        
        apm = self.total_vm()/self.total_thickness()
        
        self.apm.value = apm
        
        volfrac = self.vm_head() / (self.apm.value *
                                         self.thickness_heads.value)
        # volfrac = self.vm_head() / (self.vm_head() + self.water_vm.value * self.waters_per_head.value)
        layers[0, 4] = 1 - volfrac
        self.vh = volfrac
        # print("volfrac h", volfrac)
        # print(layers[0])
        if self.head_solvent is not None:
            # we do the solvation here, not in Structure.slabs
            layers[0] = Structure.overall_sld(layers[0], self.head_solvent_true())
            layers[0, 4] = 0

        # print(layers[0])
        # tail region
        volfrac = self.vm_tail() / (self.apm.value *
                                         self.thickness_tails.value)
        # volfrac = self.vm_tail() / (self.vm_tail() + self.water_vm.value * self.waters_per_tail.value)

        layers[1, 4] = 1 - volfrac
        self.vt = volfrac
        # print("volfrac t", volfrac)
        if self.tail_solvent is not None:
            # we do the solvation here, not in Structure.slabs
            layers[1] = Structure.overall_sld(layers[1], self.tail_solvent_true())
            layers[1, 4] = 0

        if self.reverse_monolayer:
            layers = np.flipud(layers)
            layers[:, 3] = layers[::-1, 3]
#         print("layers",layers)
        return layers

    @property
    def parameters(self):
        p = super(LipidLeafletWithWaterPer, self).parameters
        
        p.extend([
            self.water_vm, self.waters_per_head, self.waters_per_tail])
        if self.waters_per_head is not None:
            p.extend([self.waters_per_head])
        if self.waters_per_tail is not None:
            p.extend([self.waters_per_tail])
        return p

    def logp(self):
        returns = 0
        # penalise unphysical volume fractions.
#         volfrac_h = self.vm_head() / (self.apm.value *
#                                     self.thickness_heads.value)

#         # tail region
#         volfrac_t = self.vm_tail() / (self.apm.value *
#                                     self.thickness_tails.value)

#         if volfrac_h > 1 or volfrac_t > 1:
#             returns += -np.inf
        # if self.thickness_heads.value > self.thickness_tails.value:
        #     returns += -np.inf
        return returns
