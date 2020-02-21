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

    def d2o_mol_fraction_calc(self, r, i):
        r,i = r.value*1e-6, i.value*1e-6
#         print(self.name,"d2o_mol_fraction_calc", r, i)
        return (1/self.D2O-self.H2O)*((r*27.64)-self.H2O), (1/self.D2O-self.H2O)*((i*27.64)-self.H2O)
    
    def d2o_mol_fraction_head(self):
        return self.d2o_mol_fraction_calc(self.head_solvent.real, self.head_solvent.imag)

    def d2o_mol_fraction_tail(self):
        return self.d2o_mol_fraction_calc(self.tail_solvent.real, self.tail_solvent.imag)

    def wMol_calc(self, r, i):
#         print(self.name,"wMol_calc",r,i)
        return (r * self.D2O) + ((1-r)*self.H2O), (i * self.D2O) + ((1-i)*self.H2O)

    def wMol_head(self):
        d2o_mol_fraction_head_r, d2o_mol_fraction_head_i = self.d2o_mol_fraction_head()
        return self.wMol_calc(d2o_mol_fraction_head_r, d2o_mol_fraction_head_i)

    def wMol_tail(self):
        d2o_mol_fraction_tail_r, d2o_mol_fraction_tail_i = self.d2o_mol_fraction_tail()
        return self.wMol_calc(d2o_mol_fraction_tail_r, d2o_mol_fraction_tail_i)

    def calc_solvent_true(self, r, i):
#         print(self.name,"calc_solvent_true",r,i)
        return complex(r*self.waters_per_head.value/self.vm_head(), i*self.waters_per_head.value/self.vm_head())

    def head_solvent_true(self):
        wMol_r, wMol_i = self.wMol_head()
#         print("wMol_r, wMol_i:", wMol_r, wMol_i)
        return self.calc_solvent_true(wMol_r, wMol_i)

    def tail_solvent_true(self):
        wMol_tail_r, wMol_tail_i = self.wMol_tail()
#         print("wMol_tail_r, wMol_tail_i:", wMol_tail_r, wMol_tail_i)
        return complex(wMol_tail_r*self.waters_per_tail.value/self.vm_tail(), wMol_tail_i*self.waters_per_tail.value/self.vm_tail())

    def vm_head(self):
        return self.vm_heads.value + self.water_vm.value * self.waters_per_head.value

    def vm_tail(self):
        return self.vm_tails.value + self.water_vm.value * self.waters_per_tail.value

    def sld_r(self):
        head_sld = float(self.b_heads_real) / float(self.vm_head()) * 1.e6
        tail_sld = float(self.b_tails_real) / float(self.vm_tail()) * 1.e6
        return head_sld, tail_sld

    def sld_i(self):
        head_sld = float(self.b_heads_imag) / float(self.vm_head()) * 1.e6
        tail_sld = float(self.b_tails_imag) / float(self.vm_tail()) * 1.e6
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
        volfrac = self.vm_head() / (self.apm.value *
                                         self.thickness_heads.value)
        layers[0, 4] = 1 - volfrac
        if self.head_solvent is not None:
            # we do the solvation here, not in Structure.slabs
            layers[0] = Structure.overall_sld(layers[0], self.head_solvent_true())
            layers[0, 4] = 0

        # tail region
        volfrac = self.vm_tail() / (self.apm.value *
                                         self.thickness_tails.value)

        layers[1, 4] = 1 - volfrac
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
        # penalise unphysical volume fractions.
        volfrac_h = self.vm_head() / (self.apm.value *
                                    self.thickness_heads.value)

        # tail region
        volfrac_t = self.vm_tail() / (self.apm.value *
                                    self.thickness_tails.value)

        if volfrac_h > 1 or volfrac_t > 1:
            return -np.inf

        return 0

