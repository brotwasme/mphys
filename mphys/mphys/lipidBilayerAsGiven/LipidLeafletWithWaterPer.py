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
        
        self. waters_per_tail = possibly_create_parameter(
                waters_per_tail,
                name='%s - waters_per_tail' % name)


    def vm_head(self):
        return self.vm_heads.value + self.water_vm.value * self.waters_per_head.value

    def vm_tail(self):
        return self.vm_tails.value + self.water_vm.value * self.waters_per_tail.value


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
        layers[0, 1] = float(self.b_heads_real) / float(self.vm_head()) * 1.e6
        layers[0, 2] = float(self.b_heads_imag) / float(self.vm_head()) * 1.e6

        layers[1, 1] = float(self.b_tails_real) / float(self.vm_tail()) * 1.e6
        layers[1, 2] = float(self.b_tails_imag) / float(self.vm_tail()) * 1.e6

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
            layers[0] = Structure.overall_sld(layers[0], self.head_solvent)
            layers[0, 4] = 0

        # tail region
        volfrac = self.vm_tail() / (self.apm.value *
                                         self.thickness_tails.value)

        layers[1, 4] = 1 - volfrac
        if self.tail_solvent is not None:
            # we do the solvation here, not in Structure.slabs
            layers[1] = Structure.overall_sld(layers[1], self.tail_solvent)
            layers[1, 4] = 0

        if self.reverse_monolayer:
            layers = np.flipud(layers)
            layers[:, 3] = layers[::-1, 3]

        return layers

    @property
    def parameters(self):
        p = Parameters(name=self.name)
        p.extend([self.apm,
                    self.b_heads_real, self.b_heads_imag, self.vm_heads,
                    self.thickness_heads,
                    self.b_tails_real, self.b_tails_imag, self.vm_tails,
                    
                    self.waters_per_head, self.waters_per_tail,
                    
                    self.thickness_tails, self.rough_head_tail,
                    self.rough_preceding_mono])
        
        if self.head_solvent is not None:
            p.append(self.head_solvent.parameters)
        if self.tail_solvent is not None:
            p.append(self.tail_solvent.parameters)
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

