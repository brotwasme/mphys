import numpy as np
from refnx.reflect import Component, SLD, ReflectModel, Structure, LipidLeaflet
from refnx.analysis import possibly_create_parameter, Parameters, Parameter

class LipidLeafletWithProtien(LipidLeaflet):
    def __init__(self, apm, b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono,

                 water_vm, waters_per_head, waters_per_tail, b_mscl, vm_mscl, PLRatio,
                
                 head_solvent=None, tail_solvent=None,
                 reverse_monolayer=False, name=''):

        super(LipidLeafletWithProtien, self).__init__(apm, b_heads, vm_heads, thickness_heads,
                 b_tails, vm_tails, thickness_tails, rough_head_tail,
                 rough_preceding_mono, head_solvent, tail_solvent,
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

        self.water_vm = possibly_create_parameter(
                water_vm,
                name='%s - water_vm' % name)

        self.waters_per_head = possibly_create_parameter(
                waters_per_head,
                name='%s - waters_per_head' % name)
        
        self. waters_per_tail = possibly_create_parameter(
                waters_per_tail,
                name='%s - waters_per_tail' % name)


    def sld_(self, b_heads, b_tails, b_mscl):
        lipid_head_sld = float(b_heads) / self.vm_head() * 1.e6
        lipid_tail_sld = float(b_tails) / self.vm_tail() * 1.e6

        mscl_sld = float(b_mscl)/float(self.vm_mscl.value) * 1.e6

        mscl_head_sld = mscl_sld*(float(self.thickness_heads)/
                                  (2*(float(self.thickness_heads)+float(self.thickness_tails))))
        mscl_tail_sld = mscl_sld*(float(self.thickness_heads)/
                                  (2*(float(self.thickness_heads)+float(self.thickness_tails))))

        head_sld = (self.PLRatio*lipid_head_sld) + (1-self.PLRatio)*mscl_head_sld
        tail_sld = (self.PLRatio*lipid_tail_sld) + (1-self.PLRatio)*mscl_tail_sld
        return head_sld, tail_sld

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
        head_sld_real, tail_sld_real = self.sld_(self.b_heads_real, #real
                                            self.b_tails_real,
                                            self.b_mscl_real)
        head_sld_imag, tail_sld_imag = self.sld_(self.b_heads_imag, #imaginary
                                            self.b_tails_imag,
                                            self.b_mscl_imag)
        layers[0, 1] = head_sld_real
        layers[0, 2] = head_sld_imag

        layers[1, 1] = tail_sld_real
        layers[1, 2] = tail_sld_imag

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
                    self.b_mscl_real, self.b_mscl_imag, self.vm_mscl,
                    
                    self.thickness_tails, self.rough_head_tail,
                    self.rough_preceding_mono,
                    
                    self.PLRatio])
        
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

