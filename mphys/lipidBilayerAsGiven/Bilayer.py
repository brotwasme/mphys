import numpy as np
from refnx.reflect import Component, SLD, ReflectModel, Structure
from refnx.analysis import possibly_create_parameter, Parameters, Parameter

class Bilayer(Component):
    def __init__(self, Popc, Popg, name="bilayer"):
        super(Bilayer, self).__init__()
        self.name = name
        self.Popc = Popc
        self.Popg = Popg

        self.apm = possibly_create_parameter(
            value = 62,
            name = "{} Area Per Molecule Angstrom^-3".format(name)
            )
        self.apm.setp(vary=True, bounds=(60, 150))

        self.roughness_top = possibly_create_parameter(
            value = 12.0,
            name = "{} roughness top Angstrom^-1".format(name))
        self.roughness_top.setp(vary=True, bounds=(2,20))

        self.roughness_bottom = possibly_create_parameter(
            value = 10.3,
            name = "{} roughness bottom Angstrom^-1".format(name))
        self.roughness_bottom.setp(vary=True, bounds=(4,15))

        self.ratio = possibly_create_parameter(
            value = 0.75,
            name = "{} ratio of Popc to Popg".format(name))
        self.ratio.setp(vary=False)

        self.volFrac = possibly_create_parameter(
            value = 1.,
            name = "{} volume fraction of bilayer".format(name))
        self.volFrac.setp(vary=False)

    def heads_length(self):
        return (self.ratio*self.Popc.vm_heads+ (1-self.ratio)*
                        self.Popg.vm_heads)/self.apm

    def tails_length(self):
        return (self.ratio*self.Popc.vm_tails+ (1-self.ratio)*
                        self.Popg.vm_tails)/self.apm


    def heads_r(self):
        sld_popc_hr = self.Popc.b_heads_real/self.Popc.vm_heads
        sld_popg_hr = self.Popg.b_heads_real/self.Popg.vm_heads
        a=(self.ratio*sld_popc_hr + (1-self.ratio)*
                        sld_popg_hr)/self.apm#*1e6
        print(a)
        return a

    def heads_i(self):
        sld_popc_hi = self.Popc.b_heads_imag/self.Popc.vm_heads
        sld_popg_hi = self.Popg.b_heads_imag/self.Popg.vm_heads
        a = (self.ratio*sld_popc_hi + (1-self.ratio)*
                        sld_popg_hi)/self.apm#*1e6
        print(a)
        return a


    def tail_r(self):
        sld_popc_tr = self.Popc.b_tails_real/self.Popc.vm_tails
        sld_popg_tr = self.Popg.b_tails_real/self.Popg.vm_tails
        return (self.ratio*sld_popc_tr + (1-self.ratio)*
                        sld_popg_tr)/self.apm#*1e6

    def tail_i(self):
        sld_popc_ti = self.Popc.b_tails_imag/self.Popc.vm_tails
        sld_popg_ti = self.Popg.b_tails_imag/self.Popg.vm_tails
        return (self.ratio*sld_popc_ti + (1-self.ratio)*
                        sld_popg_ti)/self.apm#*1e6


    def slabs(self, structure=None):
        layers = np.zeros((4,5))
        layers[0,0] = self.heads_length()
        layers[0,1] = self.heads_r()
        layers[0,2] = self.heads_i()
        layers[0,3] = self.roughness_top
        layers[0,4] = 1-self.volFrac

        layers[1,0] = self.tails_length()
        layers[1,1] = self.tail_r()
        layers[1,2] = self.tail_i()
        layers[1,3] = self.roughness_top
        layers[1,4] = 1-self.volFrac


        layers[2,0] = self.tails_length()
        layers[1,1] = self.tail_r()
        layers[1,2] = self.tail_i()
        layers[2,3] = self.roughness_bottom
        layers[2,4] = 1-self.volFrac
                        
        layers[3,0] = self.heads_length()
        layers[0,1] = self.heads_r()
        layers[0,2] = self.heads_i()
        layers[3,3] = self.roughness_bottom
        layers[3,4] = 1-self.volFrac

        return layers
        
#     def __repr__(self):
#         d = {}
#         d.update(self.__dict__)
#         sld_bh = SLD([self.heads_r(), self.heads_i()])
#         sld_bt = SLD([self.tail_r(), self.tail_i()])
#         d['bh'] = sld_bh
#         d['bt'] = sld_bt

#         s = ("LipidLeaflet({apm!r}, {bh!r}, {vm_heads!r}, {thickness_heads!r},"
#              " {bt!r}, {vm_tails!r}, {thickness_tails!r}, {rough_head_tail!r},"
#              " {rough_preceding_mono!r}, head_solvent={head_solvent!r},"
#              " tail_solvent={tail_solvent!r},"
#              " reverse_monolayer={reverse_monolayer}, name={name!r})")
#         return s.format(**d)

        
    def __repr__(self):
        return str(self.slabs())
    @property
    def parameters(self):
        p = Parameters(name=self.name)
        p.extend([
            self.Popc.s_sld,
            self.Popc.water_per_lipid_head,
            self.Popc.water_per_lipid_tail,
            self.Popc.b_heads_real,
            self.Popc.b_heads_imag,
            self.Popc.b_tails_real,
            self.Popc.b_tails_imag,
            self.Popc.vm_heads,
            self.Popc.vm_tails,

            self.Popg.s_sld,
            self.Popg.water_per_lipid_head,
            self.Popg.water_per_lipid_tail,
            self.Popg.b_heads_real,
            self.Popg.b_heads_imag,
            self.Popg.b_tails_real,
            self.Popg.b_tails_imag,
            self.Popg.vm_heads,
            self.Popg.vm_tails,
            
            self.apm,
            self.roughness_top,
            self.roughness_bottom,
            self.ratio,
            self.volFrac
        ])
        return p
    
#     def logp(self):
#         return 0