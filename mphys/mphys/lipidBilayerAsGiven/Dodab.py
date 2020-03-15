import numpy as np
from refnx.analysis import Parameter, Parameters
from refnx.reflect import Component
# import models.helper as helper
import lipidBilayerAsGiven.SLD_maker as SLD_maker
get_scattering_length = SLD_maker.get_scattering_length

DODAB_HEAD = {"C": 2, "H": 6, "N": 1}
DODAB_TAIL = {"C": 18 * 2, "H": 74}


class Dodab(Component):
    """
    Model for DODAB monolayer.

    Parameters
    ----------
    name : str
        Name for the component (default='DODAB').
    """
    def __init__(self, solvent_sld, name="DODAB"):
        super(Dodab, self).__init__()
        d2o = (2 * 0.6671e-4 + 0.5843e-4) + 0j
        h2o = (2 * -0.3739e-4 + 0.5843e-4) + 0j
        
        self.s_sld = Parameter(
            value=solvent_sld.real,
            name="{} Solvent SLD".format(name),
            vary=False,
        )
        d2o_molfr = (1 / d2o.real - h2o.real) * ((1e-6*
            self.s_sld / 0.036182336306) - h2o.real)
        wmol_real = (d2o_molfr * d2o.real) + (
            (1 - d2o_molfr) * h2o.real)
        wmol_imag = (d2o_molfr * d2o.imag) + (
            (1 - d2o_molfr) * h2o.imag)
        
        self.water_per_lipid_head = Parameter(
            value=2.8,
            name="{} Waters Per Head".format(name),
            vary=True,
            bounds=(1, 20),
        )
        self.water_per_lipid_tail = Parameter(
            value=0.,
            name="{} Waters Per Tail".format(name),
            vary=False,
        )
        self.b_h = Parameter(
            constraint=get_scattering_length(
                DODAB_HEAD).real + self.water_per_lipid_head * wmol_real,
            name="{} Head Scattering Length/Ã…".format(name),
        )
        self.bi_h = Parameter(
            constraint=get_scattering_length(
                DODAB_HEAD).imag + self.water_per_lipid_head * wmol_imag,
            name="{} Head Scattering i Length/Ã…".format(name),
        )
        self.b_t = Parameter(
            constraint=get_scattering_length(
                DODAB_TAIL).real + self.water_per_lipid_tail * wmol_real,
            name="{} Tail Scattering Length/Ã…".format(name),
        )
        self.bi_t = Parameter(
            constraint=get_scattering_length(
                DODAB_TAIL).imag + self.water_per_lipid_tail * wmol_imag,
            name="{} Tail Scattering i Length/Ã…".format(name),
        )

        self.apm = Parameter(
            value=48,
            name="{} Area Per Molecule/Ã…^2".format(name),
            vary=True,
            bounds=(30, 120),
        )
        self.roughness = Parameter(
            value=5.8,
            name="{} Monolayer Roughness/Ã…".format(name),
            vary=True,
            bounds=(3, 8),
        )
        self.volume_heads = Parameter(
            value=214,
            name="{} Volume of Solvated Heads/Ã…^3".format(name),
            vary=False,
        )
        self.volume_tails = Parameter(
            value=997.2,
            name="{} Volume of Tails/Ã…^3".format(name),
            vary=False,
        )

    def slabs(self, structure=None):
        """
        Slab representation of Bilayer, as an array

        Parameters
        ----------
        structure : refnx.reflect.Structure
            The Structure hosting this Component
        """
        layers = np.zeros((2, 5))

        layers[1, 0] = self.volume_heads / self.apm
        layers[1, 1] = self.b_h / self.volume_heads #* 1e6
        layers[1, 2] = self.bi_h / self.volume_heads #* 1e6
        layers[1, 3] = self.roughness
        layers[1, 4] = 0

        layers[0, 0] = self.volume_tails / self.apm
        layers[0, 1] = self.b_t / self.volume_tails #* 1e6
        layers[0, 2] = self.bi_t / self.volume_tails #* 1e6
        layers[0, 3] = self.roughness
        layers[0, 4] = 0

        return layers

    @property
    def parameters(self):
        para = Parameters(name=self.name)
        para.extend(
            [
                self.water_per_lipid_head,
                self.water_per_lipid_tail,
                self.b_h,
                self.bi_h,
                self.b_t,
                self.bi_t,
                self.apm,
                self.roughness,
                self.volume_heads,
                self.volume_tails,
            ]
        )
        return para
