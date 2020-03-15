from refnx.reflect import LipidLeaflet
from refnx.analysis import possibly_create_parameter, Parameters, Parameter
# from SLD_maker import get_scattering_length
import lipidBilayerAsGiven.SLD_maker as SLD_maker
get_scattering_length = SLD_maker.get_scattering_length

deuturium = 0.6671e-4
hydrogen = -0.3739e-4
oxygen = 0.5843e-4

POPC_H = {"C":9, "H":12, "O":8, "N":1, "P":1}
POPC_T = {"C":33, "H":70}
POPG_H = {"C":7, "H":6, "O":10, "P":1}
POPG_T = {"C":33, "H":70}

# class Popc2:
#     def __init__(self, solvent_sld, name="POPC"):
#         self.name = name

#     def setupSld(self, sld, name, vary):
#         return possibly_create_parameter(value=sld,
#                 name="{} - {} SLD".format(self.name, name))

class Popc:
    """
    Parameters for the POPC component of the bilayer.

    Parameters
    ----------
    name : str
        Name for the component (default='POPC').
    """
    def __init__(self, solvent_sld, name="POPC"):
        # super(Popc, self).__init__()
        d2o = (2 * 0.6671e-4 + 0.5843e-4) + 0j
        h2o = (2 * -0.3739e-4 + 0.5843e-4) + 0j

        self.s_sld = Parameter(
            value=solvent_sld.real,
            name="{} Solvent SLD".format(name),
            vary=False,
        )
        d2o_molfr = (1 / d2o.real - h2o.real) * ((
            self.s_sld * 1e-6 / 0.036182336306) - h2o.real)
        wmol_real = (d2o_molfr * d2o.real) + (
            (1 - d2o_molfr) * h2o.real)
        wmol_imag = (d2o_molfr * d2o.imag) + (
            (1 - d2o_molfr) * h2o.imag)

        self.water_per_lipid_head = Parameter(
            value=2.8,
            name="{} Waters Per Head".format(name),
            vary=True,
            bounds=(0, 5),
        )
        self.water_per_lipid_tail = Parameter(
            value=0.,
            name="{} Waters Per Tail".format(name),
            vary=False,
        )
        self.b_heads_real = Parameter(
            constraint=get_scattering_length(
                POPC_H).real + self.water_per_lipid_head * wmol_real,
            name="{} Head Scattering Length^-1".format(name),
        )
        self.b_heads_imag = Parameter(
            value=get_scattering_length(
                POPC_H).imag + self.water_per_lipid_head * wmol_imag,
            name="{} Head Imaginary Scattering Length^-1".format(name),
        )
        self.b_tails_real = Parameter(
            value=get_scattering_length(
                POPC_T).real + self.water_per_lipid_tail * wmol_real,
            name="{} Tail Scattering Length^-1".format(name),
            vary=False,
        )
        self.b_tails_imag = Parameter(
            value=get_scattering_length(
                POPC_T).imag + self.water_per_lipid_tail * wmol_imag,
            name="{} Tail Imaginary Scattering Length^-1".format(name),
            vary=False,
        )
        self.vm_heads = Parameter(
            value=331,
            name="{} Volume of Solvated Heads^-3".format(name),
            vary=False,
        )
        self.vm_tails = Parameter(
            value=886.4,
            name="{} Volume of Tails^-3".format(name),
            vary=False,
        )