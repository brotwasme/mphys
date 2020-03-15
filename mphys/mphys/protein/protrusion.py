import numpy as np
from refnx.reflect import Slab, Component, SLD, ReflectModel, Structure
from refnx.analysis import possibly_create_parameter, Parameters, Parameter

class Protrusion(Slab):
    def __init__(self, thick, rough, vol_protein, PLratio, solventSLD,  name='', vfsolv=0, interface=None, sld=None):
        if sld == None:
            sld = 0
        super(Protrusion, self).__init__(thick, sld, rough,  name, vfsolv, interface)
        self.vol_protein = possibly_create_parameter(
                vol_protein,
                name='%s - vm_mscl, protien volume' % name)

        self.PLratio = possibly_create_parameter(
                PLratio,
                name='%s - PLratio, protien lipid ratio' % name)
        
        if isinstance(solventSLD, SLD):
            self.solventSLD = solventSLD
        elif isinstance(solventSLD, complex):
            self.solventSLD = SLD(solventSLD)
        else:
            self.solventSLD = possibly_create_parameter(
                solventSLD,
                name='%s - solventSLD' % name)
    
        bo = 0.5804e-4     #Oxygen
        bh = -0.3741e-4    #Hydrogen
        bd = 0.6671e-4     #Deuterium
        D2O = (2*bd) + (1*bo)
        H2O = (2*bh) + (1*bo)
        self.D2O = float(D2O)
        self.H2O = float(H2O)

        bc = 0.6646e-4  #Carbon
        # bo = 0.5804e-4  #Oxygen
        # bh = -0.3739e-4 #Hydrogen
        #bp = 0.513e-4   #Phosphorus
        bn = 0.936e-4   #Nitrogen
        bs = 2.847e-4   #Sulphur
        self.b_protein_part = float(2745*bc + 675*bn + 641*bo + 25*bs + 
                    (4374.5-822.5*0.9)*bh )

    def d2o_mol_fraction_calc(self, r):
        r = r.value*1e-6
#         print(self.name,"d2o_mol_fraction_calc", r, i)
        return (1/self.D2O-self.H2O)*((r*27.64)-self.H2O)#, (1/self.D2O-self.H2O)*((i*27.64)-self.H2O)

    def d2o_mol(self):
        if isinstance(self.solventSLD, SLD) or isinstance(self.solventSLD, complex):
            solventSLD = self.solventSLD.real
        else:
            solventSLD = self.solventSLD.value
            print("solvent sld not an SLD")
        return self.d2o_mol_fraction_calc(solventSLD)

    def b_mscl(self):
        bh = -0.3741e-4 #Hydrogen
        bd = 0.6671e-4  #Deuterium
        mol_frac = self.d2o_mol()
        return float(self.b_protein_part +
                    mol_frac*822.5*0.9*bd + 
                    (1-mol_frac)*822.5*0.9*bh)

    def __repr__(self):
        return (f"Slab({self.thick!r}, {self.sld!r}, {self.rough!r},"
                f" name={self.name!r}, vfsolv={self.vfsolv!r},"
                f" interface={self.interfaces!r},"
                f"plratio={self.PLratio!r},"
                f"solventSLD={self.solventSLD!r},"
                f"vol_protein={self.vol_protein!r})")

    @property
    def parameters(self):
        p = Parameters(name=self.name)
        p.extend([self.thick])
        # p.extend(self.sld.parameters)
        p.extend([self.rough, self.vol_protein, self.vfsolv])
        # p.name = self.name
        # p.extend([self.vol_protein,
        # self.PLratio])#,self.solventSLD
        # if isinstance(solventSLD, SLD):
        #     p.extend([self.solventSLD.parameters])
        # elif isinstance(solventSLD, complex):
        #     self.solventSLD = SLD(solventSLD)
        # else:
        #     p.extend([self.solventSLD
        return p
    
    def slabs(self, structure=None):
        layers = super(Protrusion, self).slabs(structure)
        proteinSLD = self.b_mscl()/self.vol_protein.value
        layers[0, 1] = (1-self.PLratio.value)*proteinSLD + self.PLratio.value*proteinSLD
        layers[0, 2] = 0
        # print("layers",layers)
        return layers