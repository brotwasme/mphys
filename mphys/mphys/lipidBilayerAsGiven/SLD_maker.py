import periodictable as pt

class SLD_maker:
    def __init__(self):
        self.POPC_H = {"C":9, "H":12, "O":8, "N":1, "P":1}
        self.POPC_T = {"C":33, "H":70}
        self.POPG_H = {"C":7, "H":6, "O":10, "P":1}
        self.POPG_T = {"C":33, "H":70}

    def __call__(self, type=""):
        if type=="POPC_H":
            self.getPOPC_H()
        
        if type=="POPC_T":
            self.getPOPC_T()
        
        if type=="POPG_H":
            self.getPOPG_H()
        
        if type=="POPG_T":
            self.getPOPG_T()

    def getPOPC_H(self):
        pass
    def getPOPC_T(self):
        pass# {"C":33, "H":70}
    def getPOPG_H(self):
        pass# = {"C":7, "H":6, "O":10, "P":1}
    def getPOPG_T(self):
        pass #{"C":33, "H":70}

    def convert(self, di):
        pass

    def C(self):
        pass
    def H(self):
        pass
    def O(self):
        pass
    def P(self):
        pass

def get_scattering_length(dictionary):
    SL = 0+0j
    for key, value in dictionary.items():
        SL+=pt.elements.symbol(key).neutron.b_c*value
        imag = pt.elements.symbol(key).neutron.b_c_i
        if imag == None:
            imag = 0
        SL+=imag*1j*value
    print(SL)
    return SL
