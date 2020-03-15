class monolayer:
    def __init__(self, contrast):
    #Define all the neutron b's
    bc = 0.6646e-4  #Carbon
    bo = 0.5804e-4  #Oxygen
    bh = -0.3739e-4 #Hydrogen
    bp = 0.513e-4   #Phosphorus
    bn = 0.936e-4   #Nitrogen
    bd = 0.6671e-4  #Deuterium

    #Formulae of the molecule sections...
    CH2 = (1*bc) + (2*bh);
    CH = (1*bc) + (1*bh);
    CH3 = (1*bc) + (3*bh);
    D2O = (2*bd) + (1*bo);
    H2O = (2*bh) + (1*bo);

    # Calculate mole fraction of D2O from the bulk SLD..
    d2o_molfr = (1/D2O-H2O)*(contrast*27.64)-H2O);
    self.wMol = (d2o_molfr * D2O) + ((1-d2o_molfr)*H2O);

    #sum b's of all the different fragments
    self.sum_b_tails = (28*CH2) + (2*CH) + (2*CH3);

    sum_popc_heads = (8*bo) + (1*bp) + (1*bn) + (2*bc) + (4*CH2) + (3*CH3) + (1*CH);
    sum_popg_heads = (10*bo) + (1*bp) + (2*bc) + (4*CH2) + (2*CH) + (2*bh);
    self.sum_b_heads_pt1 = 3*(sum_popc_heads) + (sum_popg_heads)
    #average composition of lipids in bilayer

    #monolayer compositions
    #sum_m_tails = (34*CH2) + (2*CH3);
    #sum_m_heads = (1*bn) + (2*CH3) + (Waters_per_headD * wMol);

    #volumes of each fragment
    vCH2 = 27.7;
    nCH2 = 30;
    vCH3 = 54.6;
    self.volume_tails = (nCH2 * vCH2) + (2 * vCH3);

    self.volume_heads = 331;

#monolayer volumes
#volume_m_tails = vCH2*36;
#volume_m_heads = 54.6;
    
    def APM(self, apm=None):
        if apm == None:
            self.apm = apm
        return self.apm

#Calculate thickness of each fragment
    def TailThick(self):
        return self.volume_tails / self.APM()

    def HeadThick(self):
        return self.volume_heads / self.APM()

    def sum_b_heads_t(self):
        return (1/4)*(self.sum_b_heads_pt1) + (self.Waters_per_head * self.wMol)
#monolayer thicknesses
#TailThickm = volume_m_tails / APMD;
#HeadThickm = volume_m_heads / APMD;


#and thus the SLD
    def Rho_heads(self):
        return self.sum_b_heads_t() / self.volume_heads

    def Rho_tails(self):
        return self.sum_b_tails() / self.volume_tails

#monolayer SLDs
#Rho_m_heads = sum_m_heads / volume_m_heads;
#Rho_m_tails = sum_m_tails / volume_m_tails;



