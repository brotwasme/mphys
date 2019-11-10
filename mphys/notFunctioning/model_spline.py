from refnx.dataset import ReflectDataset, Data1D
from refnx.analysis import Transform, CurveFitter, Objective, Model, Parameter, process_chain
from refnx.reflect import SLD, Slab, ReflectModel, Spline
from numpy import cumsum, array

def make_model(names, bs, thicks, roughs, fig_i, data, show=False, mcmc=False):
    extent = sum(thicks[:,0]) # (float or Parameter) – Total extent of spline region
    vs = array(bs)[:,0]#(Sequence of float/Parameter) – the real part of the SLD values of each of the knots.
    dz = cum_sum(array(thicks[:,0])) #(Sequence of float/Parameter) – the lateral offset between successive knots.
    print(dz)
    name = "number of nots "+str(len(names))  #(str) – Name of component
    component = Spline(extent, vs, dz, name)
    front = SLD(0)
    front = front(0,0)
    back = SLD(0)
    back = back(0,0)
    structure = front|component|back
    model = ReflectModel(structure, bkg=3e-6, dq=5.0)
    objective = Objective(model, data, transform=Transform('logY'))
    fitter = CurveFitter(objective)
    fitter.fit('differential_evolution');
    return structure, fitter, objective, fig_i+1

