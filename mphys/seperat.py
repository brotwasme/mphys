import os.path
import numpy as np
import matplotlib.pyplot as plt
import scipy
import refnx
from refnx.dataset import ReflectDataset, Data1D
from refnx.analysis import Transform, CurveFitter, Objective, Model, Parameter
from refnx.reflect import SLD, Slab, ReflectModel
print('refnx: %s\nscipy: %s\nnumpy: %s' % (refnx.version.version,
                   scipy.version.version, np.version.version))


def seperate(data):#not used
    outs = []

#     np.savetxt("outs.txt", [[0,0]])
    outs.append([0,0])

    air = SLD(0,name="air layer")
    airSlab = air(10,0)

    sio2 = SLD(10,name="bottem layer")
    sio2Slab = sio2(10,0)

    sld1 = SLD(5,name="first layer")
    sld1Slab = sld1(350.)

    sld1Slab.thick.setp(vary=True, bounds=(50,300))
    sld1Slab.sld.real.setp(vary=True, bounds=(4,6))

    structure = airSlab|sld1Slab|sio2Slab

    model = ReflectModel(structure, bkg=3e-6, dq=5.0)
    objective = Objective(model, data, transform=Transform('logY'))
    fitter = CurveFitter(objective)
    fitter.fit("differential_evolution")
    print(objective.parameters)
    lp = objective.logpost()
    print("log: ", lp)

#     np.savetxt("outs.txt", [[1,lp]])
    outs.append([1,lp])

    air = SLD(0,name="air layer")
    airSlab = air(10,0)

    sio2 = SLD(10,name="bottem layer")
    sio2Slab = sio2(10,0)

    sld1 = SLD(5,name="first layer")
    sld1Slab = sld1(350./2,0)

    sld2 = SLD(5,name="second layer")
    sld2Slab = sld2(350./2,0)

    sld1Slab.thick.setp(vary=True, bounds=(50,300))
    sld1Slab.sld.real.setp(vary=True, bounds=(4,6))

    sld2Slab.thick.setp(vary=True, bounds=(50,300))
    sld2Slab.sld.real.setp(vary=True, bounds=(4,6))

    structure = airSlab|sld1Slab|sld2Slab|sio2Slab

    model = ReflectModel(structure, bkg=3e-6, dq=5.0)
    objective = Objective(model, data, transform=Transform('logY'))
    fitter = CurveFitter(objective)
    fitter.fit("differential_evolution")
    print(objective.parameters)
    lp = objective.logpost()
    print("log: ", lp)

#     np.savetxt("outs.txt", [[2,lp]])
    outs.append([2,lp])


    air = SLD(0,name="air layer")
    airSlab = air(10,0)

    sio2 = SLD(10,name="bottem layer")
    sio2Slab = sio2(10,0)

    sld1 = SLD(5,name="first layer")
    sld1Slab = sld1(350./3,0)

    sld2 = SLD(5,name="second layer")
    sld2Slab = sld2(350./3,0)

    sld3 = SLD(5,name="second layer")
    sld3Slab = sld3(350./3,0)

    sld1Slab.thick.setp(vary=True, bounds=(50,250))
    sld1Slab.sld.real.setp(vary=True, bounds=(4,6))

    sld2Slab.thick.setp(vary=True, bounds=(50,250))
    sld2Slab.sld.real.setp(vary=True, bounds=(4,6))

    sld3Slab.thick.setp(vary=True, bounds=(50,250))
    sld3Slab.sld.real.setp(vary=True, bounds=(4,6))

    structure = airSlab|sld1Slab|sld2Slab|sld3Slab|sio2Slab

    model = ReflectModel(structure, bkg=3e-6, dq=5.0)
    objective = Objective(model, data, transform=Transform('logY'))
    fitter = CurveFitter(objective)
    fitter.fit("differential_evolution")
    print(objective.parameters)
    lp = objective.logpost()
    print("log: ", lp)

#     np.savetxt("outs.txt", [[3,lp]])
    outs.append([3,lp])


    air = SLD(0,name="air layer")
    airSlab = air(10,0)

    sio2 = SLD(10,name="bottem layer")
    sio2Slab = sio2(10,0)

    sld1 = SLD(5,name="first layer")
    sld1Slab = sld1(350./4,0)

    sld2 = SLD(5,name="second layer")
    sld2Slab = sld2(350./4,0)

    sld3 = SLD(5,name="second layer")
    sld3Slab = sld3(350./4,0)

    sld4 = SLD(5,name="second layer")
    sld4Slab = sld4(350./4,0)

    sld1Slab.thick.setp(vary=True, bounds=(50,250))
    sld1Slab.sld.real.setp(vary=True, bounds=(4,6))

    sld2Slab.thick.setp(vary=True, bounds=(50,250))
    sld2Slab.sld.real.setp(vary=True, bounds=(4,6))

    sld3Slab.thick.setp(vary=True, bounds=(50,250))
    sld3Slab.sld.real.setp(vary=True, bounds=(4,6))

    sld4Slab.thick.setp(vary=True, bounds=(50,250))
    sld4Slab.sld.real.setp(vary=True, bounds=(4,6))

    structure = airSlab|sld1Slab|sld2Slab|sld3Slab|sld4Slab|sio2Slab

    model = ReflectModel(structure, bkg=3e-6, dq=5.0)
    objective = Objective(model, data, transform=Transform('logY'))
    fitter = CurveFitter(objective)
    fitter.fit("differential_evolution")
    print(objective.parameters)
    lp = objective.logpost()
    print("log: ", lp)

#     np.savetxt("outs.txt", [[4,lp]])
    outs.append([4,lp])
    outs = np.array(outs)
    np.savetxt("outs.txt", outs)
    return outs



def seperateNLayer(data, nLayers, limits = None, doMCMC=False): #used

    air = SLD(0,name="air layer")
    airSlab = air(10,0)

    sio2 = SLD(10,name="bottem layer")
    sio2Slab = sio2(10,0)


    if limits is None:
        limits = [350,50,4,6]
    
#     maxThick = 350
#     lowerThick = 50
#     upperThick = maxThick - nLayers*lowerThick
#     lowerB = 4
#     upperB = 6

    maxThick = float(limits[0])
    lowerThick = limits[1]
    upperThick = maxThick - nLayers*lowerThick
    lowerB = limits[2]
    upperB = limits[3]
    if nLayers>=1:
        sld1 = SLD(5,name="first layer")
        sld1Slab = sld1(maxThick/nLayers,0)

    if nLayers>=2:
        sld2 = SLD(5,name="second layer")
        sld2Slab = sld2(maxThick/nLayers,0)

    if nLayers>=3:
        sld3 = SLD(5,name="second layer")
        sld3Slab = sld3(maxThick/nLayers,0)

    if nLayers>=4:
        sld4 = SLD(5,name="second layer")
        sld4Slab = sld4(maxThick/nLayers,0)

    if nLayers>=1:
        sld1Slab.thick.setp(vary=True, bounds=(lowerThick,upperThick))
        sld1Slab.sld.real.setp(vary=True, bounds=(lowerB,upperB))

    if nLayers>=2:
        sld2Slab.thick.setp(vary=True, bounds=(lowerThick,upperThick))
        sld2Slab.sld.real.setp(vary=True, bounds=(lowerB,upperB))

    if nLayers>=3:
        sld3Slab.thick.setp(vary=True, bounds=(lowerThick,upperThick))
        sld3Slab.sld.real.setp(vary=True, bounds=(lowerB,upperB))

    if nLayers>=4:
        sld4Slab.thick.setp(vary=True, bounds=(lowerThick,upperThick))
        sld4Slab.sld.real.setp(vary=True, bounds=(lowerB,upperB))

    if nLayers>=1:
        structure = airSlab|sld1Slab|sio2Slab
    if nLayers>=2:
        structure = airSlab|sld1Slab|sld2Slab|sio2Slab
    if nLayers>=3:
        structure = airSlab|sld1Slab|sld2Slab|sld3Slab|sio2Slab
    if nLayers>=4:
        structure = airSlab|sld1Slab|sld2Slab|sld3Slab|sld4Slab|sio2Slab

    model = ReflectModel(structure, bkg=3e-6, dq=5.0)
    objective = Objective(model, data, transform=Transform('logY'))
    fitter = CurveFitter(objective)
    if not doMCMC:
        fitter.fit("differential_evolution")
    else:
        fitter.sample(500)
        fitter.sampler.reset()
        fitter.sample(40, nthin=50)
        pass
#     print(objective.parameters)
    lp = objective.logpost()
#     print("log: ", lp)
    return lp

# #     np.savetxt("outs.txt", [[4,lp]])
#     outs.append([4,lp])
#     outs = np.array(outs)
#     np.savetxt("outs.txt", outs)
#     return outs


