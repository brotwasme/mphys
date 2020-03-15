from refnx.dataset import ReflectDataset, Data1D
from refnx.analysis import Transform, CurveFitter, Objective, Model, Parameter, process_chain
from refnx.reflect import SLD, Slab, ReflectModel

def make_model(names, bs, thicks, roughs, fig_i, data, show=False, mcmc=False):
    no_layers = len(bs)
    layers = []
    for i in range(no_layers):
        names.append('layer'+str(i))

    for i in range(no_layers):
            sld = SLD(bs[i][0],name=names[i])
            layers.append(sld(thicks[i][0],roughs[i][0]))

    layers[0].thick.setp(vary=True, bounds=(thicks[i][1],thicks[i][2]))
    layers[0].sld.real.setp(vary=True, bounds=(bs[i][1],bs[i][2]))

    for layer in layers[1:]:
        layer.thick.setp(vary=True, bounds=(thicks[i][1],thicks[i][2]))
        layer.sld.real.setp(vary=True, bounds=(bs[i][1],bs[i][2]))
        layer.rough.setp(vary=True, bounds=(roughs[i][1],roughs[i][2]))

    structure = layers[0]
    for layer in layers[1:]:
        structure|=layer
    print(structure)
    model = ReflectModel(structure, bkg=3e-6, dq=5.0)
    #model.scale.setp(bounds=(0.6, 1.2), vary=True)
    #model.bkg.setp(bounds=(1e-9, 9e-6), vary=True)
    objective = Objective(model, data, transform=Transform('logY'))
    fitter = CurveFitter(objective)
    if mcmc:
        fitter.sample(1000);
        process_chain(objective, fitter.chain, nburn=300, nthin=100)
    else:
        fitter.fit('differential_evolution');
    print(objective.parameters)

    if show:
        plt.figure(fig_i)
        plt.plot(*structure.sld_profile())
        plt.ylabel('SLD /$10^{-6} \AA^{-2}$')
        plt.xlabel('distance / $\AA$');
    return structure, fitter, objective, fig_i+1


if __name__ == "__main__":
    import numpy as np
    max_thickness = 350
    structs = []
    fitrs = []
    ln_posts = []
    fig_i = 0
    import data_in as di
    data = di.data_in("29553_54.dat")
    q, R, sim_dR = data[0],data[1],data[2]
    #print(q, R, sim_dR)
    from refnx.dataset import Data1D
    data = Data1D(data=(q, R, sim_dR))
    import make_model2 as mm
    for i in range(3,4):
        thick = round(max_thickness/(i))
        names = []
        bs = []
        thicks = []
        roughs = []
        for j in range(i):
            names.append('layers'+str(j))
            bs.append([5, 4, 6])
            thicks.append([thick, 100, 150])
            roughs.append([0, 0, 5])
        print(names, bs, thicks, roughs)
        structure, fitter, objective, fig_i = mm.make_model(names,
                                   bs, thicks, roughs, fig_i, data)
        ln_post = objective.logpost()
        print("log post out: ", ln_post, "\npost out: ", np.exp(ln_post))
        structs.append(structure)
        fitrs.append(fitter)
        ln_posts.append(ln_post)