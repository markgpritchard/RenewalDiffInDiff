
using DrWatson
@quickactivate :RenewalDiffInDiff

using CairoMakie
using Distributions
using RenewalDiD
using RenewalDiD.Plotting

id = 5; chain = 1; npriors = 1000; mapmaxtime = 60; nsamples = 25;

sampleseed = 1000 * id + chain

sim = load(simulationdir("sim$id.jld2"))["sim"]
model = renewaldid(                      
    sim, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2), 1), 
        sigma_gammaprior=Exponential(0.2),
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(8, 2),
        tauprior=Normal(0, 0.2),
        delaydistn=Exponential(1 / 0.3),
    );                          
    mu=0.2, kappa=0.5,               
)

d = priorsworkflow(model; chain, name="testprior", npriors, priorsseed=id)
priorsoutputs1 = samplerenewaldidinfections(model, d)
priorsoutputsquintiles1 = quantilerenewaldidinfections(
    priorsoutputs1, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsoutputplot1 = plotmodel(
    priorsoutputsquintiles1, sim;
    linewidth=1, interventionlinestyle=(:dot, :dense), plotproportions=true
)
