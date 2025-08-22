
using DrWatson
@quickactivate :RenewalDiffInDiff

using Distributions
using RenewalDiD

id = parse(Int, ARGS[1])
chain = parse(Int, ARGS[2])
npriors = parse(Int, ARGS[3])
mapmaxtime = parse(Int, ARGS[4])
nsamples = parse(Int, ARGS[5])

#= use line below for arguments when running in REPL 
id = 1; chain = 1; npriors = 1000; mapmaxtime = 60; nsamples = 25;
=#

sampleseed = 1000 * id + chain

sim = load(simulationdir("sim$id.jld2"))["sim"]
model = renewaldid(                      
    sim, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2), 1), 
        mu_delayprior=log(5),
        sigma_gammaprior=Exponential(0.2),
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(8, 2),
        tauprior=Normal(0, 0.2),
    );                          
    mu=0.2, kappa=0.5,               
)

analysis = analysisworkflow(
    model; 
    name="analysis$id", 
    chain, 
    npriors, 
    mapmaxtime, 
    nsamples, 
    priorsseed=id, 
    sampleseed,
)
