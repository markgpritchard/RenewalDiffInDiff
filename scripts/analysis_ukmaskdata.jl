
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

data = load(datadir("exp_pro", "covidmaskdata$id.jld2"))["data"]
model = renewaldid(                      
    data, 
    g_covid, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(0, 0.5), 
        sigma_gammaprior=Exponential(0.1),
        sigma_thetaprior=Exponential(0.025), 
        psiprior=Beta(10, 10),
        tauprior=Normal(0, 0.2),
        delaydistn=LogNormal(log(5), log(2)),
    );                          
)

analysis = analysisworkflow(
    model; 
    name="covidmaskanalysis$id", 
    chain, 
    npriors, 
    mapmaxtime, 
    nsamples, 
    priorsseed=id, 
    sampleseed,
)
