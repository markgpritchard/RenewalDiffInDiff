
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
        alphaprior=Normal(log(1.5), 0.5), 
        mu_delayprior=log(5),
        sigma_gammaprior=Exponential(0.1),
        sigma_thetaprior=Exponential(0.025), 
        psiprior=Beta(2, 2),
        tauprior=Normal(0, 0.2),
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


#=
using CairoMakie
using RenewalDiD.Plotting

asdf = priorsworkflow(model; chain, name="testprior", npriors, priorsseed=id)
priorsoutputs1 = samplerenewaldidinfections(
    g_covid, asdf, data
)
priorsoutputsquintiles1 = quantilerenewaldidinfections(
    priorsoutputs1, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsoutputplot1 = plotmodel(
    priorsoutputsquintiles1, data;
    linewidth=1, interventionlinestyle=(:dot, :dense), plotproportions=true
)
=#
