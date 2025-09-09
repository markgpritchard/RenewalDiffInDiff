
using DrWatson
@quickactivate :RenewalDiffInDiff

using CairoMakie
using Distributions
using RenewalDiD
using RenewalDiD.Plotting

id = 2; chain = 1; npriors = 1000; mapmaxtime = 60; nsamples = 25;

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

d = priorsworkflow(model; chain, name="testprior", npriors, priorsseed=id)
priorsoutputs1 = samplerenewaldidinfections(model, mapoutput[1])#d, indexesformap(d[1]5))
priorsoutputsquintiles1 = quantilerenewaldidinfections(
    priorsoutputs1, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsoutputplot1 = plotmodel(
    priorsoutputsquintiles1, data;
    linewidth=1, interventionlinestyle=(:dot, :dense), plotproportions=true, 
)


#
indexesformap(d[1], 1)


mapoutput = maximumaposterioriworkflow(
        model, d[1]; 
        chain, name="asdf", 
    )