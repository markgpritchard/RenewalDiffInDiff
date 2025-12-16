
using DrWatson
@quickactivate :RenewalDiffInDiff

import Random

using CairoMakie 
using RenewalDiD
using RenewalDiD: Plotting
using Turing


# View results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## collate chains 

#datapriors_1 = load(datadir("sims", "ukmaskdata1_chain1_thetainterval7_priors10000.jld2"))["priorsamples"]

datapriors_1 = cat(
    load(datadir("sims", "ukmaskdata1_chain1_thetainterval7_priors10000.jld2"))["priorsamples"],
    load(datadir("sims", "ukmaskdata1_chain2_thetainterval7_priors10000.jld2"))["priorsamples"],
    load(datadir("sims", "ukmaskdata1_chain3_thetainterval7_priors10000.jld2"))["priorsamples"],
    load(datadir("sims", "ukmaskdata1_chain4_thetainterval7_priors10000.jld2"))["priorsamples"];
    dims=3,
)

#datasamples_1 = load(datadir("sims", "ukmaskdata1_model1_chain1_thetainterval7_samples5000.jld2"))["samples"] 

datasamples_1 = let 
    samples = load(datadir("sims", "ukmaskdata1_model1_chain1_thetainterval7_samples5000.jld2"))["samples"] 
    for i in 2:8 
        i == 5 && continue
        samples = cat(
            samples,
            load(datadir("sims", "ukmaskdata1_model1_chain$(i)_thetainterval7_samples5000.jld2"))["samples"] ;
            dims=3,
        )
    end 
    samples
end

RenewalDiD.trplot(DataFrame(datasamples_1); ncols=5, nplots=25)

id = 1
data = load(datadir("exp_pro", "covidmaskdata$id.jld2"))["data"]

psi_beta_1 = 1
psi_beta_2 = 9

model = renewaldid(                      
    data, 
    g_covid, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2), 1), 
        sigma_gammaprior=Exponential(0.2),
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(psi_beta_1, psi_beta_2),
        tauprior=Normal(0, 0.2),
        delaydistn=Exponential(1 / 0.3),
    );                          
    thetainterval=7,
    Ns=nothing,
)

predmodel = renewaldidpredmodel(                      
    data, 
    g_covid, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2), 1), 
        sigma_gammaprior=Exponential(0.2),
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(psi_beta_1, psi_beta_2),
        tauprior=Normal(0, 0.2),
        delaydistn=Exponential(1 / 0.3),
    );                          
    thetainterval=7,
    Ns=nothing,
)

priorsamples = datapriors_1
#priorsamples = combined

predictions_priors = predict(Random.default_rng(), predmodel, priorsamples)
predictionsr0_priors = predictedR_0(predmodel, priorsamples)
predictionquantiles_priors = quantilerenewaldidinfections(
    predmodel, predictions_priors, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0_priors = quantilerenewaldidinfections(
    predmodel, predictionsr0_priors, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsplot = plotmodel(
    predictionquantilesr0_priors, predictionquantiles_priors, sim; 
    linewidth=1, interventionlinestyle=(:dot, :dense),
)

samples = datasamples_1

predictions = predictedcases(Random.default_rng(), predmodel, samples)
predictionsr0 = predictedR_0(predmodel, samples)
predictionquantiles = quantilerenewaldidinfections(
    predmodel, predictions, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0 = quantilerenewaldidinfections(
    predmodel, predictionsr0, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
fittedplot = RenewalDiD.plotmodel(
    predictionquantilesr0, predictionquantiles, data; 
    linewidth=1, interventionlinestyle=(:dot, :dense), linkyaxes=false,
)

