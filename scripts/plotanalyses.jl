
using DrWatson
@quickactivate :RenewalDiffInDiff

import Random

using CairoMakie 
using RenewalDiD
using RenewalDiD: Plotting
using Turing


# Simulation 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function loadfittedsamples_sim(
    sim, model; 
    thetainterval=7, sampletries=[5, 10, 25, 100, 500, 1000, 5000, 10_000, 50_000]
)
    return _loadfittedvalues_sim(
        sim, model, thetainterval, sampletries, "samples", "samples"
    )
end

function loadfittedpriorsampless_sim(
    sim, model; 
    thetainterval=7, sampletries=[5, 10, 25, 100, 500, 1000, 5000, 10_000, 50_000]
)
    return _loadfittedvalues_sim(
        sim, model, thetainterval, sampletries, "priors", "priorsamples"
    )
end

function _loadfittedvalues_sim(sim, model, thetai, sampletries::AbstractVector, stype, key)
    for s in sort(sampletries; rev=true)
        _fn = "sim$(sim)_model$(model)_chain1_thetainterval$(thetai)_$(stype)$(s).jld2"
        if isfile(datadir("sims", _fn))
            return _loadfittedvalues_sim(sim, model, thetai, s, stype, key)
        end
    end
end

function _loadfittedvalues_sim(sim, model, thetai, sample::Integer, stype, key)
    samples = _loadfittedvalues_sim_samples(sim, model, thetai, sample, stype, key)
end

function _loadfittedvalues_sim_samples(sim, model, thetai, sample::Integer, stype, key)
    _fn = "sim$(sim)_model$(model)_chain1_thetainterval$(thetai)_$(stype)$(sample).jld2"
    samples = load(datadir("sims", _fn))[key] 

    for i in 2:10
        _fn = "sim$(sim)_model$(model)_chain$(i)_thetainterval$(thetai)_$(stype)$(sample).jld2"
        if isfile(datadir("sims", _fn))
            samples = cat(samples, load(datadir("sims", _fn))[key]; dims=3,)
        end
    end 

    return samples
end

## ineffective intervention

sim1a = load(simulationdir("sim1.jld2"))["sim"].sima
sim1apriors = loadfittedpriorsampless_sim(1, 1)
sim1asamples = loadfittedsamples_sim(1, 1)

RenewalDiD.trplot(DataFrame(sim1asamples); ncols=5, nplots=25)

predmodel_sim1a = let 
    psi_beta_1 = 16
    psi_beta_2 = 4
    predmodel = renewaldidpredmodel(                      
        sim1a, 
        g_seir, 
        RenewalDiDPriors( ; 
            alphaprior=Normal(log(2), 1), 
            sigma_gammaprior=Exponential(0.2),
            sigma_thetaprior=Exponential(0.075), 
            psiprior=Beta(psi_beta_1, psi_beta_2),
            tauprior=Normal(0, 0.2),
            delaydistn=Exponential(1 / 0.3),
        );                          
        mu=0.2, 
        kappa=0.5,
        thetainterval=7,
        Ns=nothing
    )
    predmodel
end

predictions_priors = predictedcases(
    Random.default_rng(), predmodel_sim1a, sim1apriors; 
    mu=0.2, kappa=0.5
)
predictionsr0_priors = predictedR_0(predmodel_sim1a, sim1apriors)
predictionquantiles_priors = quantilerenewaldidinfections(
    predmodel_sim1a, predictions_priors, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0_priors = quantilerenewaldidinfections(
    predmodel_sim1a, predictionsr0_priors, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsplot = RenewalDiD.plotmodel(
    predictionquantilesr0_priors, predictionquantiles_priors, sim1a; 
    linewidth=1, interventionlinestyle=(:dot, :dense),
)

predictions = predictedcases(
    Random.default_rng(), predmodel_sim1a, sim1asamples; 
    mu=0.2, kappa=0.5
)
predictionsr0 = predictedR_0(predmodel_sim1a, sim1asamples)
predictionquantiles = quantilerenewaldidinfections(
    predmodel_sim1a, predictions, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0 = quantilerenewaldidinfections(
    predmodel_sim1a, predictionsr0, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
fittedplot = RenewalDiD.plotmodel(
    predictionquantilesr0, predictionquantiles, sim1a; 
    linewidth=1, interventionlinestyle=(:dot, :dense), linkyaxes=false,
)

## effective intervention

sim1b = load(simulationdir("sim1.jld2"))["sim"].simb

sim1bpriors = let 
    samples = load(datadir("sims", "sim1_model2_chain1_thetainterval7_priors10000.jld2"))["priorsamples"] 
    for i in 2:4 
        samples = cat(
            samples,
            load(datadir("sims", "sim1_model2_chain$(i)_thetainterval7_priors10000.jld2"))["priorsamples"];
            dims=3,
        )
    end 
    samples
end

sim1bsamples = let 
    samples = load(datadir("sims", "sim1_model2_chain1_thetainterval7_samples10000.jld2"))["samples"] 
    for i in 2:4 
        samples = cat(
            samples,
            load(datadir("sims", "sim1_model2_chain$(i)_thetainterval7_samples10000.jld2"))["samples"];
            dims=3,
        )
    end 
    samples
end

RenewalDiD.trplot(DataFrame(sim1bsamples); ncols=5, nplots=25)

predmodel_sim1b = let 
    psi_beta_1 = 16
    psi_beta_2 = 4
    predmodel = renewaldidpredmodel(                      
        sim1b, 
        g_seir, 
        RenewalDiDPriors( ; 
            alphaprior=Normal(log(2), 1), 
            sigma_gammaprior=Exponential(0.2),
            sigma_thetaprior=Exponential(0.075), 
            psiprior=Beta(psi_beta_1, psi_beta_2),
            tauprior=Normal(0, 0.2),
            delaydistn=Exponential(1 / 0.3),
        );                          
        mu=0.2, 
        kappa=0.5,
        thetainterval=7,
        Ns=nothing
    )
    predmodel
end

predictions_priors = predictedcases(Random.default_rng(), predmodel_sim1b, sim1bpriors; mu=0.2, kappa=0.5)
predictionsr0_priors = predictedR_0(predmodel_sim1b, sim1bpriors)
predictionquantiles_priors = quantilerenewaldidinfections(
    predmodel_sim1b, predictions_priors, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0_priors = quantilerenewaldidinfections(
    predmodel_sim1b, predictionsr0_priors, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsplot = RenewalDiD.plotmodel(
    predictionquantilesr0_priors, predictionquantiles_priors, sim1b; 
    linewidth=1, interventionlinestyle=(:dot, :dense),
)

predictions = predictedcases(Random.default_rng(), predmodel_sim1b, sim1bsamples; mu=0.2, kappa=0.5)
predictionsr0 = predictedR_0(predmodel_sim1b, sim1bsamples)
predictionquantiles = quantilerenewaldidinfections(
    predmodel_sim1b, predictions, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0 = quantilerenewaldidinfections(
    predmodel_sim1b, predictionsr0, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
fittedplot = RenewalDiD.plotmodel(
    predictionquantilesr0, predictionquantiles, sim1b; 
    linewidth=1, interventionlinestyle=(:dot, :dense), linkyaxes=false,
)




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
    Ns=RenewalDiD._ns(data), plotproportions=true
)



datasamples_3 = let 
    samples = load(datadir("sims", "ukmaskdata3_model1_chain1_thetainterval7_samples5000.jld2"))["samples"] 
    for i in 2:8 
        samples = cat(
            samples,
            load(datadir("sims", "ukmaskdata3_model1_chain$(i)_thetainterval7_samples5000.jld2"))["samples"] ;
            dims=3,
        )
    end 
    samples
end

RenewalDiD.trplot(DataFrame(datasamples_3); ncols=5, nplots=25)

id = 3
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

samples = datasamples_3

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
    linewidth=1, interventionlinewith=0, linkyaxes=false, Ns=RenewalDiD._ns(data), plotproportions=true

)
