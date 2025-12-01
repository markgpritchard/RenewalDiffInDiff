
using DrWatson
@quickactivate :RenewalDiffInDiff

using RenewalDiD
using Random
using LineSearches
using Optim
using Turing

id = parse(Int, ARGS[1])
modeltype = parse(Int, ARGS[2])
chain = parse(Int, ARGS[3])
npriors = parse(Int, ARGS[4])
nsamples = parse(Int, ARGS[5])

#= use one of the lines below for arguments when running in REPL 

# for the model with no effective intervention of interest:
id = 1; modeltype = 1; chain = 1; npriors = 1000; nsamples = 100;

# for the model with an intervention that reduces transmission by 20%
id = 1; modeltype = 2; chain = 1; npriors = 1000; nsamples = 100;
=#

sampleseed = 1000 * id + chain

if modeltype == 1 
    sim = load(simulationdir("sim$id.jld2"))["sim"].sima
else 
    sim = load(simulationdir("sim$id.jld2"))["sim"].simb
end

model = renewaldid(                      
    sim, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2), 1), 
        sigma_gammaprior=Exponential(0.2),
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(16, 4),
        tauprior=Normal(0, 0.2),
        delaydistn=Exponential(1 / 0.3),
    );                          
    mu=0.2, 
    kappa=0.5,
    thetainterval=7,
)

priorsamples = sample(model, Prior(), npriors)
#=_, ind = findmax(priorsamples[:lp])

paramvalues = maximum_a_posteriori(
    model, LBFGS(; linesearch=Static()); 
    adtype=AutoMooncake(), maxiters=mapmaxiters, initial_params=Array(priorsamples[ind[1]])[1, :]
)

initial_params = InitFromParams(priorsamples[ind])=#
#=
paramvalues = let 
    vec = Vector{Turing.Optimisation.ModeResult}(undef, 8)

    Threads.@threads for i in 1:8 
        @info "maximum_a_posteriori $i"
        vec[i] = maximum_a_posteriori(
            model, LBFGS(; linesearch=Static()); 
            adtype=AutoMooncake(), maxiters=mapmaxiters,
        )
    end

    vec
end
v, ind = findmax([pv.lp for pv in paramvalues])
initial_params = InitFromParams(paramvalues[ind].params)
=#
samples = sample(
    model, NUTS(nsamples, 0.65), nsamples; 
    adtype=AutoMooncake(), initial_params=InitFromPrior()
)

predmodel = renewaldidpredmodel(                      
    sim, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2), 1), 
        sigma_gammaprior=Exponential(0.2),
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(16, 4),
        tauprior=Normal(0, 0.2),
        delaydistn=Exponential(1 / 0.3),
    );                          
    mu=0.2, 
    kappa=0.5,
    thetainterval=7,
)

using CairoMakie 
using RenewalDiD.Plotting

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

predictions = predict(Random.default_rng(), predmodel, samples)
predictionsr0 = predictedR_0(predmodel, samples)
predictionquantiles = quantilerenewaldidinfections(
    predmodel, predictions, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0 = quantilerenewaldidinfections(
    predmodel, predictions, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
fittedplot = plotmodel(
    predictionquantilesr0, predictionquantiles, sim; 
    linewidth=1, interventionlinestyle=(:dot, :dense)
)
