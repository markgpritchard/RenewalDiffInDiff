
using DrWatson
@quickactivate :RenewalDiffInDiff

using RenewalDiD
using Random
#using LineSearches
#using Optim
using Turing
#using PosteriorStats
#using InferenceObjects
using AdvancedHMC

id = parse(Int, ARGS[1])
modeltype = parse(Int, ARGS[2])
chain = parse(Int, ARGS[3])
thetainterval = parse(Int, ARGS[4])
npriors = parse(Int, ARGS[5])
nsamples = parse(Int, ARGS[6])

#= use one of the lines below for arguments when running in REPL 

# for the model with no effective intervention of interest:
id = 1; modeltype = 1; chain = 1; thetainterval = 7; npriors = 1000; nsamples = 1000;

# for the model with an intervention that reduces transmission by 20%
id = 1; modeltype = 2; chain = 1; thetainterval = 7; npriors = 1000; nsamples = 1000;
=#

priorsrng = Xoshiro(id) 
samplesrng = Xoshiro(1000 * id + chain)

if modeltype == 1 
    sim = load(simulationdir("sim$id.jld2"))["sim"].sima
else 
    sim = load(simulationdir("sim$id.jld2"))["sim"].simb
end

filename = "sim$(id)_model$(modeltype)_chain$(chain)_thetainterval$(thetainterval)_samples$(nsamples).jld2"

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
    thetainterval,
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
    thetainterval,
)

priorsamples = sample(priorsrng, model, Prior(), npriors)
#=
if isfile(datadir("sims", filename))
    samples = load(datadir("sims", filename))["samples"]
    samplesrng = load(datadir("sims", filename))["rng"]
else 
    samples = sample(
        samplesrng, model, NUTS(nsamples, 0.65), 1; 
        adtype=AutoMooncake(), initial_params=InitFromPrior(), save_state=true,
    )
    safesave(datadir("sims", filename), Dict("samples" => samples, "rng" => samplesrng))
end
=#
##


samples = sample(
    samplesrng, model, externalsampler(AdvancedHMC.NUTS(0.9; )), nsamples; 
    Ns=nothing, adtype=AutoMooncake(), initial_params=InitFromPrior(), n_adapts=min(1000, nsamples), #save_state=true,
)

#=
using ProfileView



ProfileView.@profview =#
#=
asdf = sample(
        samplesrng, model, Turing.NUTS(25, 0.95), 25; 
        adtype=AutoMooncake(), initial_params=InitFromPrior(), #save_state=true,
    )
    =#
#asdfdf = DataFrame(asdf)    

#=
using AdvancedHMC, Pathfinder

# Running pathfinder
draws = 1_000
result_multi = multipathfinder(model, draws; nruns=4)

# Estimating the metric
inv_metric = result_multi.pathfinder_results[1].fit_distribution.Σ
metric = DenseEuclideanMetric(Matrix(inv_metric))

# Creating an AdvancedHMC NUTS sampler with the custom metric.
n_adapts = 1000 # Number of adaptation steps
tap = 0.9 # Large target acceptance probability to deal with the funnel structure of the posterior
nuts = AdvancedHMC.NUTS(tap; metric=metric)

# Sample
chain = sample(
    samplesrng, model, externalsampler(nuts), 10_000; 
    adtype=AutoMooncake(), initial_params=InitFromPrior(), n_adapts=1_000
)
=#
##
#=

incrementsamples = 10 
remainingsamples = nsamples - size(samples, 1)
while remainingsamples > 0 
    newsamples = sample(
        samplesrng, model, NUTS(nsamples, 0.9), min(incrementsamples, remainingsamples); 
        adtype=AutoMooncake(), 
        initial_params=InitFromPrior(), 
        save_state=true, 
        resume_from=samples,
    )
    samples = vcat(samples, newsamples)
    safesave(datadir("sims", filename), Dict("samples" => samples, "rng" => samplesrng))
    remainingsamples = nsamples - size(samples, 1)
end
=#

#=
filename = "streamed_chain__id$(id)_model$(modeltype)_chain$(chain)_samples$(nsamples).jld2"


samples_7_10 = sample(
    model_7, NUTS(nsamples, 0.65), 10; 
    adtype=AutoMooncake(), initial_params=InitFromPrior(), save_state=true,
)

safesave(datadir("sims", filename), Dict("samples" => samples_7_10))
=#


#=
using CairoMakie 
using RenewalDiD.Plotting=#
#=
predictions_priors_7 = predict(Random.default_rng(), predmodel_7, priorsamples_7)
predictionsr0_priors_7 = predictedR_0(predmodel_7, priorsamples_7)
predictionquantiles_priors_7 = quantilerenewaldidinfections(
    predmodel_7, predictions_priors_7, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0_priors_7 = quantilerenewaldidinfections(
    predmodel_7, predictionsr0_priors_7, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsplot_7 = plotmodel(
    predictionquantilesr0_priors_7, predictionquantiles_priors_7, sim; 
    linewidth=1, interventionlinestyle=(:dot, :dense),
)

predictions_priors_14 = predict(Random.default_rng(), predmodel_14, priorsamples_14)
predictionsr0_priors_14 = predictedR_0(predmodel_14, priorsamples_14)
predictionquantiles_priors_14 = quantilerenewaldidinfections(
    predmodel_14, predictions_priors_14, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0_priors_14 = quantilerenewaldidinfections(
    predmodel_14, predictionsr0_priors_14, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsplot_14 = plotmodel(
    predictionquantilesr0_priors_14, predictionquantiles_priors_14, sim; 
    linewidth=1, interventionlinestyle=(:dot, :dense),
)

predictions_priors_21 = predict(Random.default_rng(), predmodel_21, priorsamples_21)
predictionsr0_priors_21 = predictedR_0(predmodel_21, priorsamples_21)
predictionquantiles_priors_21 = quantilerenewaldidinfections(
    predmodel_21, predictions_priors_21, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0_priors_21 = quantilerenewaldidinfections(
    predmodel_21, predictionsr0_priors_21, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsplot_21 = plotmodel(
    predictionquantilesr0_priors_21, predictionquantiles_priors_21, sim; 
    linewidth=1, interventionlinestyle=(:dot, :dense),
)

predictions_priors_28 = predict(Random.default_rng(), predmodel_28, priorsamples_28)
predictionsr0_priors_28 = predictedR_0(predmodel_28, priorsamples_28)
predictionquantiles_priors_28 = quantilerenewaldidinfections(
    predmodel_28, predictions_priors_28, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0_priors_28 = quantilerenewaldidinfections(
    predmodel_28, predictionsr0_priors_28, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsplot_28 = plotmodel(
    predictionquantilesr0_priors_28, predictionquantiles_priors_28, sim; 
    linewidth=1, interventionlinestyle=(:dot, :dense),
)

predictions = predict(Random.default_rng(), predmodel, samples)
predictionsr0 = predictedR_0(predmodel, samples)
predictionquantiles = quantilerenewaldidinfections(
    predmodel, predictions, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0 = quantilerenewaldidinfections(
    predmodel, predictionsr0, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
fittedplot = plotmodel(
    predictionquantilesr0, predictionquantiles, sim; 
    linewidth=1, interventionlinestyle=(:dot, :dense)
)

predictions_14 = predict(Random.default_rng(), predmodel_14, samples_14)
predictionsr0_14 = predictedR_0(predmodel_14, samples_14)
predictionquantiles_14 = quantilerenewaldidinfections(
    predmodel_14, predictions_14, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0_14 = quantilerenewaldidinfections(
    predmodel_14, predictionsr0_14, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
fittedplot_14 = plotmodel(
    predictionquantilesr0_14, predictionquantiles_14, sim; 
    linewidth=1, interventionlinestyle=(:dot, :dense)
)

predictions_21 = predict(Random.default_rng(), predmodel_21, samples_21)
predictionsr0_21 = predictedR_0(predmodel_21, samples_21)
predictionquantiles_21 = quantilerenewaldidinfections(
    predmodel_21, predictions_21, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0_21 = quantilerenewaldidinfections(
    predmodel_21, predictionsr0_21, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
fittedplot_21 = plotmodel(
    predictionquantilesr0_21, predictionquantiles_21, sim; 
    linewidth=1, interventionlinestyle=(:dot, :dense)
)

predictions_28 = predict(Random.default_rng(), predmodel_28, samples_28)
predictionsr0_28 = predictedR_0(predmodel_28, samples_28)
predictionquantiles_28 = quantilerenewaldidinfections(
    predmodel_28, predictions_28, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
predictionquantilesr0_28 = quantilerenewaldidinfections(
    predmodel_28, predictionsr0_28, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
fittedplot_28 = plotmodel(
    predictionquantilesr0_28, predictionquantiles_28, sim; 
    linewidth=1, interventionlinestyle=(:dot, :dense)
)


psis_loo_val7 = psis_loo(model_7, samples_7)
psis_loo_val14 = psis_loo(model_14, samples_14)
psis_loo_val21 = psis_loo(model_21, samples_21)
psis_loo_val28 = psis_loo(model_28, samples_28)

pointwise_loglikelihoods(model_7, samples_7)

asdf = compare((DataFrame(samples_7), DataFrame(samples_14), DataFrame(samples_21), DataFrame(samples_28), ))

function loo_cv(model, chain; observations=:observedcases)
    loglikelihoods = pointwise_loglikelihoods(model, chain)[observations]
    return loo(loglikelihoods; var_name=observations)

end


asdf = pointwise_loglikelihoods(model_7, samples_7)[:observedcases]
loo

#

using ArviZ


    log_likelihood = pointwise_loglikelihoods(
        model_7, MCMCChains.get_sections(samples_7, :parameters)
    )
    names = sort(collect(keys(log_likelihood)); by=x->parse(Int, split(x, ",")[2][1:end-1]))
    log_likelihood_data = getindex.(Ref(log_likelihood), names)
    log_likelihood = (; data=cat(log_likelihood_data...; dims=3))
end;

idata = from_mcmcchains(samples_7; log_likelihood)


#=
InferenceData with groups:
  > posterior
  > log_likelihood
  > sample_stats
=#
loo(idata) #=
1×9 DataFrame
 Row │ elpd_loo  se       p_loo    n_samples  n_data_points  warning  ⋯
     │ Float64   Float64  Float64  Int64      Int64          Bool     ⋯
─────┼─────────────────────────────────────────────────────────────────
   1 │ -411.174  9.93684  105.107       3000            101    false  ⋯


=#


log_like = PermutedDimsArray(samples_7[:loglikelihood], (:draw, :chain, :school))




#julia> using ArviZExampleData, LogExpFunctions, MCMCDiagnosticTools

#julia> idata = load_example_data("centered_eight");

#julia> log_like = PermutedDimsArray(idata.log_likelihood.obs, (:draw, :chain, :school));

reff = ess(softmax(samples_7[:loglikelihood]; dims=(1, 2)); kind=:basic, split_chains=1, relative=true)

loo(samples_7[:loglikelihood]; reff)

#=
PSISLOOResult with estimates
 elpd  se_elpd    p  se_p
  -31      1.4  0.9  0.33

and PSISResult with 500 draws, 4 chains, and 8 parameters
Pareto shape (k) diagnostic values:
                    Count      Min. ESS
 (-Inf, 0.5]  good  4 (50.0%)  270
  (0.5, 0.7]  okay  4 (50.0%)  307

=#

using AbstractMCMC, Turing, HDF5, MCMCChains

#=
samples_7 = sample(
    model_7, NUTS(nsamples, 0.65), nsamples; 
    adtype=AutoMooncake(), initial_params=InitFromPrior()
)

=#

#id = 1; modeltype = 2; chain = 1; npriors = 1000; nsamples = 100;


nsamples = 1000



chains_reloaded = read("chain-file.jls", Chains)



chains = sample(model, NUTS(), MCMCDistributed(), 1000, 10; save_state = true)




it = AbstractMCMC.steps(model_7, NUTS(nsamples, 0.65); nsamples)  # iterator that yields samples


fid = h5open(h5fname, "w")  # read-write, destroying any existing contents (if any)
d = create_dataset(fid, datadir("sims"), Chains, (nsamples, size(priorsamples_7, 2), 1))



h5open(h5fname, "w") do f
    # prepare a (resizeable) dataset for parameters; choose dtype and initial size
    # here we assume `nparams` parameters in flat vector form
    nparams = size(priorsamples_7, 2)   # replace with the true number of parameters
    d = create_dataset(f, datadir("sims"); dims=(0, nparams), maxdims=(Inf, nparams), chunk=(1000, nparams))
    
    iter = 0
    for sample in it
        iter += 1
        vec = AbstractMCMC.sample_to_vector(sample)  # pseudo helper; extract param vector
        # append: resize then write
        size_old = size(d, 1)
        resize!(d, (size_old+1, nparams))
        d[size_old+1, :] = vec
        if iter % 1000 == 0
        flush(f)   # ensure data written to disk
        end
    end
end

##



nb_p_from_mu_phi(mu, phi) = phi / (phi + mu)

# Smoothly constrain u = exp(logI) into (0, cap), differentiable everywhere.
# - logI: real scalar (log of unconstrained positive latent)
# - cap: positive scalar upper bound (e.g. N[j] * s_prev)
# - k: positive steepness parameter (k=1 gentle, k>1 sharper)
# - eps: tiny stabilizer to avoid division by zero
function smooth_cap_from_log(logI::Real, cap::Real; k::Real=1.0, eps::Real=1e-12)
    @assert cap > 0 "cap must be positive"
    u = exp(logI)
    if k == 1.0
        return (u * cap) / (cap + u + eps)
    else
        up = u^k
        capk = cap^k
        return (up * cap) / (capk + up + eps)
    end
end


function renewaldid2(
    data::AbstractRenewalDiDData, g, priors::RenewalDiDPriors{Q, S, T, U, V, W, X}; 
    observedcases=RenewalDiD.automatic,
    interventions=RenewalDiD.automatic,
    exptdseedcases=RenewalDiD.automatic,
    Ns=RenewalDiD.automatic,
    thetainterval=RenewalDiD.automatic, 
    kwargs...
) where {
    Q <: Distribution, 
    S <: Distribution, 
    T <: Distribution, 
    U <: Distribution, 
    V <: Distribution, 
    W <: Real, 
    X <: Distribution, 
}
    return _renewaldid2(
        RenewalDiD._observedcases(observedcases, data),
        RenewalDiD._interventions(interventions, data),
        RenewalDiD._expectedseedcases(exptdseedcases, data),
        RenewalDiD._ns(Ns, data),
        g,    
        priors.alphaprior,
        priors.psiprior,
        priors.sigma_gammaprior,
        priors.sigma_thetaprior,
        priors.tauprior,
        priors.delaydistn,
        RenewalDiD._nseeds(data),
        priors.omegaprior,
        thetainterval;
        kwargs...
    )
end

function renewaldidpredmodel2(data, g, priors; kwargs...) 
    return renewaldid2(data, g, priors; observedcases=missing, kwargs...)
end


@model function _renewaldid2(
    observedcases,
    interventions,
    expectedseedcases,
    Ns,
    g,    
    alphaprior,
    psiprior,
    sigma_gammaprior,
    sigma_thetaprior,
    tauprior,
    delaydistn,
    n_seeds,
    omega,
    thetainterval;
    maxdelay=RenewalDiD.automatic,
    kwargs...
)
    ngroups = RenewalDiD._ngroups(interventions)
    ntimes = RenewalDiD._ntimes(interventions)
    ninterventions = RenewalDiD._ninterventions(interventions)
    nthetas = RenewalDiD._nthetas(ntimes, thetainterval)

    logtau ~ filldist(tauprior, ninterventions)
    alpha ~ alphaprior
    sigma_gamma ~ sigma_gammaprior
    gammas_raw ~ filldist(Normal(0, 1), ngroups - 1)
    thetas_raw ~ filldist(Normal(0, 1), nthetas)
    sigma_theta ~ sigma_thetaprior
    psi ~ psiprior
   # M_x ~ filldist(Normal(0, 1), ntimes + n_seeds, ngroups)
  #  minsigma2 ~ Beta(1, 2)

    # what's new 
    sigma_infections ~ Exponential(0.5)
    phi_observations ~ truncated(Normal(0, 1), 1e-3, Inf)

    gammavec = RenewalDiD._gammavec(gammas_raw, sigma_gamma)
    thetavec = RenewalDiD._thetavec(thetas_raw, sigma_theta, ntimes; thetainterval)
    predictedlogR_0 = RenewalDiD._predictedlogR_0(alpha, gammavec, thetavec, logtau, interventions)

    #

    expectedinfections = Array{typeof(alpha)}(undef, ntimes + n_seeds, ngroups)
    logpredictedinfections = Array{typeof(alpha)}(undef, ntimes + n_seeds, ngroups)
    predictedinfections = Array{typeof(alpha)}(undef, ntimes + n_seeds, ngroups)

    for j in 1:ngroups 
        s = one(typeof(alpha))
        for t in 1:(n_seeds + ntimes) 
            if t <= n_seeds 
                expectedinfections[t, j] = max(1e-10, expectedseedcases[t, j] * s)
            else
                logrho = alpha + gammavec[j] + thetavec[t-n_seeds] + sum([interventions[t-n_seeds, j, k] * logtau[k] for k in 1:ninterventions])
                expectedinfections[t, j] = max(1e-10, exp(logrho) * sum([predictedinfections[x, j] * g(t - x; kwargs...) for x in 1:(t - 1)]) * s)
            end

            #logpredictedinfections[t, j] = log(expectedinfections[t, j]) + M_x[t, j] * sigma_infections
            logpredictedinfections[t, j] ~ Normal(log(expectedinfections[t, j]), sigma_infections)
            #predictedinfections[t, j] = exp(logpredictedinfections[t, j])
            cap = max(s * Ns[j], 1e-9)  
            predictedinfections[t, j] = smooth_cap_from_log(logpredictedinfections[t, j], cap; k=1.0)
            s = max(0, s - predictedinfections[t, j] / Ns[j])
        end
    end

    #=
    cap = max(s_prev * float(N[j]), 1e-9)     # desired upper bound: N_j * s_{j,t-1}
I[j,t] = smooth_cap_from_log(logI, cap; k=1.0)
=#

    # delay between infection and detection
    delayedinfections = RenewalDiD._delayedinfections(
        typeof(alpha), predictedinfections, delaydistn, ngroups, ntimes, n_seeds, maxdelay
    )


    μ_Z = max.(delayedinfections[n_seeds:(n_seeds + ntimes), :] .* psi, 1e-9)

    p = nb_p_from_mu_phi.(μ_Z, phi_observations)

    if isnan(max(maximum(p), phi_observations))
        @addlogprob! (; loglikelihood=-Inf)
        return nothing
    end

    #@info "size(observedcases) = $(size(observedcases)), size(μ_Z) = $(size(μ_Z))"

    observedcases ~ arraydist(NegativeBinomial.(phi_observations, p))


#=
    if isnan(maximum(sqrtnp_1minusp))
        @addlogprob! (; loglikelihood=-Inf)
        return nothing
    end
=#
    return nothing
end

model = renewaldid2(                      
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
    thetainterval,
)



function renewaldid3(
    data::AbstractRenewalDiDData, g, priors::RenewalDiDPriors{Q, S, T, U, V, W, X}; 
    observedcases=RenewalDiD.automatic,
    interventions=RenewalDiD.automatic,
    exptdseedcases=RenewalDiD.automatic,
    Ns=RenewalDiD.automatic,
    thetainterval=RenewalDiD.automatic, 
    kwargs...
) where {
    Q <: Distribution, 
    S <: Distribution, 
    T <: Distribution, 
    U <: Distribution, 
    V <: Distribution, 
    W <: Real, 
    X <: Distribution, 
}
    return _renewaldid3(
        RenewalDiD._observedcases(observedcases, data),
        RenewalDiD._interventions(interventions, data),
        RenewalDiD._expectedseedcases(exptdseedcases, data),
        RenewalDiD._ns(Ns, data),
        g,    
        priors.alphaprior,
        priors.psiprior,
        priors.sigma_gammaprior,
        priors.sigma_thetaprior,
        priors.tauprior,
        priors.delaydistn,
        RenewalDiD._nseeds(data),
        priors.omegaprior,
        thetainterval;
        kwargs...
    )
end

function renewaldidpredmodel3(data, g, priors; kwargs...) 
    return renewaldid3(data, g, priors; observedcases=missing, kwargs...)
end


@model function _renewaldid3(
    observedcases,
    interventions,
    expectedseedcases,
    Ns,
    g,    
    alphaprior,
    psiprior,
    sigma_gammaprior,
    sigma_thetaprior,
    tauprior,
    delaydistn,
    n_seeds,
    omega,
    thetainterval;
    maxdelay=RenewalDiD.automatic,
    kwargs...
)
    ngroups = RenewalDiD._ngroups(interventions)
    ntimes = RenewalDiD._ntimes(interventions)
    ninterventions = RenewalDiD._ninterventions(interventions)
    nthetas = RenewalDiD._nthetas(ntimes, thetainterval)

    logtau ~ filldist(tauprior, ninterventions)
    alpha ~ alphaprior
    sigma_gamma ~ sigma_gammaprior
    gammas_raw ~ filldist(Normal(0, 1), ngroups - 1)
    thetas_raw ~ filldist(Normal(0, 1), nthetas)
    sigma_theta ~ sigma_thetaprior
    psi ~ psiprior
   # M_x ~ filldist(Normal(0, 1), ntimes + n_seeds, ngroups)
  #  minsigma2 ~ Beta(1, 2)

    # what's new 
    sigma_infections ~ Exponential(0.5)
    phi_observations ~ truncated(Normal(0, 1), 1e-3, Inf)

    gammavec = RenewalDiD._gammavec(gammas_raw, sigma_gamma)
    thetavec = RenewalDiD._thetavec(thetas_raw, sigma_theta, ntimes; thetainterval)
    predictedlogR_0 = RenewalDiD._predictedlogR_0(alpha, gammavec, thetavec, logtau, interventions)

    #

    expectedinfections = Array{typeof(alpha)}(undef, ntimes + n_seeds, ngroups)
    logpredictedinfections = Array{typeof(alpha)}(undef, ntimes + n_seeds, ngroups)
    predictedinfections = Array{typeof(alpha)}(undef, ntimes + n_seeds, ngroups)

    for j in 1:ngroups 
        s = one(typeof(alpha))
        for t in 1:(n_seeds + ntimes) 
            if t <= n_seeds 
                expectedinfections[t, j] = max(1e-10, expectedseedcases[t, j] * s)
            else
                logrho = alpha + gammavec[j] + thetavec[t-n_seeds] + sum([interventions[t-n_seeds, j, k] * logtau[k] for k in 1:ninterventions])
                expectedinfections[t, j] = max(1e-10, exp(logrho) * sum([predictedinfections[x, j] * g(t - x; kwargs...) for x in 1:(t - 1)]) * s)
            end

            #logpredictedinfections[t, j] = log(expectedinfections[t, j]) + M_x[t, j] * sigma_infections
            #logpredictedinfections[t, j] ~ Normal(log(expectedinfections[t, j]), sigma_infections)
            #predictedinfections[t, j] = exp(logpredictedinfections[t, j])
            #cap = max(s * Ns[j], 1e-9)  
            predictedinfections[t, j] = expectedinfections[t, j]
            s = max(0, s - predictedinfections[t, j] / Ns[j])
        end
    end

    #=
    cap = max(s_prev * float(N[j]), 1e-9)     # desired upper bound: N_j * s_{j,t-1}
I[j,t] = smooth_cap_from_log(logI, cap; k=1.0)
=#

    # delay between infection and detection
    delayedinfections = RenewalDiD._delayedinfections(
        typeof(alpha), predictedinfections, delaydistn, ngroups, ntimes, n_seeds, maxdelay
    )


    μ_Z = max.(delayedinfections[n_seeds:(n_seeds + ntimes), :] .* psi, 1e-9)

    p = nb_p_from_mu_phi.(μ_Z, phi_observations)

    if isnan(max(maximum(p), phi_observations))
        @addlogprob! (; loglikelihood=-Inf)
        return nothing
    end

    #@info "size(observedcases) = $(size(observedcases)), size(μ_Z) = $(size(μ_Z))"

    observedcases ~ arraydist(NegativeBinomial.(phi_observations, p))


#=
    if isnan(maximum(sqrtnp_1minusp))
        @addlogprob! (; loglikelihood=-Inf)
        return nothing
    end
=#
    return nothing
end

model = renewaldid3(                      
    sim, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2), 1), 
        sigma_gammaprior=Exponential(0.2),
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(4, 16),
        tauprior=Normal(0, 0.2),
        delaydistn=Exponential(1 / 0.3),
    );                          
    mu=0.2, 
    kappa=0.5,
    thetainterval,
)

import DynamicPPL
function RenewalDiD.predictedR_0(model::DynamicPPL.Model, chain::Chains)
    chaindf = DataFrame(chain)
    return RenewalDiD._modelpredictedR_0(model, chaindf)
end

RenewalDiD._defaults(m::DynamicPPL.Model) = m.defaults
RenewalDiD._delaydistn(m::DynamicPPL.Model) = m.args.delaydistn
RenewalDiD._expectedseedcases(m::DynamicPPL.Model) = m.args.expectedseedcases 
RenewalDiD._generationtimefunction(m::DynamicPPL.Model) = m.args.g
RenewalDiD._interventions(m::DynamicPPL.Model) = m.args.interventions 
RenewalDiD._ngroups(m::DynamicPPL.Model) = RenewalDiD._ngroups(RenewalDiD._interventions(m))
RenewalDiD._ninterventions(m::DynamicPPL.Model) = RenewalDiD._ninterventions(RenewalDiD._interventions(m))
RenewalDiD._ns(m::DynamicPPL.Model) = m.args.Ns 
RenewalDiD._nseeds(m::DynamicPPL.Model) = m.args.n_seeds
RenewalDiD._ntimes(m::DynamicPPL.Model) = RenewalDiD._ntimes(RenewalDiD._interventions(m))

predmodel = renewaldidpredmodel3(                      
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
    thetainterval,
)

priorsamples = sample(priorsrng, model, Prior(), npriors)
=#
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
    predmodel, predictionsr0, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
fittedplot = plotmodel(
    predictionquantilesr0, predictionquantiles, sim; 
    linewidth=1, interventionlinestyle=(:dot, :dense)
)
