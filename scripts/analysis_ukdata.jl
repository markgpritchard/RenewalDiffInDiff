
using DrWatson
@quickactivate :RenewalDiffInDiff

using RenewalDiD
using Random
using Turing
using AdvancedHMC

id = parse(Int, ARGS[1])  # simulation number (1:8)
modeltype = parse(Int, ARGS[2])  # not used with data
chain = parse(Int, ARGS[3])
thetainterval = parse(Int, ARGS[4])
nsamples = parse(Int, ARGS[5])
psi_beta_1 = parse(Int, ARGS[6])
psi_beta_2 = parse(Int, ARGS[7])

#= use one of the lines below for arguments when running in REPL 

# for the model with no effective intervention of interest:
id = 1; modeltype = 1; chain = 1; thetainterval = 7; nsamples = 1000; psi_beta_1 = 16; psi_beta_2 = 4;

# for the model with an intervention that reduces transmission by 20%
id = 1; modeltype = 2; chain = 1; thetainterval = 7; nsamples = 1000; psi_beta_1 = 16; psi_beta_2 = 4;
=#

@info "running file analysis_ukdata.jl, with parameters id = $id; \
    chain = $chain; thetainterval = $thetainterval; nsamples = $nsamples; \
    psi_beta_1 = $psi_beta_1; psi_beta_2 = $psi_beta_2"

samplesrng = Xoshiro(2000 * id + chain)

data = load(datadir("exp_pro", "covidmaskdata$id.jld2"))["data"]

filename = "ukmaskdata$(id)_model$(modeltype)_chain$(chain)_thetainterval$(thetainterval)_samples$(nsamples).jld2"

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
    thetainterval,
    Ns=nothing
)

samples = sample(
    samplesrng, model, externalsampler(AdvancedHMC.NUTS(0.9; )), nsamples; 
    adtype=AutoMooncake(), initial_params=InitFromPrior(), n_adapts=min(1000, nsamples),
)

safesave(
    datadir("sims", filename), 
    Dict("samples" => samples, "rng" => samplesrng)
)
