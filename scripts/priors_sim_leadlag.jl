
using DrWatson
@quickactivate :RenewalDiffInDiff

using RenewalDiD
using Random
using Turing
using AdvancedHMC

id = parse(Int, ARGS[1])  # simulation number (1:8)
modeltype = parse(Int, ARGS[2])  # without or with effective intervention
chain = parse(Int, ARGS[3])
thetainterval = parse(Int, ARGS[4])
npriors = parse(Int, ARGS[5])
psi_beta_1 = parse(Int, ARGS[6])
psi_beta_2 = parse(Int, ARGS[7])

#= use one of the lines below for arguments when running in REPL 

# for the model with no effective intervention of interest:
id = 1; modeltype = 1; chain = 1; thetainterval = 7; npriors = 1000; psi_beta_1 = 16; psi_beta_2 = 4;

# for the model with an intervention that reduces transmission by 20%
id = 1; modeltype = 2; chain = 1; thetainterval = 7; npriors = 1000; psi_beta_1 = 16; psi_beta_2 = 4;
=#

@info "running file priors_sim.jl, with parameters id = $id; modeltype = $modeltype; \
    chain = $chain; thetainterval = $thetainterval; npriors = $npriors; \
    psi_beta_1 = $psi_beta_1; psi_beta_2 = $psi_beta_2"

priorsrng = Xoshiro(1000 * id + chain) 

if modeltype == 1 
    sim = load(simulationdir("sim$id.jld2"))["sim"].sima
else 
    sim = load(simulationdir("sim$id.jld2"))["sim"].simb
end

filename = "sim$(id)_model$(modeltype)_chain$(chain)_thetainterval$(thetainterval)_priors$(npriors).jld2"

model = renewaldid(                      
    sim, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2), 1), 
        sigma_gammaprior=Exponential(0.2),
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(psi_beta_1, psi_beta_2),
        tauprior=Normal(0, 0.2),
        delaydistn=Exponential(1 / 0.3),
    );                     
    interventions=InterventionArray(sim.interventions; offset=-21:7:21),     
    mu=0.2, 
    kappa=0.5,
    thetainterval,
    Ns=nothing
)

priorsamples = sample(priorsrng, model, Prior(), npriors)

safesave(
    datadir("sims", filename), 
    Dict("priorsamples" => priorsamples, "rng" => priorsrng)
)
