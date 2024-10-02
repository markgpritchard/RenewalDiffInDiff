#

using DrWatson
@quickactivate :RenewalDiffInDiff
#include(srcdir("AnalysisFunctions.jl"))
#using .AnalysisFunctions
using CSV, DataFrames, Dates, Pigeons, Turing
include("setupsimulations.jl")
include("loaddata.jl")

testrun = false 

if length(ARGS) == 2 
    id = parse(Int, ARGS[1])
    n_rounds = parse(Int, ARGS[2])
else
    id = 1 
    if testrun 
        n_rounds = 4 
    else
        n_rounds = 10
    end
end

f_seirvector = let 
    u0 = [ 0.0, 1.0, 0.0, 0.0 ]
    p = SEIRParameters(0, 0.5, 0.4, 1.0)
    seirresult = seir_deterministic(u0, p, 1:1000)
    prevalence = seirresult[:, 3]
    prevalence[1:50] ./ sum(prevalence)
end

#offsettimes = [ -42, -21, -7, 7, 21, 42 ]
#offsettimes = 28 .* collect(-6:1:6)

function pol_fitparameter(config)
    @unpack model, modelname, n_rounds, n_chains, seed = config
    roundconfig = @ntuple modelname model n_rounds=1 n_chains seed
    round1 = produce_or_load(pol_fitparameter1, roundconfig, datadir("sims"))
    for n ∈ 2:(n_rounds-1)
        roundconfig = @ntuple modelname model n_rounds=n n_chains seed
        roundn = produce_or_load(pol_fitparameterelement, roundconfig, datadir("sims"))
    end
    oldfilename = "modelname=$(modelname)_n_chains=$(n_chains)_n_rounds=$(n_rounds - 1)_seed=$(seed).jld2"
    pt = load(datadir("sims", oldfilename))["pt"]
    pt = increment_n_rounds!(pt, 1)
    new_pt = pigeons(pt)
    new_chains = Chains(new_pt)

    return Dict(
        "chain" => new_chains, 
        "pt" => new_pt, 
        "modelname" => modelname, 
        "n_rounds" => n_rounds, 
        "n_chains" => n_chains,
    )
end

function pol_fitparameter1(config)
    @unpack model, modelname, n_chains, seed = config
    fitted_pt = pigeons( ;
        target=TuringLogPotential(model),
        n_rounds=1,
        n_chains=n_chains,
        multithreaded=true,
        record=[ traces; record_default() ],
        seed=(seed),
        variational=GaussianReference(),
    )
    fitted_chains = Chains(fitted_pt)
    return Dict(
        "chain" => fitted_chains, 
        "pt" => fitted_pt, 
        "modelname" => modelname, 
        "n_rounds" => n_rounds, 
        "n_chains" => n_chains,
    )
end

function pol_fitparameterelement(config)
    @unpack model, modelname, n_rounds, n_chains, seed = config
    oldfilename = "modelname=$(modelname)_n_chains=$(n_chains)_n_rounds=$(n_rounds - 1)_seed=$(seed).jld2"
    pt = load(datadir("sims", oldfilename))["pt"]
    pt = increment_n_rounds!(pt, 1)
    new_pt = pigeons(pt)
    new_chains = Chains(new_pt)

    return Dict(
        "chain" => new_chains, 
        "pt" => new_pt, 
        "modelname" => modelname, 
        "n_rounds" => n_rounds, 
        "n_chains" => n_chains,
    )
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 1 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Analysis 0 
# No effect of interventions 

W_sim1_0 = generatew_gt(
    f_seirvector, simulation1dataset.counterfactualcases, simulation1dataset.Ns
)

sim1model0 = diffindiffparameters_splinetimes(
    W_sim1_0, 
    simulation1dataset.counterfactualcases,
    simulation1dataset.interventions, 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation1dataset.Ns;
    psiprior=0.8
)

s1c0config = @ntuple modelname="sim1model0" model=sim1model0 n_rounds n_chains=8 seed=10+id
sim1chain0dict = produce_or_load(pol_fitparameter, s1c0config, datadir("sims"))

## Analysis 1 
# "Canonical" difference in differences: 2 discrete time periods 

W_sim1 = generatew_gt(f_seirvector, simulation1dataset.cases, simulation1dataset.Ns)

sim1model1 = diffindiffparameters_splinetimes(
    W_sim1, 
    simulation1dataset.cases,
    simulation1dataset.interventions, 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation1dataset.Ns;
    psiprior=0.8,
)

s1c1config = @ntuple modelname="sim1model1" model=sim1model1 n_rounds n_chains=8 seed=20+id
sim1chain1dict = produce_or_load(pol_fitparameter, s1c1config, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 2 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

W_sim2 = generatew_gt(f_seirvector, simulation2dataset.cases, simulation2dataset.Ns)

## Analysis 1 
# Changes over time modelled as discrete steps 

sim2model1 = diffindiffparameters_splinetimes(
    W_sim1, 
    simulation2dataset.cases,
    simulation2dataset.interventions, 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation2dataset.Ns;
    psiprior=0.5,
)

s2c1config = @ntuple modelname="sim2model1" model=sim2model1 n_rounds n_chains=8 seed=30+id
sim2chain1dict = produce_or_load(pol_fitparameter, s2c1config, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Covid data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Limit to those local authorities near Liverpool that were in the same tiers for control in
# 2020

## Convert DataFrame to appropriate matrices 
let 
    # check that each location has the same number of rows 
    for i ∈ 2:39 
        @assert sum(coviddf.location .== 1) == sum(coviddf.location .== i) 
    end
    
    # how many rows is it?
    covidlength = sum(coviddf.location .== 1)
    #304

    # limit analysis to those local authorities that were in the same tiers 
    # (Halton, Knowsley, Liverpool, Sefton, St Helens, Wirral)
    global allcovidcases = Matrix{Int}(undef, covidlength, 6)
    global pil1covidcases = Matrix{Int}(undef, covidlength, 6)
    #start date in Liverpool 7 November 2020 
    stl = Dates.value(Date("2020-11-07") - Date("2020-05-31"))
    stv = [ i == 3 ? stl : nothing for i ∈ 1:6 ]
    global masstesting = InterventionsMatrix(stv, covidlength)
    for i ∈ 1:6
        k = [ 15, 17, 19,28, 31, 38 ][i]
        _tdf = filter(:location => x -> x == k, coviddf)
        for j ∈ 1:covidlength
            allcovidcases[j, i] = _tdf.cases[j]
            pil1covidcases[j, i] = _tdf.pillar1cases[j]
        end
    end 
end
selectpops = [ populations[x] for x ∈ [ 15, 17, 19, 28, 31, 38 ] ]

W_allcoviddata = generatew_gt(COVIDSERIALINTERVAL, allcovidcases, selectpops)
W_pil1coviddata = generatew_gt(COVIDSERIALINTERVAL, pil1covidcases, selectpops)

## Analysis 1 
# All test results 

datamodel1 = diffindiffparameters_splinetimes(
    W_allcoviddata, 
    allcovidcases,
    masstesting, 
    collect(1.0:28:216),
    selectpops;
    psiprior=0.4,
)

datac1config = @ntuple modelname="datamodel1" model=datamodel1 n_rounds n_chains=8 seed=510+id
datachain1dict = produce_or_load(pol_fitparameter, datac1config, datadir("sims"))

dtimes = let 
    timeperiods = ones(Int, 216)
    for i ∈ 1:216 
        timeperiods[i] = 1 + round(Int, (i - 1) / 43.2, RoundDown)
    end
    timeperiods
end

datamodel1discrete = diffindiffparameters_discretetimes(
    W_allcoviddata, 
    allcovidcases,
    masstesting, 
    dtimes,
    selectpops;
    psiprior=0.4,
)

datac1configdiscrete = (
    modelname="datamodel1discrete",
    model=datamodel1discrete,
    n_rounds=n_rounds,
    n_chains=8,
    seed=515+id
)
datachain1dictdiscrete = produce_or_load(
    pol_fitparameter, datac1configdiscrete, datadir("sims")
)

## Analysis 2 
# Pillar 1 test results 

datamodel2 = diffindiffparameters_splinetimes(
    W_pil1coviddata, 
    pil1covidcases,
    masstesting, 
    collect(1.0:28:216),
    selectpops;
    psiprior=0.4,
)

datac2config = @ntuple modelname="datamodel2" model=datamodel2 n_rounds n_chains=8 seed=520+id
datachain2dict = produce_or_load(pol_fitparameter, datac2config, datadir("sims"))

datamodel2discrete = diffindiffparameters_discretetimes(
    W_pil1coviddata, 
    pil1covidcases,
    masstesting, 
    dtimes,
    selectpops;
    psiprior=0.4,
)

datac2configdiscrete = (
    modelname="datamodel2discrete",
    model=datamodel2discrete,
    n_rounds=n_rounds,
    n_chains=8,
    seed=525+id
)
datachain2dictdiscrete = produce_or_load(
    pol_fitparameter, datac2configdiscrete, datadir("sims")
)
