#

using DrWatson
@quickactivate :RenewalDiffInDiff
include(srcdir("AnalysisFunctions.jl"))
using .AnalysisFunctions
using CSV, DataFrames, Dates, Pigeons, Turing
include("simulationtransmissionparameters.jl")
include("loaddata.jl")

testrun = true 

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

# load simulations 
for i ∈ 1:4 
    name = "simulation$(i)dataset"
    var = Symbol(name)
    @eval $var = load(datadir("sims", "$($name).jld2"))
end

f_seirvector = let 
    u0 = SEIRCompartments(0.0, 1.0, 0.0, 0.0)
    p = SEIRParameters(0, 1 / 2, 1 / 2.5, 1.0)
    seirresult = runseir_deterministic(u0, p, 1000)
    prevalence = [ seirresult.compartments[i].I for i ∈ 1:1000 ]
    prevalence[1:50] ./ sum(prevalence)
end

timeperiods2 = let 
    timeperiods = ones(Int, 200)
    timeperiods[100:200] .= 2
    timeperiods
end

timeperiods5 = let 
    timeperiods = ones(Int, 200)
    for i ∈ 1:200 
        timeperiods[i] = 1 + round(Int, (i - 1) / 40, RoundDown)
    end
    timeperiods
end

#offsettimes = [ -42, -21, -7, 7, 21, 42 ]
offsettimes = 28 .* collect(-6:1:6)

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
    f_seirvector, simulation1dataset["cases_counterfactual"], simulation1dataset["Ns"]
)

sim1model0 = diffindiffparameters_discretetimes(
    W_sim1_0, 
    simulation1dataset["cases_counterfactual"],
    simulation1dataset["interventions"], 
    timeperiods2,
    simulation1dataset["Ns"];
    psiprior=0.8
)

s1c0config = @ntuple modelname="sim1model0" model=sim1model0 n_rounds n_chains=8 seed=100+id
sim1chain0dict = produce_or_load(pol_fitparameter, s1c0config, datadir("sims"))

## Analysis 1 
# "Canonical" difference in differences: 2 discrete time periods 

W_sim1 = generatew_gt(f_seirvector, simulation1dataset["cases"], simulation1dataset["Ns"])

sim1model1 = diffindiffparameters_discretetimes(
    W_sim1, 
    simulation1dataset["cases"],
    simulation1dataset["interventions"], 
    timeperiods2,
    simulation1dataset["Ns"];
    psiprior=0.8,
)

s1c1config = @ntuple modelname="sim1model1" model=sim1model1 n_rounds n_chains=8 seed=110+id
sim1chain1dict = produce_or_load(pol_fitparameter, s1c1config, datadir("sims"))

## Analysis 2 
# Assume we do not know it is the canonical situation and instead use five time periods 

sim1model2 = diffindiffparameters_discretetimes(
    W_sim1, 
    simulation1dataset["cases"],
    simulation1dataset["interventions"], 
    timeperiods5,
    simulation1dataset["Ns"];
    psiprior=0.8,
)

s1c2config = @ntuple modelname="sim1model2" model=sim1model2 n_rounds n_chains=8 seed=120+id
sim1chain2dict = produce_or_load(pol_fitparameter, s1c2config, datadir("sims"))

## Analysis 3 
# Add lead and lag times that groups can diverge 

secondaryinterventions_sim1 = [ 
    interventionsoffset(simulation1dataset["interventions"], offsettimes[i]) 
    for i ∈ 1:6
]

sim1model3 = diffindiffparameters_discretetimes(
    W_sim1, 
    simulation1dataset["cases"],
    simulation1dataset["interventions"], 
    timeperiods5,
    simulation1dataset["Ns"];
    psiprior=0.8,
    secondaryinterventions=secondaryinterventions_sim1,
)

s1c3config = @ntuple modelname="sim1model3" model=sim1model3 n_rounds n_chains=8 seed=130+id
sim1chain3dict = produce_or_load(pol_fitparameter, s1c3config, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 2 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

W_sim2 = generatew_gt(f_seirvector, simulation2dataset["cases"], simulation2dataset["Ns"])

## Analysis 1 
# Changes over time modelled as discrete steps 

sim2model1 = diffindiffparameters_discretetimes(
    W_sim2, 
    simulation2dataset["cases"],
    simulation2dataset["interventions"], 
    timeperiods5,
    simulation2dataset["Ns"];
    psiprior=0.45,
)

s2c1config = @ntuple modelname="sim2model1" model=sim2model1 n_rounds n_chains=8 seed=210+id
sim2chain1dict = produce_or_load(pol_fitparameter, s2c1config, datadir("sims"))

## Analysis 2 
# Use splines to model changes over time 

sim2model2 = diffindiffparameters_splinetimes(
    W_sim2, 
    simulation2dataset["cases"],
    simulation2dataset["interventions"], 
    [ [ 1 ]; collect(11:189/4:200) ],
    simulation2dataset["Ns"];
    psiprior=0.45,
)

s2c2config = @ntuple modelname="sim2model2" model=sim2model2 n_rounds n_chains=8 seed=220+id
sim2chain2dict = produce_or_load(pol_fitparameter, s2c2config, datadir("sims"))

## Analysis 3 
# Add lead and lag times that groups can diverge 

secondaryinterventions_sim2 = [ 
    interventionsoffset(simulation2dataset["interventions"], offsettimes[i]) 
    for i ∈ 1:6
]

sim2model3 = diffindiffparameters_splinetimes(
    W_sim2, 
    simulation2dataset["cases"],
    simulation2dataset["interventions"], 
    [ [ 1 ]; collect(11:189/4:200) ],
    simulation2dataset["Ns"];
    psiprior=0.45,
    secondaryinterventions=secondaryinterventions_sim2,
)

s2c3config = @ntuple modelname="sim2model3" model=sim2model3 n_rounds n_chains=8 seed=230+id
sim2chain3dict = produce_or_load(pol_fitparameter, s2c3config, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 3 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Analysis 0 
# No effect of interventions 

W_sim3_0 = generatew_gt(f_seirvector, simulation3dataset["cases_counterfactual"], simulation3dataset["Ns"])

sim3model0 = diffindiffparameters_splinetimes(
    W_sim3_0, 
    simulation3dataset["cases"],
    simulation3dataset["interventions"], 
    [ [ 1 ]; collect(11:189/4:200) ],
    simulation3dataset["Ns"];
    psiprior=0.3,
)

s3c0config = @ntuple modelname="sim3model0" model=sim3model0 n_rounds n_chains=8 seed=300+id
sim3chain0dict = produce_or_load(pol_fitparameter, s3c0config, datadir("sims"))

## Analysis 1 
# Do not account for "confounding" intervention 

W_sim3 = generatew_gt(f_seirvector, simulation3dataset["cases"], simulation3dataset["Ns"])

sim3model1 = diffindiffparameters_splinetimes(
    W_sim3, 
    simulation3dataset["cases"],
    simulation3dataset["interventions"], 
    [ [ 1 ]; collect(11:189/4:200) ],
    simulation3dataset["Ns"];
    psiprior=0.3,
)

s3c1config = @ntuple modelname="sim3model1" model=sim3model1 n_rounds n_chains=8 seed=310+id
sim3chain1dict = produce_or_load(pol_fitparameter, s3c1config, datadir("sims"))

## Analysis 2 
# Add secondary intervention 

sim3model2 = diffindiffparameters_splinetimes(
    W_sim3, 
    simulation3dataset["cases"],
    simulation3dataset["interventions"], 
    [ [ 1 ]; collect(11:189/4:200) ],
    simulation3dataset["Ns"];
    psiprior=0.3,
    secondaryinterventions=simulation3dataset["secondaryinterventions"]
)

s3c2config = @ntuple modelname="sim3model2" model=sim3model2 n_rounds n_chains=8 seed=320+id
sim3chain2dict = produce_or_load(pol_fitparameter, s3c2config, datadir("sims"))

## Analysis 3 
# Add lead and lag times that groups can diverge without explicit secondary intervention

secondaryinterventions_sim3 = [ 
    interventionsoffset(simulation3dataset["interventions"], offsettimes[i]) 
    for i ∈ 1:6
]

sim3model3 = diffindiffparameters_splinetimes(
    W_sim3, 
    simulation3dataset["cases"],
    simulation3dataset["interventions"], 
    [ [ 1 ]; collect(11:189/4:200) ],
    simulation3dataset["Ns"];
    psiprior=0.3,
    secondaryinterventions=secondaryinterventions_sim3,
)

s3c3config = @ntuple modelname="sim3model3" model=sim3model3 n_rounds n_chains=8 seed=330+id
sim3chain3dict = produce_or_load(pol_fitparameter, s3c3config, datadir("sims"))

## Analysis 4 
# Add lead and lag times that groups can diverge and secondary intervention

secondaryinterventions_sim3d = [ 
    secondaryinterventions_sim3; [ simulation3dataset["secondaryinterventions"] ]
]
    
sim3model4 = diffindiffparameters_splinetimes(
    W_sim3, 
    simulation3dataset["cases"],
    simulation3dataset["interventions"], 
    [ [ 1 ]; collect(11:189/4:200) ],
    simulation3dataset["Ns"];
    psiprior=0.3,
    secondaryinterventions=secondaryinterventions_sim3d,
)

s3c4config = @ntuple modelname="sim3model4" model=sim3model4 n_rounds n_chains=8 seed=340+id
sim3chain4dict = produce_or_load(pol_fitparameter, s3c4config, datadir("sims"))

## Analysis 5
# Effect of misspecifying psi 

sim3model5 = diffindiffparameters_splinetimes(
    W_sim3, 
    simulation3dataset["cases"],
    simulation3dataset["interventions"], 
    [ [ 1 ]; collect(11:189/4:200) ],
    simulation3dataset["Ns"];
    psiprior=0.6,
)

s3c5config = @ntuple modelname="sim3model5" model=sim3model5 n_rounds n_chains=8 seed=350+id
sim3chain5dict = produce_or_load(pol_fitparameter, s3c5config, datadir("sims"))

sim3model6 = diffindiffparameters_splinetimes(
    W_sim3, 
    simulation3dataset["cases"],
    simulation3dataset["interventions"], 
    [ [ 1 ]; collect(11:189/4:200) ],
    simulation3dataset["Ns"];
    psiprior=0.15,
)

s3c6config = @ntuple modelname="sim3model6" model=sim3model6 n_rounds n_chains=8 seed=360+id
sim3chain6dict = produce_or_load(pol_fitparameter, s3c6config, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 4 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

W_sim4 = generatew_gt(f_seirvector, simulation4dataset["cases"], simulation4dataset["Ns"])

## Analysis 1 
# No concern about non-parallel trends 

sim4model1 = diffindiffparameters_splinetimes(
    W_sim4, 
    simulation4dataset["cases"],
    simulation4dataset["interventions"], 
    [ [ 1 ]; collect(11:189/4:200) ],
    simulation4dataset["Ns"];
    psiprior=0.45,
)

s4c1config = @ntuple modelname="sim4model1" model=sim4model1 n_rounds n_chains=8 seed=410+id
sim4chain1dict = produce_or_load(pol_fitparameter, s4c1config, datadir("sims"))

## Analysis 2 
# Add lead and lag times that groups can diverge 

secondaryinterventions_sim4 = [ 
    interventionsoffset(simulation4dataset["interventions"], offsettimes[i]) 
    for i ∈ 1:6
]

sim4model2 = diffindiffparameters_splinetimes(
    W_sim4, 
    simulation4dataset["cases"],
    simulation4dataset["interventions"], 
    [ [ 1 ]; collect(11:189/4:200) ],
    simulation4dataset["Ns"];
    psiprior=0.45,
    secondaryinterventions=secondaryinterventions_sim4,
)

s4c2config = @ntuple modelname="sim4model2" model=sim4model2 n_rounds n_chains=8 seed=420+id
sim4chain2dict = produce_or_load(pol_fitparameter, s4c2config, datadir("sims"))

## Analysis 3
# Effect of misspecifying psi 

sim4model3 = diffindiffparameters_splinetimes(
    W_sim4, 
    simulation4dataset["cases"],
    simulation4dataset["interventions"], 
    [ [ 1 ]; collect(11:189/4:200) ],
    simulation4dataset["Ns"];
    psiprior=0.9,
)

s4c3config = @ntuple modelname="sim4model3" model=sim4model3 n_rounds n_chains=8 seed=430+id
sim4chain3dict = produce_or_load(pol_fitparameter, s4c3config, datadir("sims"))


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
