
using DrWatson
@quickactivate :RenewalDiffInDiff
using CSV, DataFrames, Dates, Pigeons, Turing
include("setupsimulations.jl")

include("loaddata.jl")
include("loaduscoviddata.jl")

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

function fseir(t, mu=0.5, gamma=0.4)
    return mu * gamma * (exp(-gamma * t) - exp(-mu * t)) / (mu - gamma) 
end

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
        "n_rounds" => 1, 
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

W_sim1_0 = generatew_gt(
    fseir, simulation1dataset["cases_counterfactual"], simulation1dataset["Ns"]
)

W_sim1 = generatew_gt(fseir, simulation1dataset["cases"], simulation1dataset["Ns"])

## Discrete-time analyses

timeperiods2 = let 
    timeperiods = ones(Int, 100)
    timeperiods[50:100] .= 2
    timeperiods
end

sim1model0discrete = diffindiffparameters_twodiscretetimes(
    W_sim1_0, 
    simulation1dataset["cases_counterfactual"],
    simulation1dataset["interventions"], 
    timeperiods2,
    simulation1dataset["Ns"];
    psiprior=0.8
)

s1c0configdiscrete = (
    modelname="sim1model0discrete",
    model=sim1model0discrete,
    n_rounds=n_rounds,
    n_chains=8,
    seed=100+id,
)  
sim1chain0dictdiscrete = produce_or_load(pol_fitparameter, s1c0configdiscrete, datadir("sims"))

sim1model1discrete = diffindiffparameters_twodiscretetimes(
    W_sim1, 
    simulation1dataset["cases"],
    simulation1dataset["interventions"], 
    timeperiods2,
    simulation1dataset["Ns"];
    psiprior=0.8,
)

s1c1configdiscrete = (
    modelname="sim1model1discrete",
    model=sim1model1discrete, 
    n_rounds=n_rounds, 
    n_chains=8, 
    seed=110+id,
)  
sim1chain1dictdiscrete = produce_or_load(pol_fitparameter, s1c1configdiscrete, datadir("sims"))

## Analysis 0 
# No effect of interventions 

sim1model0 = diffindiffparameters_splinetimes(
    W_sim1_0, 
    simulation1dataset["cases_counterfactual"],
    simulation1dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation1dataset["Ns"];
    psiprior=0.8
)

s1c0config = @ntuple modelname="sim1model0" model=sim1model0 n_rounds n_chains=8 seed=100+id
sim1chain0dict = produce_or_load(pol_fitparameter, s1c0config, datadir("sims"))

## Analysis 1 
# "Canonical" difference in differences

sim1model1 = diffindiffparameters_splinetimes(
    W_sim1, 
    simulation1dataset["cases"],
    simulation1dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation1dataset["Ns"];
    psiprior=0.8,
)

s1c1config = @ntuple modelname="sim1model1" model=sim1model1 n_rounds n_chains=8 seed=110+id
sim1chain1dict = produce_or_load(pol_fitparameter, s1c1config, datadir("sims"))


## Analysis 2

# Add lag and lead times to explore non-parallel trends 

sim1model2 = diffindiffparameters_splinetimes(
    W_sim1, 
    simulation1dataset["cases"],
    simulation1dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation1dataset["Ns"];
    psiprior=0.3,
    secondaryinterventions=[
        InterventionsMatrix([ nothing, 36 ], 100),
        InterventionsMatrix([ nothing, 64 ], 100)
    ],
)

s1c2config = @ntuple modelname="sim1model2" model=sim1model2 n_rounds n_chains=8 seed=120+id
sim1chain2dict = produce_or_load(pol_fitparameter, s1c2config, datadir("sims"))


## Analysis a1 
# "fit to curve"

sim1modela1 = diffindiffparameters_fittocurve_splinetimes(
    W_sim1, 
    simulation1dataset["cases"],
    simulation1dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation1dataset["Ns"];
    psiprior=0.8,
)

s1ca1config = @ntuple modelname="sim1modela1" model=sim1modela1 n_rounds n_chains=8 seed=10110+id
sim1chaina1dict = produce_or_load(pol_fitparameter, s1ca1config, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 1a 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

W_sim1a = generatew_gt(fseir, simulation1a_dataset["cases"], simulation1a_dataset["Ns"])

## Analysis 1 

sim1amodel1 = diffindiffparameters_splinetimes(
    W_sim1a, 
    simulation1a_dataset["cases"],
    simulation1a_dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation1a_dataset["Ns"];
    psiprior=0.6,
)

s1ac1config = @ntuple modelname="sim1amodel1" model=sim1amodel1 n_rounds n_chains=8 seed=160+id
sim1achain1dict = produce_or_load(pol_fitparameter, s1ac1config, datadir("sims"))


## Analysis 2

# Add lag and lead times to explore non-parallel trends but do not account for known confounder

sim1amodel2 = diffindiffparameters_splinetimes(
    W_sim1a, 
    simulation1a_dataset["cases"],
    simulation1a_dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation1a_dataset["Ns"];
    psiprior=0.6,
    secondaryinterventions=[
        InterventionsMatrix([ nothing, 36 ], 100),
        InterventionsMatrix([ nothing, 64 ], 100),
    ],
)

s1ac2config = @ntuple modelname="sim1amodel2" model=sim1amodel2 n_rounds n_chains=8 seed=170+id
sim1achain2dict = produce_or_load(pol_fitparameter, s1ac2config, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 2 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

W_sim2 = generatew_gt(fseir, simulation2dataset["cases"], simulation2dataset["Ns"])

## Analysis 0

# Do not account for known confounder

sim2model0 = diffindiffparameters_splinetimes(
    W_sim2, 
    simulation2dataset["cases"],
    simulation2dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation2dataset["Ns"];
    psiprior=0.5,
)

s2c0config = @ntuple modelname="sim2model0" model=sim2model0 n_rounds n_chains=8 seed=200+id
sim2chain0dict = produce_or_load(pol_fitparameter, s2c0config, datadir("sims"))


## Analysis 1 

sim2model1 = diffindiffparameters_splinetimes(
    W_sim2, 
    simulation2dataset["cases"],
    simulation2dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation2dataset["Ns"];
    psiprior=0.5,
    secondaryinterventions=InterventionsMatrix([ nothing, nothing, 70 ], 100),
)

s2c1config = @ntuple modelname="sim2model1" model=sim2model1 n_rounds n_chains=8 seed=210+id
sim2chain1dict = produce_or_load(pol_fitparameter, s2c1config, datadir("sims"))


## Analysis 2

# Add lag and lead times to explore non-parallel trends but do not account for known confounder

sim2model2 = diffindiffparameters_splinetimes(
    W_sim2, 
    simulation2dataset["cases"],
    simulation2dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation2dataset["Ns"];
    psiprior=0.5,
    secondaryinterventions=[
        InterventionsMatrix([ nothing, 36, nothing ], 100),
        InterventionsMatrix([ nothing, nothing, 16 ], 100),
        InterventionsMatrix([ nothing, 64, nothing ], 100),
        InterventionsMatrix([ nothing, nothing, 54 ], 100),
    ],
)

s2c2config = @ntuple modelname="sim2model2" model=sim2model2 n_rounds n_chains=8 seed=220+id
sim2chain2dict = produce_or_load(pol_fitparameter, s2c2config, datadir("sims"))


## Analysis 3

# Add lag and lead times to explore non-parallel trends plus the known confounder

sim2model3 = diffindiffparameters_splinetimes(
    W_sim2, 
    simulation2dataset["cases"],
    simulation2dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation2dataset["Ns"];
    psiprior=0.5,
    secondaryinterventions=[
        InterventionsMatrix([ nothing, nothing, 70 ], 100),
        InterventionsMatrix([ nothing, 36, 16 ], 100),
        InterventionsMatrix([ nothing, 64, 54 ], 100)
    ],
)

s2c3config = @ntuple modelname="sim2model3" model=sim2model3 n_rounds n_chains=8 seed=230+id
sim2chain3dict = produce_or_load(pol_fitparameter, s2c3config, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 3 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

W_sim3_0 = generatew_gt(
    fseir, simulation3dataset["cases_counterfactual"], simulation3dataset["Ns"]
)

W_sim3 = generatew_gt(fseir, simulation3dataset["cases"], simulation3dataset["Ns"])


## Analysis 0 
# No effect of interventions 

sim3model0 = diffindiffparameters_splinetimes(
    W_sim3_0, 
    simulation3dataset["cases_counterfactual"],
    simulation3dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation3dataset["Ns"];
    psiprior=0.3,
)

s3c0config = @ntuple modelname="sim3model0" model=sim3model0 n_rounds n_chains=8 seed=300+id
sim3chain0dict = produce_or_load(pol_fitparameter, s3c0config, datadir("sims"))

# Add lag and lead times to explore non-parallel trends 

sim3model0leadlag = diffindiffparameters_splinetimes(
    W_sim3_0, 
    simulation3dataset["cases_counterfactual"],
    simulation3dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation3dataset["Ns"];
    psiprior=0.3,
    secondaryinterventions=[
        InterventionsMatrix([ nothing, 36 ], 100),
        InterventionsMatrix([ nothing, 64 ], 100)
    ],
)

s3c0configleadlag = @ntuple modelname="sim3model0leadlag" model=sim3model0leadlag n_rounds n_chains=8 seed=305+id
sim3chain0dictleadlag = produce_or_load(pol_fitparameter, s3c0configleadlag, datadir("sims"))


## Analysis 1 

sim3model1 = diffindiffparameters_splinetimes(
    W_sim3, 
    simulation3dataset["cases"],
    simulation3dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation3dataset["Ns"];
    psiprior=0.3,
)

s3c1config = @ntuple modelname="sim3model1" model=sim3model1 n_rounds n_chains=8 seed=310+id
sim3chain1dict = produce_or_load(pol_fitparameter, s3c1config, datadir("sims"))


## Analysis 2

# Add lag and lead times to explore non-parallel trends 

sim3model2 = diffindiffparameters_splinetimes(
    W_sim3, 
    simulation3dataset["cases"],
    simulation3dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation3dataset["Ns"];
    psiprior=0.3,
    secondaryinterventions=[
        InterventionsMatrix([ nothing, 36 ], 100),
        InterventionsMatrix([ nothing, 64 ], 100)
    ],
)

s3c2config = @ntuple modelname="sim3model2" model=sim3model2 n_rounds n_chains=8 seed=320+id
sim3chain2dict = produce_or_load(pol_fitparameter, s3c2config, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 4 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

W_sim4_0 = generatew_gt(
    fseir, simulation4dataset["cases_counterfactual"], simulation4dataset["Ns"]
)

W_sim4 = generatew_gt(fseir, simulation4dataset["cases"], simulation4dataset["Ns"])

## Analysis 0 
# No effect of interventions 

sim4model0 = diffindiffparameters_splinetimes(
    W_sim4_0, 
    simulation4dataset["cases_counterfactual"],
    simulation4dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation4dataset["Ns"];
    psiprior=0.5
)

s4c0config = @ntuple modelname="sim4model0" model=sim4model0 n_rounds n_chains=8 seed=400+id
sim4chain0dict = produce_or_load(pol_fitparameter, s4c0config, datadir("sims"))

## Analysis 1 
# Changes over time modelled as discrete steps 

sim4model1 = diffindiffparameters_splinetimes(
    W_sim4, 
    simulation4dataset["cases"],
    simulation4dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation4dataset["Ns"];
    psiprior=0.5,
)

s4c1config = @ntuple modelname="sim4model1" model=sim4model1 n_rounds n_chains=8 seed=410+id
sim4chain1dict = produce_or_load(pol_fitparameter, s4c1config, datadir("sims"))


## Analysis 2

# Add lag and lead times to explore non-parallel trends 

sim4model2 = diffindiffparameters_splinetimes(
    W_sim4, 
    simulation4dataset["cases"],
    simulation4dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation4dataset["Ns"];
    psiprior=0.5,
    secondaryinterventions=[
        InterventionsMatrix([ nothing, 36 ], 100),
        InterventionsMatrix([ nothing, 64 ], 100)
    ],
)

s4c2config = @ntuple modelname="sim4model2" model=sim4model2 n_rounds n_chains=8 seed=420+id
sim4chain2dict = produce_or_load(pol_fitparameter, s4c2config, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Covid data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

datamodelus1 = diffindiffparameters_splinetimes(
    W_uscoviddata, 
    incidence,
    maskday, 
    [ collect(1.0:28:113); [ 123 ] ],
    populations;
    psiprior=0.8,
)

datac1config = @ntuple modelname="datamodelus1" model=datamodelus1 n_rounds n_chains=12 seed=100+id
datachain1dict = produce_or_load(pol_fitparameter, datac1config, datadir("sims"))

