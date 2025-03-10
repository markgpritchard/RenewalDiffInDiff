
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
# Simulation 1 "Canonical" difference in differences
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

W_sim1_0 = generatew_gt(
    fseir, simulation1dataset["cases_counterfactual"], simulation1dataset["Ns"]
)
W_sim1 = generatew_gt(fseir, simulation1dataset["cases"], simulation1dataset["Ns"])

## No effect of interventions 

sim1model0 = diffindiffparameters_splinetimes(
    W_sim1_0, 
    simulation1dataset["cases_counterfactual"],
    simulation1dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation1dataset["Ns"];
    psiprior=0.8,
)

s1c0config = @ntuple modelname="sim1model0" model=sim1model0 n_rounds n_chains=8 seed=100+id
sim1chain0dict = produce_or_load(pol_fitparameter, s1c0config, datadir("sims"))

### with lag and lead times 

sim1model0laglead = diffindiffparameters_splinetimes(
    W_sim1_0, 
    simulation1dataset["cases_counterfactual"],
    simulation1dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation1dataset["Ns"];
    psiprior=0.8,
    secondaryinterventions=lagleadinterventionsmatrix(
        simulation1dataset["interventions"], -21:7:21
    ),
)

s1c0lagleadconfig = (
    modelname="sim1model0laglead",
    model=sim1model0laglead, 
    n_rounds=n_rounds, 
    n_chains=8, 
    seed=105+id,
)
sim1chain0lagleaddict = produce_or_load(pol_fitparameter, s1c0lagleadconfig, datadir("sims"))

## Effective intervention

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

### with lag and lead times 

sim1model1laglead = diffindiffparameters_splinetimes(
    W_sim1, 
    simulation1dataset["cases"],
    simulation1dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation1dataset["Ns"];
    psiprior=0.8,
    secondaryinterventions=lagleadinterventionsmatrix(
        simulation1dataset["interventions"], -21:7:21
    ),
)

s1c1lagleadconfig = (
    modelname="sim1model1laglead",
    model=sim1model1laglead, 
    n_rounds=n_rounds, 
    n_chains=8, 
    seed=115+id,
)
sim1chain1lagleaddict = produce_or_load(pol_fitparameter, s1c1lagleadconfig, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 2: Known confounding variables 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

W_sim2_0 = generatew_gt(
    fseir, simulation2dataset["cases_counterfactual"], simulation2dataset["Ns"]
)
W_sim2 = generatew_gt(fseir, simulation2dataset["cases"], simulation2dataset["Ns"])

## No effect of interventions 

### Do not account for known confounder

sim2model1_0 = diffindiffparameters_splinetimes(
    W_sim2_0, 
    simulation2dataset["cases_counterfactual"],
    simulation2dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation2dataset["Ns"];
    psiprior=0.5,
)

s2c1_0config = @ntuple modelname="sim2model1_0" model=sim2model1_0 n_rounds n_chains=8 seed=250+id
sim2chain1_0dict = produce_or_load(pol_fitparameter, s2c1_0config, datadir("sims"))

### with lag and lead times 

sim2model1_0laglead = diffindiffparameters_splinetimes(
    W_sim2_0, 
    simulation2dataset["cases_counterfactual"],
    simulation2dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation2dataset["Ns"];
    psiprior=0.5,
    secondaryinterventions=lagleadinterventionsmatrix(
        simulation2dataset["interventions"], -21:7:21
    ),
)

s2c1_0lagleadconfig = (
    modelname="sim2model1_0laglead", 
    model=sim2model1_0laglead, 
    n_rounds=n_rounds, 
    n_chains=8, 
    seed=255+id
)
sim2chain1_0lagleaddict = produce_or_load(
    pol_fitparameter, s2c1_0lagleadconfig, datadir("sims")
)

### with the known confounder

sim2model2_0 = diffindiffparameters_splinetimes(
    W_sim2_0, 
    simulation2dataset["cases_counterfactual"],
    simulation2dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation2dataset["Ns"];
    psiprior=0.5,
    secondaryinterventions=InterventionsMatrix([ nothing, nothing, 70 ], 100),
)

s2c2_0config = @ntuple modelname="sim2model2_0" model=sim2model2_0 n_rounds n_chains=8 seed=260+id
sim2chain2_0dict = produce_or_load(pol_fitparameter, s2c2_0config, datadir("sims"))

### with the known confounder, and lag and lead times 

sim2model2_0laglead = diffindiffparameters_splinetimes(
    W_sim2_0, 
    simulation2dataset["cases_counterfactual"],
    simulation2dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation2dataset["Ns"];
    psiprior=0.5,
    secondaryinterventions=[
        [ InterventionsMatrix([ nothing, nothing, 70 ], 100) ];
        lagleadinterventionsmatrix(simulation2dataset["interventions"], -21:7:21)
    ],
)

s2c2_0lagleadconfig = (
    modelname="sim2model2_0laglead", 
    model=sim2model2_0laglead, 
    n_rounds=n_rounds, 
    n_chains=8,
    seed=265+id
)
sim2chain2_0lagleaddict = produce_or_load(
    pol_fitparameter, s2c2_0lagleadconfig, datadir("sims")
)


## With effective intervention

### Do not account for known confounder

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

### with lag and lead times 

sim2model0laglead = diffindiffparameters_splinetimes(
    W_sim2, 
    simulation2dataset["cases"],
    simulation2dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation2dataset["Ns"];
    psiprior=0.5,
    secondaryinterventions=lagleadinterventionsmatrix(
        simulation2dataset["interventions"], -21:7:21
    ),
)

s2c0lagleadconfig = (
    modelname="sim2model0laglead", 
    model=sim2model0laglead,
    n_rounds=n_rounds, 
    n_chains=8, 
    seed=205+id
)
sim2lagleadchain0dict = produce_or_load(pol_fitparameter, s2c0lagleadconfig, datadir("sims"))

### with the known confounder

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

### with the known confounder, and lag and lead times 

sim2model1laglead = diffindiffparameters_splinetimes(
    W_sim2, 
    simulation2dataset["cases"],
    simulation2dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation2dataset["Ns"];
    psiprior=0.5,
    secondaryinterventions=[
        [ InterventionsMatrix([ nothing, nothing, 70 ], 100) ];
        lagleadinterventionsmatrix(simulation2dataset["interventions"], -21:7:21)
    ],)

s2c1lagleadconfig = (
    modelname="sim2model1laglead",
    model=sim2model1laglead,
    n_rounds=n_rounds,
    n_chains=8,
    seed=215+id
)
sim2chain1lagleaddict = produce_or_load(pol_fitparameter, s2c1lagleadconfig, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 3: Non-parallel transmission trends
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

W_sim3_0 = generatew_gt(
    fseir, simulation3dataset["cases_counterfactual"], simulation3dataset["Ns"]
)
W_sim3 = generatew_gt(fseir, simulation3dataset["cases"], simulation3dataset["Ns"])


## No effect of interventions 

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

### with lag and lead times 

sim3model0leadlag = diffindiffparameters_splinetimes(
    W_sim3_0, 
    simulation3dataset["cases_counterfactual"],
    simulation3dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation3dataset["Ns"];
    psiprior=0.3,
    secondaryinterventions=lagleadinterventionsmatrix(
        simulation3dataset["interventions"], -21:7:21
    )
)

s3c0configleadlag = (
    modelname="sim3model0leadlag", 
    model=sim3model0leadlag, 
    n_rounds=n_rounds, 
    n_chains=8, 
    seed=355+id
)
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
# with lag and lead times 

sim3model2 = diffindiffparameters_splinetimes(
    W_sim3, 
    simulation3dataset["cases"],
    simulation3dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation3dataset["Ns"];
    psiprior=0.3,
    secondaryinterventions=lagleadinterventionsmatrix(
        simulation3dataset["interventions"], -21:7:21
    )
)

s3c2config = @ntuple modelname="sim3model2" model=sim3model2 n_rounds n_chains=8 seed=315+id
sim3chain2dict = produce_or_load(pol_fitparameter, s3c2config, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 4 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

W_sim4_0 = generatew_gt(
    fseir, simulation4dataset["cases_counterfactual"], simulation4dataset["Ns"]
)
W_sim4 = generatew_gt(fseir, simulation4dataset["cases"], simulation4dataset["Ns"])

## No effect of interventions 

sim4model0 = diffindiffparameters_splinetimes(
    W_sim4_0, 
    simulation4dataset["cases_counterfactual"],
    simulation4dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation4dataset["Ns"];
    psiprior=0.5,
)

s4c0config = @ntuple modelname="sim4model0" model=sim4model0 n_rounds n_chains=8 seed=400+id
sim4chain0dict = produce_or_load(pol_fitparameter, s4c0config, datadir("sims"))

### with lag and lead times 

sim4model0leadlag = diffindiffparameters_splinetimes(
    W_sim4_0, 
    simulation4dataset["cases_counterfactual"],
    simulation4dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation4dataset["Ns"];
    psiprior=0.5,
    secondaryinterventions=lagleadinterventionsmatrix(
        simulation4dataset["interventions"], -21:7:21
    ),
)

s4c0leadlagconfig = (
    modelname="sim4model0leadlag", 
    model=sim4model0leadlag, 
    n_rounds=n_rounds, 
    n_chains=8,
    seed=405+id
)
sim4chain0leadlagdict = produce_or_load(pol_fitparameter, s4c0leadlagconfig, datadir("sims"))


## Analysis 1 

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
# with lag and lead times 

sim4model2 = diffindiffparameters_splinetimes(
    W_sim4, 
    simulation4dataset["cases"],
    simulation4dataset["interventions"], 
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation4dataset["Ns"];
    psiprior=0.5,
    secondaryinterventions=lagleadinterventionsmatrix(
        simulation4dataset["interventions"], -21:7:21
    ),
)

s4c2config = @ntuple modelname="sim4model2" model=sim4model2 n_rounds n_chains=8 seed=415+id
sim4chain2dict = produce_or_load(pol_fitparameter, s4c2config, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# UK masking data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Convert DataFrame to appropriate matrices 
let 
    # check that each location has the same number of rows 
    @assert sum(maskcoviddf.RegionId .== 1) == sum(maskcoviddf.RegionId .== 2) 
    @assert sum(maskcoviddf.RegionId .== 1) == sum(maskcoviddf.RegionId .== 3) 
    @assert sum(maskcoviddf.RegionId .== 1) == sum(maskcoviddf.RegionId .== 4)
    # how many rows is it?
    covidlength = sum(maskcoviddf.RegionId .== 1)

    # function that is only used here
    function _findfirstendofrestriction(stayathomevec, lowvalue=0)
        calcvec = Vector{Int}(undef, length(stayathomevec))
        for i ∈ eachindex(stayathomevec)
            if i == 1 
                calcvec[i] = 0 
            elseif stayathomevec[i] == lowvalue && stayathomevec[i-1] > lowvalue 
                calcvec[i] = 1
            else 
                calcvec[i] = 0
            end
        end
        return findfirst(x -> x == 1, calcvec)
        # there is a subsequent increase in stay-at-home guidance in England but this
        # applied only to Leicester so is not included as a further intervention here
    end

    global maskcovidcases = Matrix{Int}(undef, covidlength, 4)
    global facialcoveringsrecommended = Matrix{Int}(undef, covidlength, 4)
    global facialcoveringsrequired = Matrix{Int}(undef, covidlength, 4)
    endstayathometimesvec = Vector{Union{Int, Nothing}}(undef, 4)
    somebusinessreopenvec = Vector{Union{Int, Nothing}}(undef, 4)
    global gri_nofc = Matrix{Float64}(undef, covidlength, 4)
    for i ∈ 1:4 
        _tdf = filter(:RegionId => x -> x == i, maskcoviddf)
        for j ∈ 1:covidlength
            if ismissing( _tdf.ConfirmedCases[j])
                maskcovidcases[j, i] = 0
            elseif j == 1 || ismissing( _tdf.ConfirmedCases[j-1])
                maskcovidcases[j, i] = _tdf.ConfirmedCases[j]
            else
                maskcovidcases[j, i] = _tdf.ConfirmedCases[j] - _tdf.ConfirmedCases[j-1]
            end
        end        
        facialcoveringsrecommended[:, i] .= _tdf.FacialCoveringRecommended
        facialcoveringsrequired[:, i] .= _tdf.FacialCoveringRequired
        endstayathometimesvec[i] = _findfirstendofrestriction(_tdf.C6E_Stayathome)
        somebusinessreopenvec[i] = _findfirstendofrestriction(_tdf.C2E_Workplaceclosing, 2)
        gri_nofc[:, i] .= _tdf.Gri_nofc
    end 
    global endstayathometimes = InterventionsMatrix(Int, endstayathometimesvec, covidlength)
    global somebusinessreopen = InterventionsMatrix(Int, somebusinessreopenvec, covidlength)
    # constructed manually to allow leads and lags 
    global facialcoveringsrequired_IM = InterventionsMatrix(
        Int, [ 167, 192, 174, nothing ], 257
    )
end

W_maskcoviddata = generatew_gt(COVIDSERIALINTERVAL, maskcovidcases, POPULATION2020)

## Analysis 2 
# Effect of mask requirements. No other considerations of confounding 

maskingdatamodel2 = diffindiffparameters_splinetimes(
    W_maskcoviddata, 
    maskcovidcases,
    facialcoveringsrequired, 
    [ 1.0; collect(56.0:28:224); 257 ],
    POPULATION2020;
    psiprior=0.4,
)

maskdatac2config = (
    modelname="maskingdatamodel2", 
    model=maskingdatamodel2, 
    n_rounds=n_rounds, 
    n_chains=8, 
    seed=1120+id
)
maskingdatamodel2dict = produce_or_load(pol_fitparameter, maskdatac2config, datadir("sims"))

### with lag and lead times 

maskingdatamodel2leadlag = diffindiffparameters_splinetimes(
    W_maskcoviddata, 
    maskcovidcases,
    facialcoveringsrequired, 
    [ 1.0; collect(56.0:28:224); 257 ],
    POPULATION2020;
    psiprior=0.4,
    secondaryinterventions=lagleadinterventionsmatrix(facialcoveringsrequired_IM, -21:7:21),
)

maskdatac2leadlagconfig = (
    modelname="maskingdatamodel2leadlag", 
    model=maskingdatamodel2leadlag, 
    n_rounds=n_rounds, 
    n_chains=8, 
    seed=1125+id
)
maskingdatamodel2leadlagdict = produce_or_load(
    pol_fitparameter, maskdatac2leadlagconfig, datadir("sims")
)


## Analysis 3 
# Effect of mask requirements with secondary interventions of end of stay-at-home and some
# businesses reopening

maskingdatamodel3 = diffindiffparameters_splinetimes(
    W_maskcoviddata, 
    maskcovidcases,
    facialcoveringsrequired, 
    [ 1.0; collect(56.0:28:224); 257 ],
    POPULATION2020;
    psiprior=0.4,
    secondaryinterventions=[ endstayathometimes, somebusinessreopen ],
)

maskdatac3config = (
    modelname="maskingdatamodel3", 
    model=maskingdatamodel3, 
    n_rounds=n_rounds, 
    n_chains=8, 
    seed=1130+id
)
maskingdatamodel3dict = produce_or_load(pol_fitparameter, maskdatac3config, datadir("sims"))


## Analysis 4 
# with lead and lag times  

maskingdatamodel4 = diffindiffparameters_splinetimes(
    W_maskcoviddata, 
    maskcovidcases, 
    facialcoveringsrequired, 
    [ 1.0; collect(56.0:28:224); 257 ],
    POPULATION2020;
    secondaryinterventions=[ 
        [ endstayathometimes, somebusinessreopen ]; 
        lagleadinterventionsmatrix(facialcoveringsrequired_IM, -21:7:21) 
    ],
)

maskdatac4config = (
    modelname="maskingdatamodel4", 
    model=maskingdatamodel4, 
    n_rounds=n_rounds, 
    n_chains=8, 
    seed=1135+id
)
maskingdatamodel4dict = produce_or_load(pol_fitparameter, maskdatac4config, datadir("sims"))

## Analysis 5 
# Effect of mask requirements with secondary interventions of end of stay-at-home and some
# businesses reopening and mask recommendations

maskingdatamodel5 = diffindiffparameters_splinetimes(
    W_maskcoviddata, 
    maskcovidcases,
    facialcoveringsrequired, 
    [ 1.0; collect(56.0:28:224); 257 ],
    POPULATION2020;
    psiprior=0.4,
    secondaryinterventions=[ 
        endstayathometimes, somebusinessreopen, facialcoveringsrecommended 
    ],
)

maskdatac5config = (
    modelname="maskingdatamodel5", 
    model=maskingdatamodel5, 
    n_rounds=n_rounds, 
    n_chains=8, 
    seed=1150+id
)
maskingdatamodel5dict = produce_or_load(pol_fitparameter, maskdatac5config, datadir("sims"))


## Analysis 6 
# Effect of mask requirements with secondary interventions of end of stay-at-home and some
# businesses reopening and mask recommendations plus lead and lag

maskingdatamodel6 = diffindiffparameters_splinetimes(
    W_maskcoviddata, 
    maskcovidcases,
    facialcoveringsrequired, 
    [ 1.0; collect(56.0:28:224); 257 ],
    POPULATION2020;
    psiprior=0.4,
    secondaryinterventions=[
        [ endstayathometimes, somebusinessreopen, facialcoveringsrecommended ];
        lagleadinterventionsmatrix(facialcoveringsrequired_IM, -21:7:21) 
    ],
)

maskdatac6config = (
    modelname="maskingdatamodel6", 
    model=maskingdatamodel6, 
    n_rounds=n_rounds, 
    n_chains=8, 
    seed=1155+id
)
maskingdatamodel6dict = produce_or_load(pol_fitparameter, maskdatac6config, datadir("sims"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# US data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

const usnrounds = n_rounds > 8 ? 8 : n_rounds

datamodelus1 = diffindiffparameters_splinetimes(
    W_uscoviddata, 
    incidence,
    maskday, 
    [ collect(1.0:28:113); [ 123 ] ],
    populations;
    psiprior=0.8,
)

datac1config = (
    modelname="datamodelus1", 
    model=datamodelus1, 
    n_rounds=usnrounds, 
    n_chains=12,
    seed=101+id
)
datachain1dict = produce_or_load(pol_fitparameter, datac1config, datadir("sims"))

### with lag and lead times 

datamodelus1leadlag = diffindiffparameters_splinetimes(
    W_uscoviddata, 
    incidence,
    maskday, 
    [ collect(1.0:28:113); [ 123 ] ],
    populations;
    psiprior=0.8,
    secondaryinterventions=lagleadinterventionsmatrix(maskday, -21:7:21),
)

datac1leadlagconfig = (
    modelname="datamodelus1leadlag", 
    model=datamodelus1leadlag, 
    n_rounds=usnrounds, 
    n_chains=12,
    seed=105+id
)
datachain1leadlagdict = produce_or_load(
    pol_fitparameter, datac1leadlagconfig, datadir("sims")
)

## With confounding interventions 

datamodelus2 = diffindiffparameters_splinetimes(
    W_uscoviddata, 
    incidence,
    maskday, 
    [ collect(1.0:28:113); [ 123 ] ],
    populations;
    psiprior=0.8,
    secondaryinterventions = [
        relaxshelterinplace, reopenbusiness, reopenrestaurants, reopengyms, reopencinemas 
    ]
)

datac2config = (
    modelname="datamodelus2", 
    model=datamodelus2, 
    n_rounds=usnrounds, 
    n_chains=12,
    seed=110+id
)
datachain2dict = produce_or_load(pol_fitparameter, datac2config, datadir("sims"))

### with confounding interventions, and lag and lead times 

datamodelus2leadlag = diffindiffparameters_splinetimes(
    W_uscoviddata, 
    incidence,
    maskday, 
    [ collect(1.0:28:113); [ 123 ] ],
    populations;
    psiprior=0.8,
    secondaryinterventions = [
        [ relaxshelterinplace, reopenbusiness, reopenrestaurants, reopengyms, reopencinemas ];
        lagleadinterventionsmatrix(maskday, -21:7:21)
    ]
)

datac2leadlagconfig = (
    modelname="datamodelus2leadlag", 
    model=datamodelus2leadlag, 
    n_rounds=usnrounds, 
    n_chains=12,
    seed=115+id
)
datachain2leadlagdict = produce_or_load(
    pol_fitparameter, datac2leadlagconfig, datadir("sims")
)
