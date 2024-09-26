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
        n_rounds = 12
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
#=
timeperiods5 = let 
    timeperiods = ones(Int, 200)
    for i ∈ 1:200 
        timeperiods[i] = 1 + round(Int, (i - 1) / 40, RoundDown)
    end
    timeperiods
end
=#
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
# Functions for analyses 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 1 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

W_sim1 = generatew_gt(f_seirvector, simulation1dataset["cases"], simulation1dataset["Ns"])

## Analysis 0 
# "Canonical" difference in differences: 2 discrete time periods 

sim1model0 = diffindiffparameters_discretetimes(
    W_sim1, 
    simulation1dataset["cases"],
    simulation1dataset["interventions"], 
    timeperiods2,
    simulation1dataset["Ns"];
    psiprior=0.8,
)
s1c0config = @ntuple modelname="sim1model0" model=sim1model0 n_rounds n_chains=8 seed=100+id
sim1chain0dict = produce_or_load(pol_fitparameter, s1c0config, datadir("sims"))

## Analysis 1
# With spline for time-varying parameter 

sim1model1 = diffindiffparameters_splinetimes(
    W_sim1, 
    simulation1dataset["cases"],
    simulation1dataset["interventions"], 
    [ [ 1 ]; collect(11:189/4:200) ],
    simulation1dataset["Ns"];
    psiprior=0.8,
)
s1c1config = @ntuple modelname="sim1model1" model=sim1model1 n_rounds n_chains=8 seed=110+id
sim1chain1dict = produce_or_load(pol_fitparameter, s1c1config, datadir("sims"))

## Analysis 2
# With secondary intervention.
# This simulation does not have a secondary intervention, but we state that group has one on
# day 50 with an effect size of 0.



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

## Convert DataFrame to appropriate matrices 
let 
    # check that each location has the same number of rows 
    @assert sum(coviddf.RegionId .== 1) == sum(coviddf.RegionId .== 2) 
    @assert sum(coviddf.RegionId .== 1) == sum(coviddf.RegionId .== 3) 
    @assert sum(coviddf.RegionId .== 1) == sum(coviddf.RegionId .== 4)
    # how many rows is it?
    covidlength = sum(coviddf.RegionId .== 1)

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

    global covidcases = Matrix{Int}(undef, covidlength, 4)
    global facialcoveringsrecommended = Matrix{Int}(undef, covidlength, 4)
    global facialcoveringsrequired = Matrix{Int}(undef, covidlength, 4)
    endstayathometimesvec = Vector{Union{Int, Nothing}}(undef, 4)
    somebusinessreopenvec = Vector{Union{Int, Nothing}}(undef, 4)
    global gri_nofc = Matrix{Float64}(undef, covidlength, 4)
    for i ∈ 1:4 
        _tdf = filter(:RegionId => x -> x == i, coviddf)
        for j ∈ 1:covidlength
            if j == 1 
                if ismissing(_tdf.ConfirmedCases[j])
                    covidcases[j, i] = 0
                else
                    covidcases[j, i] = _tdf.ConfirmedCases[j]
                end
            else
                if ismissing(_tdf.ConfirmedCases[j])
                    covidcases[j, i] = 0
                else
                    if ismissing(_tdf.ConfirmedCases[j-1])
                        covidcases[j, i] = _tdf.ConfirmedCases[j]
                    else
                        covidcases[j, i] = _tdf.ConfirmedCases[j] - _tdf.ConfirmedCases[j-1]
                    end
                end
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
end

W_coviddata = generatew_gt(COVIDSERIALINTERVAL, covidcases, POPULATION2020; blankn=10)


#asdf = [ RenewalDiffInDiff._generatew_gtrow_susceptibles(covidcases, .04, POPULATION2020[1], t, 1) for t in 1:133 ]


## Analysis 1 
# Effect of mask recommendations. No other considerations of confounding 

# everywhere has a recommendation by day 134 so limit analysis to first 133 days 

datamodel1 = diffindiffparameters_splinetimes(
    W_coviddata[1:133, :], 
    covidcases,
    facialcoveringsrecommended[1:133, :], 
    [ [ 1 ]; collect(75:182/4:257) ],
    POPULATION2020;
    psiprior=0.4,
)

datac1config = @ntuple modelname="datamodel1" model=datamodel1 n_rounds n_chains=8 seed=510+id
datachain1dict = produce_or_load(pol_fitparameter, datac1config, datadir("sims"))

## Analysis 2 
# Effect of mask requirements. No other considerations of confounding 

datamodel2 = diffindiffparameters_splinetimes(
    W_coviddata, 
    covidcases,
    facialcoveringsrequired, 
    [ [ 1 ]; collect(75:182/4:257) ],
    POPULATION2020;
    psiprior=0.4,
)

datac2config = @ntuple modelname="datamodel2" model=datamodel2 n_rounds n_chains=8 seed=520+id
datachain2dict = produce_or_load(pol_fitparameter, datac2config, datadir("sims"))

## Analysis 3 
# Effect of mask requirements with secondary interventions of end of stay-at-home and some
# businesses reopening

datamodel3 = diffindiffparameters_splinetimes(
    W_coviddata, 
    covidcases,
    facialcoveringsrequired, 
    [ [ 1 ]; collect(75:182/4:257) ],
    POPULATION2020;
    psiprior=0.4,
    secondaryinterventions=[ endstayathometimes, somebusinessreopen ],
)

datac3config = @ntuple modelname="datamodel3" model=datamodel3 n_rounds n_chains=8 seed=530+id
datachain3dict = produce_or_load(pol_fitparameter, datac3config, datadir("sims"))

## Analysis 4 
# Effect of mask requirements with the government response index as a secondary intervention

datamodel4 = diffindiffparameters_splinetimes(
    W_coviddata, 
    covidcases,
    facialcoveringsrequired, 
    [ [ 1 ]; collect(75:182/4:257) ],
    POPULATION2020;
    psiprior=0.4,
    secondaryinterventions=gri_nofc,
)

datac4config = @ntuple modelname="datamodel4" model=datamodel4 n_rounds n_chains=8 seed=540+id
datachain4dict = produce_or_load(pol_fitparameter, datac4config, datadir("sims"))

## Analysis 5 
# Add lead and lag times that groups can diverge 

secondaryinterventions_data = [ 
    interventionsoffset(facialcoveringsrequired, offsettimes[i]) 
    for i ∈ 1:6
]

datamodel5 = diffindiffparameters_splinetimes(
    W_coviddata, 
    covidcases,
    facialcoveringsrequired, 
    [ [ 1 ]; collect(75:182/4:257) ],
    POPULATION2020;
    psiprior=0.4,
    secondaryinterventions=[ 
        [ endstayathometimes, somebusinessreopen ]; secondaryinterventions_data 
    ],
)

datac5config = @ntuple modelname="datamodel5" model=datamodel5 n_rounds n_chains=8 seed=550+id
datachain5dict = produce_or_load(pol_fitparameter, datac5config, datadir("sims"))

datamodel6 = diffindiffparameters_splinetimes(
    W_coviddata, 
    covidcases,
    facialcoveringsrequired, 
    [ [ 1 ]; collect(75:182/4:257) ],
    POPULATION2020;
    psiprior=0.4,
    secondaryinterventions=[ 
        [ gri_nofc ]; secondaryinterventions_data 
    ],
)

datac6config = @ntuple modelname="datamodel6" model=datamodel6 n_rounds n_chains=8 seed=560+id
datachain6dict = produce_or_load(pol_fitparameter, datac6config, datadir("sims"))

## Analysis 7
# With recommendation as an alternative exposure 

datamodel7 = diffindiffparameters_splinetimes(
    W_coviddata, 
    covidcases,
    facialcoveringsrequired, 
    [ [ 1 ]; collect(75:182/4:257) ],
    POPULATION2020;
    psiprior=0.4,
    secondaryinterventions=[ 
        facialcoveringsrecommended, endstayathometimes, somebusinessreopen  
    ],
)

datac7config = @ntuple modelname="datamodel7" model=datamodel7 n_rounds n_chains=8 seed=570+id
datachain7dict = produce_or_load(pol_fitparameter, datac7config, datadir("sims"))

datamodel8 = diffindiffparameters_splinetimes(
    W_coviddata, 
    covidcases,
    facialcoveringsrequired, 
    [ [ 1 ]; collect(75:182/4:257) ],
    POPULATION2020;
    psiprior=0.4,
    secondaryinterventions=[ 
        facialcoveringsrecommended, gri_nofc 
    ],
)

datac8config = @ntuple modelname="datamodel8" model=datamodel8 n_rounds n_chains=8 seed=580+id
datachain8dict = produce_or_load(pol_fitparameter, datac8config, datadir("sims"))