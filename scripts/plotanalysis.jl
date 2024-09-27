
using DrWatson
@quickactivate :RenewalDiffInDiff

using CairoMakie 
include("analysis.jl")
include(srcdir("plottingfunctions.jl"))

maxrounds = 12

###########################
###########################

#=
fig = Figure()
axs = [ Axis(fig[i, j]) for i ∈ 1:4, j ∈ 1:4 ]

for g ∈ 1:4
    _tdf = filter(:RegionId => x -> x == g, coviddf)
    for (i, v) ∈ enumerate([
        :C1E_Schoolclosing,
        :C2E_Workplaceclosing,
        :C3E_Cancelpublicevents,
        :C4E_Restrictionsongatherings,
        :C5E_Closepublictransport,
        :C6E_Stayathome,
        :C7E_Restrictionsoninternalmovement,
        :C8E_Internationaltravelcontrols,
        :E1E_Incomesupport,
        :E2E_Debtcontractrelief,
        :H1E_Publicinformationcampaigns,
        :H2E_Testingpolicy,
        :H3E_Contacttracing,
        :H6E_FacialCoverings,
        :H7E_Vaccinationpolicy,
        :H8E_Protectionofelderlypeople,
    ])
        lines!(axs[i], axes(_tdf, 1), getproperty(_tdf, v); color=COLOURVECTOR[g])
    end
end

fig
=#

###########################
###########################


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function loadanalysisdictsasdf(modelname, n_chains, maxrounds, seedstart; nsims=4)
    fn1 = _findanalysisfilename(modelname, n_chains, maxrounds, seedstart + 1)
    chain1 = load(fn1)["chain"]
    df = DataFrame(chain1)
    for i ∈ 2:nsims 
        fn = _findanalysisfilename(modelname, n_chains, maxrounds, seedstart + i)
        isnothing(fn) && continue
        chain = load(fn)["chain"]
        _tdf = DataFrame(chain)
        _tdf.chain = [ i for _ ∈ axes(_tdf, 1) ]
        df = vcat(df, _tdf) 
    end
    return df
end

function _findanalysisfilename(modelname, n_chains, maxrounds, seed)
    n_rounds = maxrounds
    while true 
        filename = "modelname=$(modelname)_n_chains=$(n_chains)_n_rounds=$(n_rounds)_seed=$(seed).jld2"
        if isfile(datadir("sims", filename))
            return datadir("sims", filename)
        else 
            n_rounds += -1 
            if n_rounds <= 0 return nothing end
        end
    end
end

#=
sim1chain0 = loadanalysisdictsasdf("sim1model0", 8, maxrounds, 100)
plotchains(sim1chain0)
sim1fit0 = samplerenewalequation_2sets(
    f_seirvector, sim1chain0, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases_counterfactual"][1:10, :], Ns=simulation1dataset["Ns"], 
    timeperiods=timeperiods2
)
sim1fit0plot = plotrenewalequationsamples(
    simulation1dataset, W_sim1_0, sim1fit0; 
    betafunctions=betafunctions1_counterfactual, betafunctions_counterfactual=betafunctions1_counterfactual,
    infectiousduration=2.5,
)
=#

function keyvalues(fitteddf, fittedvaluesset)
    deltamean = exp(mean(fitteddf.logdelta))
    deltap05_95 = exp.(quantile(fitteddf.logdelta, [ 0.05, 0.95 ]))
    totalcases = Matrix{Int}(
        undef, 
        length(fittedvaluesset.y_matrix_poisson_vec), 
        size(fittedvaluesset.y_matrix_poisson_vec[1], 2)
    )
    peakcases = Matrix{Int}(
        undef, 
        length(fittedvaluesset.y_matrix_poisson_vec),
        size(fittedvaluesset.y_matrix_poisson_vec[1], 2)
    )
    peakcasesdate = Matrix{Int}(
        undef, 
        length(fittedvaluesset.y_matrix_poisson_vec), 
        size(fittedvaluesset.y_matrix_poisson_vec[1], 2)
    )
    for i ∈ eachindex(fittedvaluesset.y_matrix_poisson_vec)
        totalcases[i, :] = [
            sum(fittedvaluesset.y_matrix_poisson_vec[i][:, j]) - 
                sum(fittedvaluesset.y_matrix_poisson_vec_counterfactual[i][:, j])
            for j ∈ axes(fittedvaluesset.y_matrix_poisson_vec[i], 2)
        ]
        for j ∈  axes(fittedvaluesset.y_matrix_poisson_vec[i], 2)
            x, ind = findmax(fittedvaluesset.y_matrix_poisson_vec[i][:, j])
            x_cf, ind_cf = findmax(
                fittedvaluesset.y_matrix_poisson_vec_counterfactual[i][:, j]
            )
            peakcases[i, j] = x - x_cf 
            peakcasesdate[i, j] = ind - ind_cf
        end
    end
    totalcasesmeaneffect = [ mean(totalcases[:, i]) for i ∈ axes(totalcases, 2) ]
    totalcasesp05_90 = [ 
        quantile(totalcases[:, i], [ 0.05, 0.95 ]) 
        for i ∈ axes(totalcases, 2) 
    ]
    peakcasesmeaneffect = [ mean(peakcases[:, i]) for i ∈ axes(peakcases, 2) ]
    peakcasesp05_90 = [ quantile(peakcases[:, i], [ 0.05, 0.95 ]) for i ∈ axes(peakcases, 2) ]
    peakcasesdatemeaneffect = [ mean(peakcasesdate[:, i]) for i ∈ axes(peakcasesdate, 2) ]
    peakcasesdatep05_90 = [ 
        quantile(peakcasesdate[:, i], [ 0.05, 0.95 ]) 
        for i ∈ axes(peakcasesdate, 2) 
    ]
    return (
        deltamean=deltamean,
        deltap05_95=deltap05_95,
        totalcases=totalcases,
        peakcases=peakcases,
        peakcasesdate=peakcasesdate,
        totalcasesmeaneffect=totalcasesmeaneffect,
        totalcasesp05_90=totalcasesp05_90,
        peakcasesmeaneffect=peakcasesmeaneffect,
        peakcasesp05_90=peakcasesp05_90,
        peakcasesdatemeaneffect=peakcasesdatemeaneffect,
        peakcasesdatep05_90=peakcasesdatep05_90,
    )
end

sim1chain1 = loadanalysisdictsasdf("sim1model1", 8, maxrounds, 110)
plotchains(sim1chain1)
sim1fit1 = samplerenewalequation_2sets(
    f_seirvector, sim1chain1, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], Ns=simulation1dataset["Ns"], 
    psi=0.8, timeperiods=timeperiods2
)
sim1fit1kv = keyvalues(sim1chain1, sim1fit1)
sim1fit1plot = plotrenewalequationsamples(
    simulation1dataset, W_sim1, sim1fit1; 
    betafunctions=betafunctions1, betafunctions_counterfactual=betafunctions1_counterfactual,
    infectiousduration=2.5,
)

sim1chain2 = loadanalysisdictsasdf("sim1model2", 8, maxrounds, 120)
plotchains(sim1chain2)
sim1fit2 = samplerenewalequation_2sets(
    f_seirvector, sim1chain2, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], Ns=simulation1dataset["Ns"], 
    psi=0.8, timeperiods=timeperiods5,
)
sim1fit2plot = plotrenewalequationsamples(
    simulation1dataset, W_sim1, sim1fit2; 
    betafunctions=betafunctions1, betafunctions_counterfactual=betafunctions1_counterfactual,
    infectiousduration=2.5,
)

sim1chain3 = loadanalysisdictsasdf("sim1model3", 8, maxrounds, 130)
plotchains(sim1chain3)
sim1fit3 = samplerenewalequation_2sets(
    f_seirvector, sim1chain3, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], Ns=simulation1dataset["Ns"], 
    psi=0.8, timeperiods=timeperiods5,
    secondaryinterventions=secondaryinterventions_sim1,
)
sim1fit3plot = plotrenewalequationsamples(
    simulation1dataset, W_sim1, sim1fit3; 
    betafunctions=betafunctions1, betafunctions_counterfactual=betafunctions1_counterfactual,
    infectiousduration=2.5,
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 2 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim2chain1 = loadanalysisdictsasdf("sim2model1", 8, maxrounds, 210)
plotchains(sim2chain1)
sim2fit1 = samplerenewalequation_2sets(
    f_seirvector, sim2chain1, simulation2dataset["interventions"]; 
    initialvalues=simulation2dataset["cases"][1:10, :], Ns=simulation2dataset["Ns"], 
    psi=0.45, timeperiods=timeperiods5,
)
sim2fit1plot = plotrenewalequationsamples(
    simulation2dataset, W_sim2, sim2fit1; 
    betafunctions=betafunctions2, betafunctions_counterfactual=betafunctions2_counterfactual,
    infectiousduration=2.5,
)

sim2chain2 = loadanalysisdictsasdf("sim2model2", 8, maxrounds, 220)
plotchains(sim2chain2)
sim2fit2 = samplerenewalequation_2sets(
    f_seirvector, sim2chain2, simulation2dataset["interventions"]; 
    initialvalues=simulation2dataset["cases"][1:10, :], Ns=simulation2dataset["Ns"], 
    psi=0.45, timeknots=[ [ 1 ]; collect(11:189/4:200) ],
)
sim2fit2plot = plotrenewalequationsamples(
    simulation2dataset, W_sim2, sim2fit2; 
    betafunctions=betafunctions2, betafunctions_counterfactual=betafunctions2_counterfactual,
    infectiousduration=2.5,
)

sim2chain3 = loadanalysisdictsasdf("sim2model3", 8, maxrounds, 230)
plotchains(sim2chain3)
sim2fit3 = samplerenewalequation_2sets(
    f_seirvector, sim2chain3, simulation2dataset["interventions"]; 
    initialvalues=simulation2dataset["cases"][1:10, :], Ns=simulation2dataset["Ns"], 
    psi=0.45, timeknots=[ [ 1 ]; collect(11:189/4:200) ],
)
sim2fit3plot = plotrenewalequationsamples(
    simulation2dataset, W_sim2, sim2fit3; 
    betafunctions=betafunctions2, betafunctions_counterfactual=betafunctions2_counterfactual,
    infectiousduration=2.5,
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 3 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#=
sim3chain0 = loadanalysisdictsasdf("sim3model0", 8, maxrounds, 300)
plotchains(sim3chain0)
sim3fit0 = samplerenewalequation_2sets(
    f_seirvector, sim3chain0, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases_counterfactual"][1:10, :], Ns=simulation3dataset["Ns"], 
    psi=0.3, timeknots=[ [ 1 ]; collect(11:189/4:200) ],
)
sim3fit0plot = plotrenewalequationsamples(
    simulation3dataset, W_sim3_0, sim3fit0; 
    betafunctions=betafunctions3_counterfactual, betafunctions_counterfactual=betafunctions3_counterfactual,
    infectiousduration=2.5,
)
=#
sim3chain1 = loadanalysisdictsasdf("sim3model1", 8, maxrounds, 310)
plotchains(sim3chain1)
sim3fit1 = samplerenewalequation_2sets(
    f_seirvector, sim3chain1, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases"][1:10, :], Ns=simulation3dataset["Ns"], 
    psi=0.3, timeknots=[ [ 1 ]; collect(11:189/4:200) ],
)
sim3fit1plot = plotrenewalequationsamples(
    simulation3dataset, W_sim3, sim3fit1; 
    betafunctions=betafunctions3, betafunctions_counterfactual=betafunctions3_counterfactual,
    infectiousduration=2.5,
)

sim3chain2 = loadanalysisdictsasdf("sim3model2", 8, maxrounds, 320)
plotchains(sim3chain2)
sim3fit2 = samplerenewalequation_2sets(
    f_seirvector, sim3chain2, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases"][1:10, :], Ns=simulation3dataset["Ns"], 
    psi=0.3, timeknots=[ [ 1 ]; collect(11:189/4:200) ],
    secondaryinterventions=simulation3dataset["secondaryinterventions"],
)
sim3fit2plot = plotrenewalequationsamples(
    simulation3dataset, W_sim3, sim3fit2; 
    betafunctions=betafunctions3, betafunctions_counterfactual=betafunctions3_counterfactual,
    infectiousduration=2.5,
)

sim3chain3 = loadanalysisdictsasdf("sim3model3", 8, maxrounds, 330)
plotchains(sim3chain3)
sim3fit3 = samplerenewalequation_2sets(
    f_seirvector, sim3chain3, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases"][1:10, :], Ns=simulation3dataset["Ns"], 
    psi=0.3, timeknots=[ [ 1 ]; collect(11:189/4:200) ],
    secondaryinterventions=secondaryinterventions_sim3,
)
sim3fit3plot = plotrenewalequationsamples(
    simulation3dataset, W_sim3, sim3fit3; 
    betafunctions=betafunctions3, betafunctions_counterfactual=betafunctions3_counterfactual,
    infectiousduration=2.5,
)

sim3chain4 = loadanalysisdictsasdf("sim3model4", 8, maxrounds, 340)
plotchains(sim3chain4)
sim3fit4 = samplerenewalequation_2sets(
    f_seirvector, sim3chain4, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases"][1:10, :], Ns=simulation3dataset["Ns"], 
    psi=0.3, timeknots=[ [ 1 ]; collect(11:189/4:200) ],
    secondaryinterventions=secondaryinterventions_sim3d,
)
sim3fit4plot = plotrenewalequationsamples(
    simulation3dataset, W_sim3, sim3fit4; 
    betafunctions=betafunctions3, betafunctions_counterfactual=betafunctions3_counterfactual,
    infectiousduration=2.5,
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 4 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim4chain1 = loadanalysisdictsasdf("sim4model1", 8, maxrounds, 410)
plotchains(sim4chain1)
sim4fit1 = samplerenewalequation_2sets(
    f_seirvector, sim4chain1, simulation4dataset["interventions"]; 
    initialvalues=simulation4dataset["cases"][1:10, :], Ns=simulation4dataset["Ns"], 
    psi=0.45, timeknots=[ [ 1 ]; collect(11:189/4:200) ],
)
sim4fit1plot = plotrenewalequationsamples(
    simulation4dataset, W_sim4, sim4fit1; 
    betafunctions=betafunctions4, betafunctions_counterfactual=betafunctions4_counterfactual,
    infectiousduration=2.5,
)

sim4chain2 = loadanalysisdictsasdf("sim4model2", 8, maxrounds, 420)
plotchains(sim4chain2)
sim4fit2 = samplerenewalequation_2sets(
    f_seirvector, sim4chain2, simulation4dataset["interventions"]; 
    initialvalues=simulation4dataset["cases"][1:10, :], Ns=simulation4dataset["Ns"], 
    psi=0.45, timeknots=[ [ 1 ]; collect(11:189/4:200) ],
    secondaryinterventions=secondaryinterventions_sim4,
)
sim4fit2plot = plotrenewalequationsamples(
    simulation4dataset, W_sim4, sim4fit2; 
    betafunctions=betafunctions4, betafunctions_counterfactual=betafunctions4_counterfactual,
    infectiousduration=2.5,
)

sim4chain3 = loadanalysisdictsasdf("sim4model3", 8, maxrounds, 430)
plotchains(sim4chain3)
sim4fit3 = samplerenewalequation_2sets(
    f_seirvector, sim4chain3, simulation4dataset["interventions"]; 
    initialvalues=simulation4dataset["cases"][1:10, :], Ns=simulation4dataset["Ns"], 
    psi=0.9, timeknots=[ [ 1 ]; collect(11:189/4:200) ],
)
sim4fit3plot = plotrenewalequationsamples(
    simulation4dataset, W_sim4, sim4fit3; 
    betafunctions=betafunctions4, betafunctions_counterfactual=betafunctions4_counterfactual,
    infectiousduration=2.5,
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Covid data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

datachain1 = loadanalysisdictsasdf("datamodel1", 8, maxrounds, 510)
plotchains(datachain1)
datafit1 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain1, masstesting; 
    initialvalues=allcovidcases[1:99, :], Ns=selectpops,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=collect(1.0:28:216),
)
datafit1plot = plotrenewalequationsamples(allcovidcases, W_allcoviddata, selectpops, datafit1)

datachain2 = loadanalysisdictsasdf("datamodel2", 8, maxrounds, 520)
plotchains(datachain2)
datafit2 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain2, masstesting; 
    initialvalues=pil1covidcases[1:99, :], Ns=selectpops,
    timeknots=collect(1.0:28:216),
    #timeknots=collect(1:303/10:304),
    #secondaryinterventions=secondaryinterventions_data
)
datafit2plot = plotrenewalequationsamples(pil1covidcases, W_pil1coviddata, selectpops, datafit2)

datachain3 = loadanalysisdictsasdf("datamodel3", 8, maxrounds, 530)
plotchains(datachain3)
datafit3 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain3, nw_masstesting; 
    initialvalues=covidcases[1:99, :], Ns=populations,
    psi=0.4, timeknots=collectcollect(1:215/10:216),
)
datafit3plot = plotrenewalequationsamples(nw_covidcases, nw_W_coviddata, populations, datafit3)

datachain4 = loadanalysisdictsasdf("datamodel4", 8, maxrounds, 540)
plotchains(datachain4)
datafit4 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain4, facialcoveringsrequired; 
    initialvalues=covidcases[1:75, :], Ns=POPULATION2020, 
    psi=0.4, timeknots=[ [ 1 ]; collect(75:182/4:257) ],
    secondaryinterventions=gri_nofc,
)
datafit4plot = plotrenewalequationsamples(covidcases, W_coviddata, POPULATION2020, datafit4)

datachain5 = loadanalysisdictsasdf("datamodel5", 8, maxrounds, 550)
plotchains(datachain5)
datafit5 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain5, facialcoveringsrequired; 
    initialvalues=covidcases[1:75, :], Ns=POPULATION2020, 
    psi=0.4, timeknots=[ [ 1 ]; collect(75:182/4:257) ],
    secondaryinterventions=[ 
        [ endstayathometimes, somebusinessreopen ]; secondaryinterventions_data 
    ],
)
datafit5plot = plotrenewalequationsamples(covidcases, W_coviddata, POPULATION2020, datafit5)

datachain6 = loadanalysisdictsasdf("datamodel6", 8, maxrounds, 560)
plotchains(datachain6)
datafit6 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain6, facialcoveringsrequired; 
    initialvalues=covidcases[1:75, :], Ns=POPULATION2020, 
    psi=0.4, timeknots=[ [ 1 ]; collect(75:182/4:257) ],
    secondaryinterventions=[ 
        [ gri_nofc ]; secondaryinterventions_data 
    ],
)
datafit6plot = plotrenewalequationsamples(covidcases, W_coviddata, POPULATION2020, datafit6)

datachain7 = loadanalysisdictsasdf("datamodel7", 8, maxrounds, 570)
plotchains(datachain7)
datafit7 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain7, facialcoveringsrequired; 
    initialvalues=covidcases[1:75, :], Ns=POPULATION2020, 
    psi=0.4, timeknots=[ [ 1 ]; collect(75:182/4:257) ],
    secondaryinterventions=[ 
        facialcoveringsrecommended, endstayathometimes, somebusinessreopen  
    ],
)
datafit7plot = plotrenewalequationsamples(covidcases, W_coviddata, POPULATION2020, datafit7)

datachain8 = loadanalysisdictsasdf("datamodel8", 8, maxrounds, 580)
plotchains(datachain8)
datafit8 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain8, facialcoveringsrequired; 
    initialvalues=covidcases[1:75, :], Ns=POPULATION2020, 
    psi=0.4, timeknots=[ [ 1 ]; collect(75:182/4:257) ],
    secondaryinterventions=[ 
        facialcoveringsrecommended, gri_nofc 
    ],
)
datafit8plot = plotrenewalequationsamples(covidcases, W_coviddata, POPULATION2020, datafit8)

datachain9 = loadanalysisdictsasdf("datamodel9", 8, maxrounds, 590)
plotchains(datachain9)
datafit9 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain9, facialcoveringsrequired; 
    initialvalues=covidcases[1:75, :], Ns=POPULATION2020, 
    #initialvalues=covidcases[1:80, :], Ns=POPULATION2020, 
    psi=0.4, timeknots=[ [ 1 ]; collect(75:182/4:257) ],
    secondaryinterventions=[ 
        facialcoveringsrecommended, 
        MC1E_1, MC1E_2, MC1E_3, MC6E_1, MC6E_2, MC7E_1, MC7E_2, MH3E_1, MH8E_1
    ],
)
datafit9plot = plotrenewalequationsamples(covidcases, W_coviddata, POPULATION2020, datafit9)

datachain10 = loadanalysisdictsasdf("datamodel10", 8, maxrounds, 600)
plotchains(datachain10)
datafit10 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain10, facialcoveringsrecommended[1:191, :]; 
    initialvalues=covidcases[1:75, :], Ns=POPULATION2020, 
    psi=0.4, timeknots=[ [ 1 ]; collect(75:116/4:191) ],
    secondaryinterventions=[ 
        MC1E_1, MC1E_2, MC1E_3, MC6E_1, MC6E_2, MC7E_1, MC7E_2, MH3E_1, MH8E_1
    ],
)
datafit10plot = plotrenewalequationsamples(
    covidcases[1:191, :], W_coviddata191[1:191, :], POPULATION2020, datafit10
)
