
using DrWatson
@quickactivate :RenewalDiffInDiff

using CairoMakie 
include("analysis.jl")
include(srcdir("plottingfunctions.jl"))

maxrounds = 12


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim1chain0 = loadanalysisdictsasdf("sim1model0", 8, maxrounds, 10)
plotchains(sim1chain0)
sim1fit0 = samplerenewalequation_2sets(
    f_seirvector, sim1chain0, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases_counterfactual"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim1fit0kv = keyvalues(sim1chain0, sim1fit0)
sim1fit0plot = plotrenewalequationsamples(
    simulation1dataset, W_sim1_0, sim1fit0; 
    betafunctions=[ beta1a, beta1bcounterfactual ], 
    betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
)

safesave(plotsdir("sim1fit0plot.svg"), sim1fit0plot)

sim1chain1 = loadanalysisdictsasdf("sim1model1", 8, maxrounds, 20)
plotchains(sim1chain1)
sim1fit1 = samplerenewalequation_2sets(
    f_seirvector, sim1chain1, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim1fit1kv = keyvalues(sim1chain1, sim1fit1)
sim1fit1plot = plotrenewalequationsamples(
    simulation1dataset, W_sim1, sim1fit1; 
    betafunctions=[ beta1a, beta1b ], 
    betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    rhoclip = 2
)

safesave(plotsdir("sim1fit1plot.svg"), sim1fit1plot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 2 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim2chain1 = loadanalysisdictsasdf("sim2model1", 8, maxrounds, 30)
plotchains(sim2chain1)
sim2fit1 = samplerenewalequation_2sets(
    f_seirvector, sim2chain1, simulation2dataset["interventions"]; 
    initialvalues=simulation2dataset["cases"][1:20, :], 
    Ns=simulation2dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
    secondaryinterventions=InterventionsMatrix([ nothing, nothing, 30 ], 100),
)
sim2fit1plot = plotrenewalequationsamples(
    simulation2dataset, W_sim2, sim2fit1; 
    betafunctions=[ beta2a, beta2b, beta2c ], 
    betafunctions_counterfactual=[ beta2a, beta2bcounterfactual, beta2ccounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    rhoclip=2.5,
)

safesave(plotsdir("sim2fit1plot.svg"), sim2fit1plot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Covid data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

datachain1 = loadanalysisdictsasdf("datamodel1", 8, maxrounds, 510)
plotchains(datachain1)
datafit1 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain1, masstesting; 
    initialvalues=allcovidcases[1:100, :], 
    Ns=selectpops,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=collect(1.0:28:216),
)
datafit1plot = plotrenewalequationsamples(
    allcovidcases, W_allcoviddata, selectpops, datafit1;
    plotsize=( 400, 400 ),
    columntitles=[ "Halton", "Knowsley", "Liverpool", "Sefton", "St Helens", "Wirral" ],
    columntitlefontsize=10,
    xticklabelrotation=-π/4,
    xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
    xtitle="Date, 2020–2021",
)

safesave(plotsdir("datafit1plot.svg"), datafit1plot)


greater = zeros(Int, size(datafit1.y_matrix_poisson_vec[1]))
lesser = zeros(Int, size(datafit1.y_matrix_poisson_vec[1]))
for i ∈ eachindex(datafit1.y_matrix_poisson_vec)
    vs = cumsum(datafit1.y_matrix_poisson_vec[i]; dims=1)
    cfvs = cumsum(datafit1.y_matrix_poisson_vec_counterfactual[i]; dims=1)
    for t ∈ axes(greater, 1), g ∈ axes(greater, 2)
        if vs[t, g] < cfvs[t, g]
            lesser[t, g] += 1 
        elseif vs[t, g] > cfvs[t, g]
            greater[t, g] += 1 
        end
    end
end

axs =[ Axis(datafit1plot[5, i]) for i ∈ 1:6 ]
for i ∈ 1:6 
    band!(
        axs[i],
        axes(lesser, 1),
        zeros(size(lesser, 1)),
        greater[:, i] ./ length(datafit1.y_matrix_poisson_vec);
        color=( :red, 0.5 ),
    )
    band!(
        axs[i],
        axes(lesser, 1),
        -lesser[:, i] ./ length(datafit1.y_matrix_poisson_vec),
        zeros(size(lesser, 1));
        color=( :green, 0.5 ),
    )
end

datafit1plot




datachain1discrete = loadanalysisdictsasdf("datamodel1discrete", 8, maxrounds, 515)
plotchains(datachain1discrete)
datafit1discrete = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain1discrete, masstesting; 
    initialvalues=allcovidcases[1:49, :], 
    Ns=selectpops,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeperiods=dtimes,
)
datafit1discreteplot = plotrenewalequationsamples(
    allcovidcases, W_allcoviddata, selectpops, datafit1discrete
)

datachain2 = loadanalysisdictsasdf("datamodel2", 8, maxrounds, 520)
plotchains(datachain2)
datafit2 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain2, masstesting; 
    initialvalues=pil1covidcases[1:124, :], 
    Ns=selectpops,
    timeknots=collect(1.0:28:216),
    #timeknots=collect(1:303/10:304),
    #secondaryinterventions=secondaryinterventions_data
)
datafit2plot = plotrenewalequationsamples(
    pil1covidcases, W_pil1coviddata, selectpops, datafit2
)

datachain2discrete = loadanalysisdictsasdf("datamodel2discrete", 8, maxrounds, 525)
plotchains(datachain2discrete)
datafit2discrete = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain2discrete, masstesting; 
    initialvalues=pil1covidcases[1:124, :], 
    Ns=selectpops,
    timeperiods=dtimes,
    #timeknots=collect(1:303/10:304),
    #secondaryinterventions=secondaryinterventions_data
)
datafit2discreteplot = plotrenewalequationsamples(
    pil1covidcases, W_pil1coviddata, selectpops, datafit2discrete
)

datachain3 = loadanalysisdictsasdf("datamodel3", 8, maxrounds, 530)
plotchains(datachain3)
datafit3 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain3, masstesting2; 
    initialvalues=allcovidcases2[1:100, :], 
    Ns=selectpops,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=collect(1.0:28:317),
    secondaryinterventions=universaltesting,
)
datafit3plot = plotrenewalequationsamples(
    allcovidcases2, W_allcoviddata2, selectpops, datafit3;
    plotsize=( 400, 400 ),
    columntitles=[ "Halton", "Knowsley", "Liverpool", "Sefton", "St Helens", "Wirral" ],
    columntitlefontsize=10,
    xticklabelrotation=-π/4,
    xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
    xtitle="Date, 2020–2021",
)

safesave(plotsdir("datafit1plot.svg"), datafit1plot)
