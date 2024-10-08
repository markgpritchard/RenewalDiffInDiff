
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
# Simulation 3 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim3chain1 = loadanalysisdictsasdf("sim3model1", 8, maxrounds, 40)
plotchains(sim3chain1)
sim3fit1 = samplerenewalequation_2sets(
    f_seirvector, sim3chain1, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases"][1:20, :], 
    Ns=simulation3dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim3fit1plot = plotrenewalequationsamples(
    simulation3dataset, W_sim3, sim3fit1; 
    betafunctions=[ beta3a, beta3b ], 
    betafunctions_counterfactual=[ beta3a, beta3bcounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    rhoclip=2.5,
)

safesave(plotsdir("sim3fit1plot.svg"), sim3fit1plot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 4 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim4chain1 = loadanalysisdictsasdf("sim4model1", 8, maxrounds, 50)
plotchains(sim4chain1)
sim4fit1 = samplerenewalequation_2sets(
    f_seirvector, sim4chain1, simulation4dataset["interventions"]; 
    initialvalues=simulation4dataset["cases"][1:20, :], 
    Ns=simulation4dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim4fit1plot = plotrenewalequationsamples(
    simulation4dataset, W_sim4, sim4fit1; 
    betafunctions=[ beta4a, beta4b ], 
    betafunctions_counterfactual=[ beta4a, beta4bcounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    rhoclip=2.5,
)

safesave(plotsdir("sim4fit1plot.svg"), sim4fit1plot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Covid data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

datachain1 = loadanalysisdictsasdf("datamodel1", 8, maxrounds, 3010)
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
    plotsize=( 600, 400 ),
    columntitles=[ 
        "Halton", 
        "Knowsley", 
        "Liverpool", 
        "Sefton", 
        "St Helens", 
        "Warrington", 
        "W. Lancs", 
        "Wigan", 
        "Wirral" 
    ],
    columntitlefontsize=10,
    xticklabelrotation=-π/4,
    xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
    xtitle="Date, 2020–2021",
)

safesave(plotsdir("datafit1plot.svg"), datafit1plot)

# force 28% reduction 

reduceddatafit1 = deepcopy(datafit1)
for g ∈ 1:9, t ∈ 1:216 
    if masstesting[t, g] == 1 
        for x ∈ eachindex(reduceddatafit1.y_matrix_poisson_vec)
            reduceddatafit1.y_matrix_poisson_vec[x][t, g] = round(
                Int, reduceddatafit1.y_matrix_poisson_vec[x][t, g] / 1.28
            )
        end
    end 
end

datafit1plotreduced = plotrenewalequationsamples(
    allcovidcases, W_allcoviddata, selectpops, reduceddatafit1;
    plotsize=( 600, 400 ),
    columntitles=[ 
        "Halton", 
        "Knowsley", 
        "Liverpool", 
        "Sefton", 
        "St Helens", 
        "Warrington", 
        "W. Lancs", 
        "Wigan", 
        "Wirral" 
    ],
    columntitlefontsize=10,
    xticklabelrotation=-π/4,
    xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
    xtitle="Date, 2020–2021",
)





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

safesave(plotsdir("datafit3plot.svg"), datafit3plot)

datachain5 = loadanalysisdictsasdf("datamodel5", 8, maxrounds, 540)
plotchains(datachain5)
datafit5 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain5, masstesting; 
    initialvalues=allcovidcases[1:49, :], 
    Ns=selectpops,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=collect(1.0:28:216),
    secondaryinterventions=[ 
        InterventionsMatrix([ 217, 217, 146, 217, 217, 217], 261), 
        InterventionsMatrix([ 217, 217, 174, 217, 217, 217], 261), 
    ]
)
datafit5plot = plotrenewalequationsamples(
    allcovidcases, W_allcoviddata2, selectpops, datafit5;
    plotsize=( 400, 400 ),
    columntitles=[ "Halton", "Knowsley", "Liverpool", "Sefton", "St Helens", "Wirral" ],
    columntitlefontsize=10,
    xticklabelrotation=-π/4,
    xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
    xtitle="Date, 2020–2021",
)

safesave(plotsdir("datafit5plot.svg"), datafit5plot)
