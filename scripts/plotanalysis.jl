
using DrWatson
@quickactivate :RenewalDiffInDiff

using CairoMakie 
include("analysis.jl")
include(srcdir("plottingfunctions.jl"))

maxrounds = 12


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim1chain1 = loadanalysisdictsasdf("sim1model1", 8, maxrounds, 110)
plotchains(sim1chain1)
sim1fit1 = samplerenewalequation_2sets(
    f_seirvector, sim1chain1, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    #psi=0.8, 
    timeperiods=timeperiods2
)
sim1fit1kv = keyvalues(sim1chain1, sim1fit1)
sim1fit1plot = plotrenewalequationsamples(
    simulation1dataset, W_sim1, sim1fit1; 
    betafunctions=betafunctions1, 
    betafunctions_counterfactual=betafunctions1_counterfactual,
    infectiousduration=2.5,
    plotsize=( 400, 400 )
)

safesave(plotsdir("sim1fit1plot.svg"), sim1fit1plot)

sim1chain2 = loadanalysisdictsasdf("sim1model2", 8, maxrounds, 120)
plotchains(sim1chain2)
sim1fit2 = samplerenewalequation_2sets(
    f_seirvector, sim1chain2, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    #psi=0.8, 
    timeperiods=timeperiods5,
)
sim1fit2plot = plotrenewalequationsamples(
    simulation1dataset, W_sim1, sim1fit2; 
    betafunctions=betafunctions1, 
    betafunctions_counterfactual=betafunctions1_counterfactual,
    infectiousduration=2.5,
)

sim1chain3 = loadanalysisdictsasdf("sim1model3", 8, maxrounds, 130)
plotchains(sim1chain3)
sim1fit3 = samplerenewalequation_2sets(
    f_seirvector, sim1chain3, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    #psi=0.8, 
    timeperiods=timeperiods5,
    secondaryinterventions=secondaryinterventions_sim1,
)
sim1fit3plot = plotrenewalequationsamples(
    simulation1dataset, W_sim1, sim1fit3; 
    betafunctions=betafunctions1, 
    betafunctions_counterfactual=betafunctions1_counterfactual,
    infectiousduration=2.5,
)

sim1chain4 = loadanalysisdictsasdf("sim1model4", 8, maxrounds, 140)
plotchains(sim1chain4)
sim1fit4 = samplerenewalequation_2sets(
    f_seirvector, sim1chain4, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    #psi=0.8, 
    timeknots=[ [ 1 ]; collect(11:189/4:200) ],
)
sim1fit4plot = plotrenewalequationsamples(
    simulation1dataset, W_sim1, sim1fit4; 
    betafunctions=betafunctions1, 
    betafunctions_counterfactual=betafunctions1_counterfactual,
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    rhoclip=2,
)

safesave(plotsdir("sim1fit4plot.svg"), sim1fit4plot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 2 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim2chain1 = loadanalysisdictsasdf("sim2model1", 8, maxrounds, 210)
plotchains(sim2chain1)
sim2fit1 = samplerenewalequation_2sets(
    f_seirvector, sim2chain1, simulation2dataset["interventions"]; 
    initialvalues=simulation2dataset["cases"][1:10, :], 
    Ns=simulation2dataset["Ns"], 
    #psi=0.45, 
    timeperiods=timeperiods5,
)
sim2fit1plot = plotrenewalequationsamples(
    simulation2dataset, W_sim2, sim2fit1; 
    betafunctions=betafunctions2, 
    betafunctions_counterfactual=betafunctions2_counterfactual,
    infectiousduration=2.5,
)

sim2chain2 = loadanalysisdictsasdf("sim2model2", 8, maxrounds, 220)
plotchains(sim2chain2)
sim2fit2 = samplerenewalequation_2sets(
    f_seirvector, sim2chain2, simulation2dataset["interventions"]; 
    initialvalues=simulation2dataset["cases"][1:10, :], 
    Ns=simulation2dataset["Ns"], 
    #psi=0.45, 
    timeknots=[ [ 1 ]; collect(11:189/4:200) ],
)
sim2fit2plot = plotrenewalequationsamples(
    simulation2dataset, W_sim2, sim2fit2; 
    betafunctions=betafunctions2, 
    betafunctions_counterfactual=betafunctions2_counterfactual,
    infectiousduration=2.5,
)

sim2chain3 = loadanalysisdictsasdf("sim2model3", 8, maxrounds, 230)
plotchains(sim2chain3)
sim2fit3 = samplerenewalequation_2sets(
    f_seirvector, sim2chain3, simulation2dataset["interventions"]; 
    initialvalues=simulation2dataset["cases"][1:10, :], 
    Ns=simulation2dataset["Ns"], 
    #psi=0.45, 
    timeknots=[ [ 1 ]; collect(11:189/4:200) ],
)
sim2fit3plot = plotrenewalequationsamples(
    simulation2dataset, W_sim2, sim2fit3; 
    betafunctions=betafunctions2, 
    betafunctions_counterfactual=betafunctions2_counterfactual,
    infectiousduration=2.5,
    plotsize=( 400, 400 )
)

safesave(plotsdir("sim2fit3plot.svg"), sim2fit3plot)


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
    initialvalues=simulation3dataset["cases"][1:10, :], 
    Ns=simulation3dataset["Ns"], 
    #psi=0.3, 
    timeknots=[ [ 1 ]; collect(11:189/4:200) ],
)
sim3fit1plot = plotrenewalequationsamples(
    simulation3dataset, W_sim3, sim3fit1; 
    betafunctions=betafunctions3, 
    betafunctions_counterfactual=betafunctions3_counterfactual,
    infectiousduration=2.5,
)

sim3chain2 = loadanalysisdictsasdf("sim3model2", 8, maxrounds, 320)
plotchains(sim3chain2)
sim3fit2 = samplerenewalequation_2sets(
    f_seirvector, sim3chain2, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases"][1:10, :], 
    Ns=simulation3dataset["Ns"], 
    #psi=0.3, 
    timeknots=[ [ 1 ]; collect(11:189/4:200) ],
    secondaryinterventions=simulation3dataset["secondaryinterventions"],
)
sim3fit2plot = plotrenewalequationsamples(
    simulation3dataset, W_sim3, sim3fit2; 
    betafunctions=betafunctions3, 
    betafunctions_counterfactual=betafunctions3_counterfactual,
    infectiousduration=2.5,
)

sim3chain3 = loadanalysisdictsasdf("sim3model3", 8, maxrounds, 330)
plotchains(sim3chain3)
sim3fit3 = samplerenewalequation_2sets(
    f_seirvector, sim3chain3, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases"][1:10, :], 
    Ns=simulation3dataset["Ns"], 
    #psi=0.3, 
    timeknots=[ [ 1 ]; collect(11:189/4:200) ],
    secondaryinterventions=secondaryinterventions_sim3,
)
sim3fit3plot = plotrenewalequationsamples(
    simulation3dataset, W_sim3, sim3fit3; 
    betafunctions=betafunctions3, 
    betafunctions_counterfactual=betafunctions3_counterfactual,
    infectiousduration=2.5,
)

sim3chain4 = loadanalysisdictsasdf("sim3model4", 8, maxrounds, 340)
plotchains(sim3chain4)
sim3fit4 = samplerenewalequation_2sets(
    f_seirvector, sim3chain4, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases"][1:10, :], 
    Ns=simulation3dataset["Ns"], 
    #psi=0.3, 
    timeknots=[ [ 1 ]; collect(11:189/4:200) ],
    secondaryinterventions=secondaryinterventions_sim3d,
)
sim3fit4plot = plotrenewalequationsamples(
    simulation3dataset, W_sim3, sim3fit4; 
    betafunctions=betafunctions3, 
    betafunctions_counterfactual=betafunctions3_counterfactual,
    infectiousduration=2.5,
    plotsize=( 400, 400 )
)

safesave(plotsdir("sim2fit3plot.svg"), sim2fit3plot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 4 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim4chain1 = loadanalysisdictsasdf("sim4model1", 8, maxrounds, 410)
plotchains(sim4chain1)
sim4fit1 = samplerenewalequation_2sets(
    f_seirvector, sim4chain1, simulation4dataset["interventions"]; 
    initialvalues=simulation4dataset["cases"][1:10, :], 
    Ns=simulation4dataset["Ns"], 
    #psi=0.45, 
    timeknots=[ [ 1 ]; collect(11:189/4:200) ],
)
sim4fit1plot = plotrenewalequationsamples(
    simulation4dataset, W_sim4, sim4fit1; 
    betafunctions=betafunctions4, 
    betafunctions_counterfactual=betafunctions4_counterfactual,
    infectiousduration=2.5,
)

sim4chain2 = loadanalysisdictsasdf("sim4model2", 8, maxrounds, 420)
plotchains(sim4chain2)
sim4fit2 = samplerenewalequation_2sets(
    f_seirvector, sim4chain2, simulation4dataset["interventions"]; 
    initialvalues=simulation4dataset["cases"][1:10, :], 
    Ns=simulation4dataset["Ns"], 
    #psi=0.45, 
    timeknots=[ [ 1 ]; collect(11:189/4:200) ],
    secondaryinterventions=secondaryinterventions_sim4,
)
sim4fit2plot = plotrenewalequationsamples(
    simulation4dataset, W_sim4, sim4fit2; 
    betafunctions=betafunctions4, 
    betafunctions_counterfactual=betafunctions4_counterfactual,
    infectiousduration=2.5,
)

sim4chain3 = loadanalysisdictsasdf("sim4model3", 8, maxrounds, 430)
plotchains(sim4chain3)
sim4fit3 = samplerenewalequation_2sets(
    f_seirvector, sim4chain3, simulation4dataset["interventions"]; 
    initialvalues=simulation4dataset["cases"][1:10, :], 
    Ns=simulation4dataset["Ns"], 
    #psi=0.9, 
    timeknots=[ [ 1 ]; collect(11:189/4:200) ],
)
sim4fit3plot = plotrenewalequationsamples(
    simulation4dataset, W_sim4, sim4fit3; 
    betafunctions=betafunctions4, 
    betafunctions_counterfactual=betafunctions4_counterfactual,
    infectiousduration=2.5,
)

safesave(plotsdir("sim4fit3plot.svg"), sim4fit3plot)



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
    plotsize=( 400, 400 )
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
