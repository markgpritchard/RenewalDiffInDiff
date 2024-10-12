
using DrWatson
@quickactivate :RenewalDiffInDiff

using CairoMakie, StatsBase
include("analysis.jl")
include(srcdir("plottingfunctions.jl"))

maxrounds = 12


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot of generation interval
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

generationintervalplot = let
    fig = Figure(; size=( 500, 250 ))
    ax = Axis(fig[1, 1])
    scatter!(ax, 0:20, [ fseir(t) for t ∈ 0:20 ]; color=:black)
    formataxis!(ax)

    Label(fig[1, 0], "generation interval, f(τ)"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(fig[2, 1], "time, τ"; fontsize=11.84, tellwidth=false)

    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)
    fig
end

safesave(plotsdir("generationintervalplot.pdf"), generationintervalplot)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Discrete-time analyses

### Priors 

sim1model0discretepriorsample = sample(sim1model0discrete, Prior(), 16384)
sim1model0discretepriorsampledf = DataFrame(sim1model0discretepriorsample)
plotchains(sim1model0discretepriorsampledf)
sim1fit0discretepriorsample = samplerenewalequation_2sets(
    fseir, sim1model0discretepriorsampledf, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases_counterfactual"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    timeperiods=timeperiods2,
)
sim1fit0kvdiscretepriorsample = keyvalues(sim1model0discretepriorsampledf, sim1fit0discretepriorsample)
println(sim1fit0kvdiscretepriorsample)

sim1model1discretepriorsample = sample(sim1model1discrete, Prior(), 16384)
sim1model1discretepriorsampledf = DataFrame(sim1model1discretepriorsample)
plotchains(sim1model1discretepriorsampledf)
sim1fit1discretepriorsample = samplerenewalequation_2sets(
    fseir, sim1model1discretepriorsampledf, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    timeperiods=timeperiods2,
)
sim1fit1kvdiscretepriorsample = keyvalues(sim1model1discretepriorsampledf, sim1fit1discretepriorsample)
#println(sim1fit1kvdiscretepriorsample)

sim1discretepriorsampleplot = let
    fig = Figure(; size=( 500, 350 ))

    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])

    plotrenewalequationsamples!(
        ga,
        simulation1dataset["cases_counterfactual"],
        simulation1dataset["cases_counterfactual"], 
        W_sim1_0, 
        simulation1dataset["Ns"], 
        sim1fit0discretepriorsample,
        fitws(simulation1dataset["cases_counterfactual"], simulation1dataset["Ns"], sim1fit0discretepriorsample); 
        betafunctions=[ beta1a, beta1bcounterfactual ], 
        betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
        infectiousduration=2.5,
        #rhoclip = 2.5,
        columntitles=[ 
            "Group 1", 
            "Group 2" 
        ],
        columntitlefontsize=10,
        xtitle="Time (days)",
    )

    plotrenewalequationsamples!(
        gb, simulation1dataset, W_sim1, sim1fit1discretepriorsample; 
        betafunctions=[ beta1a, beta1b ], 
        betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
        infectiousduration=2.5,
        #rhoclip = 2.5,
        columntitles=[ 
            "Group 1", 
            "Group 2" 
        ],
        columntitlefontsize=10,
        xtitle="Time (days)",
    )

    labelplots!([ "A", "B" ], [ ga, gb ])

    fig
end

safesave(plotsdir("sim1discretepriorsampleplot.pdf"), sim1discretepriorsampleplot)

sim1chain0discrete = loadanalysisdictsasdf("sim1model0discrete", 8, maxrounds, 100)
plotchains(sim1chain0discrete)
sim1fit0discrete = samplerenewalequation_2sets(
    fseir, sim1chain0discrete, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases_counterfactual"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    timeperiods=timeperiods2,
)
sim1fit0kvdiscrete = keyvalues(sim1chain0discrete, sim1fit0discrete)
println(sim1fit0kvdiscrete)

sim1chain1discrete = loadanalysisdictsasdf("sim1model1discrete", 8, maxrounds, 110)
plotchains(sim1chain1discrete)
sim1fit1discrete = samplerenewalequation_2sets(
    fseir, sim1chain1discrete, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    timeperiods=timeperiods2,
)
sim1fit1kvdiscrete = keyvalues(sim1chain1discrete, sim1fit1discrete)
println(sim1fit1kvdiscrete)




#=
sim1fit1discreteplot = plotrenewalequationsamples(
    simulation1dataset["cases_counterfactual"],
    simulation1dataset["cases_counterfactual"], 
    W_sim1_0, 
    simulation1dataset["Ns"], 
    sim1fit0discrete,
    fitws(simulation1dataset["cases_counterfactual"], simulation1dataset["Ns"], sim1fit0discrete); 
    betafunctions=[ beta1a, beta1bcounterfactual ], 
    betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    #rhoclip = 2,
   # secondaryinterventions=[
   #     InterventionsMatrix([ nothing, 36 ], 100),
   #     InterventionsMatrix([ nothing, 64 ], 100)
   # ],
)
=#
sim1discreteplot = let
    fig = Figure(; size=( 500, 350 ))

    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])

    plotrenewalequationsamples!(
        ga,
        simulation1dataset["cases_counterfactual"],
        simulation1dataset["cases_counterfactual"], 
        W_sim1_0, 
        simulation1dataset["Ns"], 
        sim1fit0discrete,
        fitws(simulation1dataset["cases_counterfactual"], simulation1dataset["Ns"], sim1fit0discrete); 
        betafunctions=[ beta1a, beta1bcounterfactual ], 
        betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
        infectiousduration=2.5,
        rhoclip = 2.5,
        columntitles=[ 
            "Group 1", 
            "Group 2" 
        ],
        columntitlefontsize=10,
        xtitle="Time (days)",
    )

    plotrenewalequationsamples!(
        gb, simulation1dataset, W_sim1, sim1fit1discrete; 
        betafunctions=[ beta1a, beta1b ], 
        betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
        infectiousduration=2.5,
        rhoclip = 2.5,
        columntitles=[ 
            "Group 1", 
            "Group 2" 
        ],
        columntitlefontsize=10,
        xtitle="Time (days)",
    )

    labelplots!([ "A", "B" ], [ ga, gb ])

    fig
end

safesave(plotsdir("sim1discreteplot.pdf"), sim1discreteplot)


sim1chain0 = loadanalysisdictsasdf("sim1model0", 8, maxrounds, 100)
chainsplot1_0 = plotchains(
    sim1chain0; 
    size=( 500, 700 ),
    ylabels = [
        "log density",
        "ln ζ mean",
        "σ ζ",
        "ln ζ₁",
        "ln ζ₂",
        "ln η mean",
        "σ η",
        "ln η₁",
        "ln η₃",
        "ln η₄",
        "ln η₅",
        "ln η₆",
        "ln δ",
        "σ²",
        "θ",
    ],
)
safesave(plotsdir("chainsplot1_0.pdf"), chainsplot1_0)

sim1fit0 = samplerenewalequation_2sets(
    fseir, sim1chain0, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases_counterfactual"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim1fit0kv = keyvalues(sim1chain0, sim1fit0)
println(sim1fit0kv)

println("mean=$(exp(mean(sim1chain0.logdelta))); CrI=$(exp.(quantile(sim1chain0.logdelta, [ 0.05, 0.95 ]))))")

sim1chain1 = loadanalysisdictsasdf("sim1model1", 8, maxrounds, 110)
plotchains(sim1chain1)
sim1fit1 = samplerenewalequation_2sets(
    fseir, sim1chain1, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim1fit1kv = keyvalues(sim1chain1, sim1fit1)
println(sim1fit1kv)

println("mean=$(exp(mean(sim1chain1.logdelta))); CrI=$(exp.(quantile(sim1chain1.logdelta, [ 0.05, 0.95 ]))))")

sim1chain2 = loadanalysisdictsasdf("sim1model2", 8, maxrounds, 120)
plotchains(sim1chain2)
sim1fit2 = samplerenewalequation_2sets(
    fseir, sim1chain2, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
    secondaryinterventions=[
        InterventionsMatrix([ nothing, 36 ], 100),
        InterventionsMatrix([ nothing, 64 ], 100)
    ],
)
sim1fit2kv = keyvalues(sim1chain2, sim1fit2)
sim1fit1plot = plotrenewalequationsamples(
    simulation1dataset, W_sim1, sim1fit2; 
    betafunctions=[ beta1a, beta1b ], 
    betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    #rhoclip = 2,
   # secondaryinterventions=[
   #     InterventionsMatrix([ nothing, 36 ], 100),
   #     InterventionsMatrix([ nothing, 64 ], 100)
   # ],
)

sim1chaina1 = loadanalysisdictsasdf("sim1modela1", 8, maxrounds, 10110)
plotchains(sim1chaina1)
sim1fita1 = samplerenewalequation_2sets(
    fseir, sim1chaina1, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim1fita1kv = keyvalues(sim1chaina1, sim1fita1)
sim1fit1plot = plotrenewalequationsamples(
    simulation1dataset, W_sim1, sim1fita1; 
    betafunctions=[ beta1a, beta1b ], 
    betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    #rhoclip = 2,
   # secondaryinterventions=[
   #     InterventionsMatrix([ nothing, 36 ], 100),
   #     InterventionsMatrix([ nothing, 64 ], 100)
   # ],
)

sim1plot = let
    fig = Figure(; size=( 500, 350 ))

    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])

    plotrenewalequationsamples!(
        ga,
        simulation1dataset["cases_counterfactual"],
        simulation1dataset["cases_counterfactual"], 
        W_sim1_0, 
        simulation1dataset["Ns"], 
        sim1fit0,
        fitws(simulation1dataset["cases_counterfactual"], simulation1dataset["Ns"], sim1fit0); 
        betafunctions=[ beta1a, beta1bcounterfactual ], 
        betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
        infectiousduration=2.5,
        rhoclip = 2.5,
        columntitles=[ 
            "Group 1", 
            "Group 2" 
        ],
        columntitlefontsize=10,
        xtitle="Time (days)",
    )

    plotrenewalequationsamples!(
        gb, simulation1dataset, W_sim1, sim1fit1; 
        betafunctions=[ beta1a, beta1b ], 
        betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
        infectiousduration=2.5,
        rhoclip = 2.5,
        columntitles=[ 
            "Group 1", 
            "Group 2" 
        ],
        columntitlefontsize=10,
        xtitle="Time (days)",
    )

    labelplots!([ "A", "B" ], [ ga, gb ])

    fig
end

safesave(plotsdir("sim1plot.pdf"), sim1plot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 1a 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim1achain1 = loadanalysisdictsasdf("sim1amodel1", 8, maxrounds, 160)
plotchains(sim1achain1)
sim1afit1 = samplerenewalequation_2sets(
    fseir, sim1achain1, simulation1a_dataset["interventions"]; 
    initialvalues=simulation1a_dataset["cases"][1:10, :], 
    Ns=simulation1a_dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim1afit1kv = keyvalues(sim1achain1, sim1afit1)
println(sim1afit1kv)

sim1afit1plot = plotrenewalequationsamples(
    simulation1a_dataset, W_sim1a, sim1afit1; 
    betafunctions=[ beta1a_a, beta1a_b ], 
    betafunctions_counterfactual=[ beta1a_a, beta1a_bcounterfactual ],
    infectiousduration=2.5,
    plotsize=( 500, 350 ),
    rhoclip = 2.5,
   # secondaryinterventions=[
   #     InterventionsMatrix([ nothing, 36 ], 100),
   #     InterventionsMatrix([ nothing, 64 ], 100)
   # ],
)

safesave(plotsdir("sim1afit1plot.pdf"), sim1afit1plot)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 2 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim2chain0 = loadanalysisdictsasdf("sim2model0", 8, maxrounds, 200)
plotchains(sim2chain0)
sim2fit0 = samplerenewalequation_2sets(
    fseir, sim2chain0, simulation2dataset["interventions"]; 
    initialvalues=simulation2dataset["cases"][1:20, :], 
    Ns=simulation2dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim2fit0plot = plotrenewalequationsamples(
    simulation2dataset, W_sim2, sim2fit0; 
    betafunctions=[ beta2a, beta2b, beta2c ], 
    betafunctions_counterfactual=[ beta2a, beta2bcounterfactual, beta2ccounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    rhoclip=2.5,
)
sim2fit0kv = keyvalues(sim2chain0, sim2fit0)
println(sim2fit0kv)

sim2chain1 = loadanalysisdictsasdf("sim2model1", 8, maxrounds, 210)
plotchains(sim2chain1)
sim2fit1 = samplerenewalequation_2sets(
    fseir, sim2chain1, simulation2dataset["interventions"]; 
    initialvalues=simulation2dataset["cases"][1:20, :], 
    Ns=simulation2dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
    secondaryinterventions=InterventionsMatrix([ nothing, nothing, 70 ], 100),
)
sim2fit1plot = plotrenewalequationsamples(
    simulation2dataset, W_sim2, sim2fit1; 
    betafunctions=[ beta2a, beta2b, beta2c ], 
    betafunctions_counterfactual=[ beta2a, beta2bcounterfactual, beta2ccounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    rhoclip=2.5,
)
sim2fit1kv = keyvalues(sim2chain1, sim2fit1)
println(sim2fit1kv)

sim2plot = let
    fig = Figure(; size=( 500, 350 ))

    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])

    plotrenewalequationsamples!(
        ga, simulation2dataset, W_sim2, sim2fit0; 
        betafunctions=[ beta2a, beta2b, beta2c ], 
        betafunctions_counterfactual=[ beta2a, beta2bcounterfactual, beta2ccounterfactual ],
        infectiousduration=2.5,
        rhoclip = 2.5,
        columntitles=[ 
            "Group 1", 
            "Group 2",
            "Group 3" 
        ],
        columntitlefontsize=10,
        xtitle="Time (days)",
    )

    plotrenewalequationsamples!(
        gb, simulation2dataset, W_sim2, sim2fit1; 
        betafunctions=[ beta2a, beta2b, beta2c ], 
        betafunctions_counterfactual=[ beta2a, beta2bcounterfactual, beta2ccounterfactual ],
        infectiousduration=2.5,
        rhoclip = 2.5,
        columntitles=[ 
            "Group 1", 
            "Group 2",
            "Group 3" 
        ],
        columntitlefontsize=10,
        xtitle="Time (days)",
    )

    labelplots!([ "A", "B" ], [ ga, gb ])

    fig
end

safesave(plotsdir("sim2plot.pdf"), sim2plot)

sim2chain2 = loadanalysisdictsasdf("sim2model2", 8, maxrounds, 220)
plotchains(sim2chain2)
sim2fit2 = samplerenewalequation_2sets(
    fseir, sim2chain2, simulation2dataset["interventions"]; 
    initialvalues=simulation2dataset["cases"][1:20, :], 
    Ns=simulation2dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
    secondaryinterventions=[
        InterventionsMatrix([ nothing, 36, nothing ], 100),
        InterventionsMatrix([ nothing, nothing, 56 ], 100),
        InterventionsMatrix([ nothing, 64, nothing ], 100),
        InterventionsMatrix([ nothing, nothing, 84 ], 100),
    ],
)
sim2fit2plot = plotrenewalequationsamples(
    simulation2dataset, W_sim2, sim2fit2; 
    betafunctions=[ beta2a, beta2b, beta2c ], 
    betafunctions_counterfactual=[ beta2a, beta2bcounterfactual, beta2ccounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    rhoclip=2.5,
)

sim2chain3 = loadanalysisdictsasdf("sim2model3", 8, maxrounds, 230)
plotchains(sim2chain3)
sim2fit3 = samplerenewalequation_2sets(
    fseir, sim2chain3, simulation2dataset["interventions"]; 
    initialvalues=simulation2dataset["cases"][1:20, :], 
    Ns=simulation2dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
    secondaryinterventions=[
        InterventionsMatrix([ nothing, nothing, 30 ], 100),
        InterventionsMatrix([ nothing, 36, 56 ], 100),
        InterventionsMatrix([ nothing, 64, 84 ], 100)
    ],
)
sim2fit3plot = plotrenewalequationsamples(
    simulation2dataset, W_sim2, sim2fit3; 
    betafunctions=[ beta2a, beta2b, beta2c ], 
    betafunctions_counterfactual=[ beta2a, beta2bcounterfactual, beta2ccounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    rhoclip=2.5,
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 3 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim3chain0 = loadanalysisdictsasdf("sim3model0", 8, maxrounds, 300)
plotchains(sim3chain0)
sim3fit0 = samplerenewalequation_2sets(
    fseir, sim3chain0, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases"][1:20, :], 
    Ns=simulation3dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim3fit0plot = plotrenewalequationsamples(
    simulation3dataset["cases_counterfactual"],
    simulation3dataset["cases_counterfactual"], 
    W_sim3_0, 
    simulation3dataset["Ns"], 
    sim3fit0,
    fitws(simulation3dataset["cases_counterfactual"], simulation3dataset["Ns"], sim3fit0); 
    betafunctions=[ beta3a, beta3bcounterfactual ], 
    betafunctions_counterfactual=[ beta3a, beta3bcounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    rhoclip=2.5,
)

sim3chain1 = loadanalysisdictsasdf("sim3model1", 8, maxrounds, 310)
plotchains(sim3chain1)
sim3fit1 = samplerenewalequation_2sets(
    fseir, sim3chain1, simulation3dataset["interventions"]; 
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

sim3chain2 = loadanalysisdictsasdf("sim3model2", 8, maxrounds, 320)
plotchains(sim3chain2)
sim3fit2 = samplerenewalequation_2sets(
    fseir, sim3chain2, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases"][1:20, :], 
    Ns=simulation3dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
    secondaryinterventions=[
        InterventionsMatrix([ nothing, 36 ], 100),
        InterventionsMatrix([ nothing, 64 ], 100)
    ],
)
sim3fit2plot = plotrenewalequationsamples(
    simulation3dataset, W_sim3, sim3fit2; 
    betafunctions=[ beta3a, beta3b ], 
    betafunctions_counterfactual=[ beta3a, beta3bcounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    rhoclip=5,
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 4 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim4chain0 = loadanalysisdictsasdf("sim4model0", 8, maxrounds, 400)
plotchains(sim4chain0)
sim4fit0 = samplerenewalequation_2sets(
    fseir, sim4chain0, simulation4dataset["interventions"]; 
    initialvalues=simulation4dataset["cases"][1:20, :], 
    Ns=simulation4dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim4fit0plot = plotrenewalequationsamples(
    simulation4dataset["cases_counterfactual"],
    simulation4dataset["cases_counterfactual"], 
    W_sim4_0, 
    simulation4dataset["Ns"], 
    sim4fit0,
    fitws(simulation4dataset["cases_counterfactual"], simulation4dataset["Ns"], sim4fit0); 
    betafunctions=[ beta4a, beta4bcounterfactual ], 
    betafunctions_counterfactual=[ beta4a, beta4bcounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    rhoclip=2.5,
)

sim4chain1 = loadanalysisdictsasdf("sim4model1", 8, maxrounds, 410)
plotchains(sim4chain1)
sim4fit1 = samplerenewalequation_2sets(
    fseir, sim4chain1, simulation4dataset["interventions"]; 
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

sim4chain2 = loadanalysisdictsasdf("sim4model2", 8, maxrounds, 420)
plotchains(sim4chain2)
sim4fit2 = samplerenewalequation_2sets(
    fseir, sim4chain2, simulation4dataset["interventions"]; 
    initialvalues=simulation4dataset["cases"][1:20, :], 
    Ns=simulation4dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
    secondaryinterventions=[
        InterventionsMatrix([ nothing, 36 ], 100),
        InterventionsMatrix([ nothing, 64 ], 100)
    ],
)
sim4fit2plot = plotrenewalequationsamples(
    simulation4dataset, W_sim4, sim4fit2; 
    betafunctions=[ beta4a, beta4b ], 
    betafunctions_counterfactual=[ beta4a, beta4bcounterfactual ],
    infectiousduration=2.5,
    plotsize=( 400, 400 ),
    rhoclip=5,
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Covid data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

datachain1 = loadanalysisdictsasdf("datamodel1", 8, maxrounds, 1010)
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
    plotsize=( 500, 350 ),
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

safesave(plotsdir("datafit1plot.pdf"), datafit1plot)

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

## Analysis 2 
# Pillar 1 test results 

datachain2 = loadanalysisdictsasdf("datamodel2", 8, maxrounds, 1020)
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
    pil1covidcases, W_pil1coviddata, selectpops, datafit2;
    #plotsize=( 600, 400 ),
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


## Analysis 3

# with lead and lag 

datachain3 = loadanalysisdictsasdf("datamodel3", 8, maxrounds, 1030)
plotchains(datachain3)
datafit3 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain3, masstesting; 
    initialvalues=allcovidcases[1:100, :], 
    Ns=selectpops,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=collect(1.0:28:216),
    secondaryinterventions=[ 
        InterventionsMatrix([ 172, 172, 146, 172, 172, 217, 217, 217, 172 ], 216), 
        InterventionsMatrix([ 200, 200, 174, 200, 200, 217, 217, 217, 200 ], 216), 
    ],
)
datafit3plot = plotrenewalequationsamples(
    allcovidcases, W_allcoviddata, selectpops, datafit3;
    #plotsize=( 600, 400 ),
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

safesave(plotsdir("datafit3plot.svg"), datafit3plot)

#=
## Analysis 4

# with separate lead and lag for each area 

datachain4 = loadanalysisdictsasdf("datamodel4", 8, maxrounds, 1040)
plotchains(datachain4)
datafit4 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain4, masstesting; 
    initialvalues=allcovidcases[1:100, :], 
    Ns=selectpops,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=collect(1.0:28:216),
    secondaryinterventions=[ 
        InterventionsMatrix([ 172, 217, 217, 217, 217, 217, 217, 217, 217 ], 216), 
        InterventionsMatrix([ 200, 217, 217, 217, 217, 217, 217, 217, 217 ], 216), 
        InterventionsMatrix([ 217, 172, 217, 217, 217, 217, 217, 217, 217 ], 216), 
        InterventionsMatrix([ 217, 200, 217, 217, 217, 217, 217, 217, 217 ], 216),
        InterventionsMatrix([ 217, 217, 146, 217, 217, 217, 217, 217, 217 ], 216), 
        InterventionsMatrix([ 217, 217, 174, 217, 217, 217, 217, 217, 217 ], 216),        
        InterventionsMatrix([ 217, 217, 217, 172, 217, 217, 217, 217, 217 ], 216), 
        InterventionsMatrix([ 217, 217, 217, 200, 172, 217, 217, 217, 217 ], 216),    
        InterventionsMatrix([ 217, 217, 217, 217, 200, 217, 217, 217, 217 ], 216), 
        InterventionsMatrix([ 217, 217, 217, 217, 217, 217, 217, 217, 217 ], 216),    
        InterventionsMatrix([ 217, 217, 217, 217, 217, 217, 217, 217, 172 ], 216), 
        InterventionsMatrix([ 217, 217, 217, 217, 217, 217, 217, 217, 200 ], 216),         
    ],
)
datafit4plot = plotrenewalequationsamples(
    allcovidcases, W_allcoviddata, selectpops, datafit4;
    #plotsize=( 600, 400 ),
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
=#


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Masking data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Analysis 1 
# Effect of mask recommendations. No other considerations of confounding 

maskingdatachain1 = loadanalysisdictsasdf("maskingdatamodel1", 8, maxrounds, 1110)
plotchains(maskingdatachain1)
maskingdatafit1 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, maskingdatachain1, facialcoveringsrecommended[1:133, :]; 
    initialvalues=maskcovidcases[1:28, :], 
    Ns=POPULATION2020,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=collect(1.0:28:133),
)
maskingdatafit1plot = plotrenewalequationsamples(
    maskcovidcases[1:133, :], W_maskcoviddata[1:133, :], POPULATION2020, maskingdatafit1;
    #plotsize=( 600, 400 ),
    #=columntitles=[ 
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
    xtitle="Date, 2020–2021", =#
)

safesave(plotsdir("maskingdatafit1plot.pdf"), maskingdatafit1plot)




