
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
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 250 ))
        ax = Axis(fig[1, 1])
        scatter!(ax, 0:20, [ fseir(t) for t ∈ 0:20 ]; color=:black)
        formataxis!(ax)
    
        Label(
            fig[1, 0], L"generation interval, $g(\tau)$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(fig[2, 1], L"time, $\tau$"; fontsize=11.84, tellwidth=false)
    
        colgap!(fig.layout, 1, 5)
        rowgap!(fig.layout, 1, 5)
        fig
    end
    fig
end

safesave(plotsdir("generationintervalplot.pdf"), generationintervalplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Discrete-time analyses

### Priors 
#=
sim1model0discretepriorsample = sample(sim1model0discrete, Prior(), 16384)
sim1model0discretepriorsampledf = DataFrame(sim1model0discretepriorsample)
plotchains(sim1model0discretepriorsampledf)
sim1fit0discretepriorsample = samplerenewalequation_2sets(
    fseir, sim1model0discretepriorsampledf, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases_counterfactual"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    timeperiods=timeperiods2,
)
sim1fit0kvdiscretepriorsample = keyvalues(
    sim1model0discretepriorsampledf, sim1fit0discretepriorsample
)
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
sim1fit1kvdiscretepriorsample = keyvalues(
    sim1model1discretepriorsampledf, sim1fit1discretepriorsample
)
println(sim1fit1kvdiscretepriorsample)
=#
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

sim1fit0 = samplerenewalequation_2sets(
    fseir, sim1chain0, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases_counterfactual"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim1fit0kv = keyvalues(sim1chain0, sim1fit0)
println(sim1fit0kv)

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

sim1chaina1 = loadanalysisdictsasdf("sim1modela1", 8, maxrounds, 10110)
plotchains(sim1chaina1)
sim1fita1 = samplerenewalequation_2sets(
    fseir, sim1chaina1, simulation1dataset["interventions"]; 
    initialvalues=simulation1dataset["cases"][1:10, :], 
    Ns=simulation1dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim1fita1kv = keyvalues(sim1chaina1, sim1fita1)

sim1plot = let
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 587, 700 ))

        ga = GridLayout(fig[1, 1])
        gb = GridLayout(fig[1, 2])        
        gc = GridLayout(fig[2, 1])
        gd = GridLayout(fig[2, 2])

        plotrenewalequationsamples!(
            ga,
            simulation1dataset["cases_counterfactual"],
            simulation1dataset["cases_counterfactual"], 
            W_sim1_0, 
            simulation1dataset["Ns"], 
            sim1fit0discrete,
            fitws(
                simulation1dataset["cases_counterfactual"], 
                simulation1dataset["Ns"], 
                sim1fit0discrete
            ); 
            betafunctions=[ beta1a, beta1bcounterfactual ], 
            betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
            infectiousduration=2.5,
            rhoclip = 2.5,
            simcolour=COLOURVECTOR[3],
            columntitles=[ "Group 1", "Group 2" ],
            columntitlefontsize=10,
            markersize=2,
            xtitle="Time, days",
        )

        plotrenewalequationsamples!(
            gb, simulation1dataset, W_sim1, sim1fit1discrete; 
            betafunctions=[ beta1a, beta1b ], 
            betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
            infectiousduration=2.5,
            rhoclip = 2.5,
            simcolour=COLOURVECTOR[3],
            columntitles=[ "Group 1", "Group 2" ],
            columntitlefontsize=10,
            markersize=2,
            xtitle="Time, days",
        )

        plotrenewalequationsamples!(
            gc,
            simulation1dataset["cases_counterfactual"],
            simulation1dataset["cases_counterfactual"], 
            W_sim1_0, 
            simulation1dataset["Ns"], 
            sim1fit0,
            fitws(
                simulation1dataset["cases_counterfactual"], 
                simulation1dataset["Ns"], 
                sim1fit0
            ); 
            betafunctions=[ beta1a, beta1bcounterfactual ], 
            betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
            infectiousduration=2.5,
            rhoclip = 2.5,
            simcolour=COLOURVECTOR[3],
            columntitles=[ "Group 1", "Group 2" ],
            columntitlefontsize=10,
            markersize=2,
            xtitle="Time, days",
        )

        plotrenewalequationsamples!(
            gd, simulation1dataset, W_sim1, sim1fit1; 
            betafunctions=[ beta1a, beta1b ], 
            betafunctions_counterfactual=[ beta1a, beta1bcounterfactual ],
            infectiousduration=2.5,
            rhoclip = 2.5,
            simcolour=COLOURVECTOR[3],
            columntitles=[ "Group 1", "Group 2" ],
            columntitlefontsize=10,
            markersize=2,
            xtitle="Time, days",
        )

        labelplots!([ "A", "B", "C", "D" ], [ ga, gb, gc, gd ])
        fig
    end
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

sim4chain1 = loadanalysisdictsasdf("sim4model1", 8, maxrounds, 410)
plotchains(sim4chain1)
sim4fit1 = samplerenewalequation_2sets(
    fseir, sim4chain1, simulation4dataset["interventions"]; 
    initialvalues=simulation4dataset["cases"][1:20, :], 
    Ns=simulation4dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)

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

sim2plot = let
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 587, 700 ))

        gu = GridLayout(fig[1, 1]) 
        gl = GridLayout(fig[2, 1])
        ga = GridLayout(gu[1, 1])
        gb = GridLayout(gu[1, 2])        
        gc = GridLayout(gl[1, 1])
        gd = GridLayout(gl[1, 2])

        plotrenewalequationsamples!(
            ga, simulation1a_dataset, W_sim1a, sim1afit1; 
            betafunctions=[ beta1a, beta1b ], 
            betafunctions_counterfactual=[ beta1a_a, beta1a_bcounterfactual ],
            infectiousduration=2.5,
            #rhoclip = 2.5,
            simcolour=COLOURVECTOR[3],
            columntitles=[ "Group 1", "Group 2" ],
            columntitlefontsize=10,
            markersize=2,
            xtitle="Time, days",
        )

        plotrenewalequationsamples!(
            gb, simulation2dataset, W_sim2, sim2fit2; 
            betafunctions=[ beta2a, beta2b, beta2c ], 
            betafunctions_counterfactual=[ beta2a, beta2bcounterfactual, beta2ccounterfactual ],
            infectiousduration=2.5,
            #rhoclip = 2.5,
            simcolour=COLOURVECTOR[3],
            columntitles=[ "Group 1", "Group 2", "Group 3" ],
            columntitlefontsize=10,
            markersize=2,
            xtitle="Time, days",
        )

        plotrenewalequationsamples!(
            gc, simulation3dataset, W_sim3, sim3fit2; 
            betafunctions=[ beta3a, beta3b ], 
            betafunctions_counterfactual=[ beta3a, beta3bcounterfactual ],
            infectiousduration=2.5,
            #rhoclip = 2.5,
            simcolour=COLOURVECTOR[3],
            columntitles=[ "Group 1", "Group 2" ],
            columntitlefontsize=10,
            markersize=2,
            xtitle="Time, days",
        )

        plotrenewalequationsamples!(
            gd, simulation4dataset, W_sim4, sim4fit2; 
            betafunctions=[ beta4a, beta4b ], 
            betafunctions_counterfactual=[ beta4a, beta4bcounterfactual ],        
            infectiousduration=2.5,
            #rhoclip = 2.5,
            simcolour=COLOURVECTOR[3],
            columntitles=[ "Group 1", "Group 2" ],
            columntitlefontsize=10,
            markersize=2,
            xtitle="Time, days",
        )

        colsize!(gu, 1, Auto(2))
        colsize!(gu, 2, Auto(3))
        labelplots!([ "A", "B", "C", "D" ], [ ga, gb, gc, gd ])
        fig
    end
    fig
end

safesave(plotsdir("sim2plot.pdf"), sim2plot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Covid data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

datachain1 = loadanalysisdictsasdf("datamodel1", 8, maxrounds, 1010)
plotchains(datachain1)
datafit1 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain1, masstesting; 
    initialvalues=allcovidcases[1:100, :], 
    Ns=selectpops,
    timeknots=[ collect(1.0:28:216); [216] ],
)
datafit1plot = plotrenewalequationsamples(
    allcovidcases, W_allcoviddata, selectpops, datafit1;
    plotsize=( 587, 411 ),
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
    markersize=2,
    xticklabelrotation=-π/4,
    xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
    xtitle="Date, 2020–2021",
)

safesave(plotsdir("datafit1plot.pdf"), datafit1plot)

datafit1kv = keyvalues(datachain1, datafit1)
println(datafit1kv)


## Analysis 2 
# Pillar 1 test results 

datachain2 = loadanalysisdictsasdf("datamodel2", 8, maxrounds, 1020)
plotchains(datachain2)
datafit2 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain2, masstesting; 
    initialvalues=pil1covidcases[1:124, :], 
    Ns=selectpops,
    timeknots=[ collect(1.0:28:216); [216] ],
)
datafit2plot = plotrenewalequationsamples(
    pil1covidcases, W_pil1coviddata, selectpops, datafit2;
    plotsize=( 587, 411 ),
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
    markersize=2,
    xticklabelrotation=-π/4,
    xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
    xtitle="Date, 2020–2021",
)

datafit2kv = keyvalues(datachain2, datafit2)
println(datafit2kv)


## Analysis 3

# with lead and lag 

datachain3 = loadanalysisdictsasdf("datamodel3", 8, maxrounds, 1030)
plotchains(datachain3)
datafit3 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain3, masstesting; 
    initialvalues=allcovidcases[1:100, :], 
    Ns=selectpops,
    timeknots=[ collect(1.0:28:216); [216] ],
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


## Analysis 4

# with separate lead and lag for each area 

datachain4 = loadanalysisdictsasdf("datamodel4", 8, maxrounds, 1040)
plotchains(datachain4)
datafit4 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain4, masstesting; 
    initialvalues=allcovidcases[1:100, :], 
    Ns=selectpops,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=[ collect(1.0:28:216); [216] ],
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

## Analysis 5 
# Pillar 1 test results with lead and lag

datachain5 = loadanalysisdictsasdf("datamodel5", 8, maxrounds, 1050)
plotchains(datachain5)
datafit2 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain5, masstesting; 
    initialvalues=pil1covidcases[1:124, :], 
    Ns=selectpops,
    timeknots=[ collect(1.0:28:216); [216] ],
    secondaryinterventions=[ 
        InterventionsMatrix([ 172, 172, 146, 172, 172, 217, 217, 217, 172 ], 216), 
        InterventionsMatrix([ 200, 200, 174, 200, 200, 217, 217, 217, 200 ], 216), 
    ],
)
datafit5plot = plotrenewalequationsamples(
    pil1covidcases, W_pil1coviddata, selectpops, datafit5;
    plotsize=( 587, 411 ),
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
    markersize=2,
    xticklabelrotation=-π/4,
    xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
    xtitle="Date, 2020–2021",
)

datafit5kv = keyvalues(datachain5, datafit5)
println(datafit5kv)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Masking data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Analysis 1 
# Effect of mask recommendations. No other considerations of confounding 

maskingdatachain1 = loadanalysisdictsasdf("maskingdatamodel1", 8, maxrounds, 1110)
plotchains(maskingdatachain1)
maskingdatafit1 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, maskingdatachain1, facialcoveringsrecommended[1:191, :]; 
    initialvalues=maskcovidcases[1:90, :], 
    Ns=POPULATION2020,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=[ 1.0; collect(56.0:28:191); 191 ],
)
maskingdatafit1plot = plotrenewalequationsamples(
    maskcovidcases[1:191, :], W_maskcoviddata[1:191, :], POPULATION2020, maskingdatafit1;
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

## Analysis 2 
# Effect of mask requirements. No other considerations of confounding 

maskingdatachain2 = loadanalysisdictsasdf("maskingdatamodel2", 8, maxrounds, 1120)
plotchains(maskingdatachain2)
maskingdatafit2 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, maskingdatachain2, facialcoveringsrequired; 
    initialvalues=maskcovidcases[1:90, :], 
    Ns=POPULATION2020,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=[ 1.0; collect(56.0:28:224); 257 ],
)
maskingdatafit2plot = plotrenewalequationsamples(
    maskcovidcases, W_maskcoviddata, POPULATION2020, maskingdatafit2;
)


## Analysis 3 
# Effect of mask requirements with secondary interventions of end of stay-at-home and some
# businesses reopening

maskingdatachain3 = loadanalysisdictsasdf("maskingdatamodel3", 8, maxrounds, 1130)
plotchains(maskingdatachain3)
maskingdatafit3 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, maskingdatachain3, facialcoveringsrequired; 
    initialvalues=maskcovidcases[1:100, :], 
    Ns=POPULATION2020,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=[ 1.0; collect(56.0:28:224); 257 ],
    secondaryinterventions=[ endstayathometimes, somebusinessreopen ],
)
maskingdatafit3plot = plotrenewalequationsamples(
    maskcovidcases, W_maskcoviddata, POPULATION2020, maskingdatafit3;
)


## Analysis 4 
# Add lead and lag times  

maskingdatachain4 = loadanalysisdictsasdf("maskingdatamodel4", 8, maxrounds, 1140)
plotchains(maskingdatachain4)
maskingdatafit4 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, maskingdatachain4, facialcoveringsrequired; 
    initialvalues=maskcovidcases[1:100, :], 
    Ns=POPULATION2020,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=[ 1.0; collect(56.0:28:224); 257 ],
    secondaryinterventions=[ 
        [ endstayathometimes, somebusinessreopen ]; secondaryinterventions_data 
    ],
)
maskingdatafit4plot = plotrenewalequationsamples(
    maskcovidcases, W_maskcoviddata, POPULATION2020, maskingdatafit4;
)



