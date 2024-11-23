
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
        formataxis!(ax; trimspines=true, hidespines=( :t, :r ))
    
        Label(
            fig[1, 0], L"generation interval, $g(x)$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(fig[2, 1], L"time, $x$"; fontsize=11.84, tellwidth=false)
    
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

subsetsim1plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    g0 = GridLayout(fig[1, 0])
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    
    let
        axs1 = plotrenewalequationsamples_w!(
            ga, 
            simulation1dataset["cases_counterfactual"], 
            W_sim1_0, sim1fit0discrete, 
            fitws(
                simulation1dataset["cases_counterfactual"], 
                simulation1dataset["Ns"], 
                sim1fit0discrete
            ), 
            1;
            markersize=2,
            hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
        )
        axs2 = plotrenewalequationsamples_r0!(
            ga, simulation1dataset["cases_counterfactual"], sim1fit0discrete, 2;
            betafunctions=[ beta1a, beta1bcounterfactual ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            yticks=[ 1.4, 1.8, 2.2 ], ytitle=L"$\mathcal{R}_0$",
        )
        setvalue!(axs2[1], 1.4)
        setvalue!(axs2[1], 2.2)
        axs3 = plotrenewalequationsamples_cases!(
            ga, 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["Ns"], 
            sim1fit0discrete, 
            3;
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle="No\nintervetion",
        )
        axs4 = plotrenewalequationsamples_cases!(
            ga, 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["Ns"],
            sim1fit0discrete,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle="Intervention",
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            ga, simulation1dataset["cases_counterfactual"], simulation1dataset["cases_counterfactual"], simulation1dataset["Ns"], sim1fit0discrete, 5;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle="Cumulative\ndifference",
        )
        for axs ∈ [ axs1, axs2, axs3, axs4, axs5 ]
            if axs === axs5 
                formataxis!(
                    axs[1]; 
                    hidespines=( :r, :t ), trimspines=true,
                )
                formataxis!(
                    axs[2]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            else
                formataxis!(
                    axs[1]; 
                    hidex=true, hidexticks=true, 
                    hidespines=( :r, :t, :b ), trimspines=true,
                )
                formataxis!(
                    axs[2]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    
        interventionax = Axis(ga[1:5, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )

        linkaxes!(axs3..., axs4...)
        linkxaxes!(axs1[2], axs2[2], axs3[2], axs4[2], axs5[2], interventionax)

        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
   
        for r ∈ [ 1, 6 ] rowgap!(ga, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_w!(
            gb, 
            simulation1dataset["cases"], 
            W_sim1, 
            sim1fit1discrete, 
            fitws(
                simulation1dataset["cases"], 
                simulation1dataset["Ns"], 
                sim1fit1discrete
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=nothing,
            yticks=[ -1, 0, 1, 2 ], 
        )
        setvalue!(axs1[1], -1)
        axs2 = plotrenewalequationsamples_r0!(
            gb, simulation1dataset["cases"], sim1fit1discrete, 2;
            betafunctions=[ beta1a, beta1b ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=nothing,
            yticks=[ 1.4, 1.8, 2.2 ], 
        )
        setvalue!(axs2[1], 1.4)
        setvalue!(axs2[1], 2.2)
        axs3 = plotrenewalequationsamples_cases!(
            gb, 
            simulation1dataset["cases"], 
            simulation1dataset["Ns"], 
            sim1fit1discrete, 
            3;
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=nothing,
        )
        axs4 = plotrenewalequationsamples_cases!(
            gb, 
            simulation1dataset["cases"], 
            simulation1dataset["Ns"],
            sim1fit1discrete,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle=nothing,
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            gb, simulation1dataset["cases"], simulation1dataset["cases_counterfactual"], simulation1dataset["Ns"], sim1fit1discrete, 5;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle=nothing,
        )
        for axs ∈ [ axs1, axs2, axs3, axs4, axs5 ]
            if axs === axs5 
                formataxis!(
                    axs[1]; 
                    hidespines=( :r, :t ), trimspines=true,
                )
                formataxis!(
                    axs[2]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            else
                formataxis!(
                    axs[1]; 
                    hidex=true, hidexticks=true, 
                    hidespines=( :r, :t, :b ), trimspines=true,
                )
                formataxis!(
                    axs[2]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    
        interventionax = Axis(gb[1:5, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )

        linkaxes!(axs3..., axs4...)
        linkxaxes!(axs1[2], axs2[2], axs3[2], axs4[2], axs5[2], interventionax)

        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                gb[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
    
        for r ∈ [ 1, 6 ] rowgap!(gb, r, 5) end
    end

    Label(
        g0[3:5, 0], L"Diagnoses, per $100\,000$"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    colgap!(ga, 1, 5) 

    labelplots!([ "A", "B" ], [ ga, gb ]; cols=1,)

    colgap!(fig.layout, 1, -5)
    colsize!(fig.layout, 0, Auto(0.05))
    colsize!(fig.layout, 2, Auto(0.75))

    fig
end

safesave(plotsdir("subsetsim1plot.pdf"), subsetsim1plot)


## Continuous time analysis

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

subsetsim1plotc = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 587, 250 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    
    let
        axs1 = plotrenewalequationsamples_r0!(
            ga, simulation1dataset["cases_counterfactual"], sim1fit0, 1;
            betafunctions=[ beta1a, beta1bcounterfactual ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=L"$\mathcal{R}_0$",
        )
        axs2 = plotrenewalequationsamples_causaleffect!(
            ga, 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["Ns"], 
            sim1fit0, 
            2;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle=L"Cumulative \\ incidence$$",
        )
        
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
     
        colgap!(ga, 1, 5)  
        for r ∈ [ 1, 3 ] rowgap!(ga, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_r0!(
            gb, simulation1dataset["cases"], sim1fit1, 1;
            betafunctions=[ beta1a, beta1b ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=L"$\mathcal{R}_0$",
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            gb, 
            simulation1dataset["cases"], 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["Ns"], 
            sim1fit1,
            2;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle=L"Cumulative \\ incidence$$",
        )
        
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                gb[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
    
        colgap!(gb, 1, 5)  
        for r ∈ [ 1, 3 ] rowgap!(gb, r, 5) end
    end

    labelplots!([ "A", "B" ], [ ga, gb ])

    fig
end

safesave(plotsdir("subsetsim1plotc.pdf"), subsetsim1plotc)

subsetsim1plotcsuppl = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 587, 350 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    
    let
        axs1 = plotrenewalequationsamples_w!(
            ga, 
            simulation1dataset["cases_counterfactual"], 
            W_sim1_0, sim1fit0, 
            fitws(
                simulation1dataset["cases_counterfactual"], 
                simulation1dataset["Ns"], 
                sim1fit0
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=L"$w_{jt}$",
        )
        axs2 = plotrenewalequationsamples_cases!(
            ga, 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["Ns"], 
            sim1fit0, 
            2;
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=L"Without \\ intervetion$$",
        )
        axs3 = plotrenewalequationsamples_cases!(
            ga, 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["Ns"],
            sim1fit0,
            3;
            hidex=false,
            markersize=2, fittedparameter=:y_matrix_det_vec,
            xtitle="Time, days",
            ytitle=L"With \\ intervetion$$",
        )

    
        linkaxes!(axs2..., axs3...)
    
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
     
        colgap!(ga, 1, 5)  
        for r ∈ [ 1, 4 ] rowgap!(ga, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_w!(
            gb, 
            simulation1dataset["cases"], 
            W_sim1, 
            sim1fit1discrete, 
            fitws(
                simulation1dataset["cases"], 
                simulation1dataset["Ns"], 
                sim1fit1
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=L"$w_{jt}$",
        )
        axs2 = plotrenewalequationsamples_cases!(
            gb, 
            simulation1dataset["cases"], 
            simulation1dataset["Ns"], 
            sim1fit1, 
            2;
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=L"Without \\ intervetion$$",
        )
        axs3 = plotrenewalequationsamples_cases!(
            gb, 
            simulation1dataset["cases"], 
            simulation1dataset["Ns"],
            sim1fit1,
            3;
            hidex=false,
            markersize=2, fittedparameter=:y_matrix_det_vec,
            xtitle="Time, days",
            ytitle=L"With \\ intervetion$$",
        )

        linkaxes!(axs2..., axs3...)
    
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                gb[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
    
        colgap!(gb, 1, 5)  
        for r ∈ [ 1, 4 ] rowgap!(gb, r, 5) end    end

    labelplots!([ "A", "B" ], [ ga, gb ])

    fig
end

safesave(plotsdir("subsetsim1plotcsuppl.pdf"), subsetsim1plotcsuppl)


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

subsetsim1aplot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 250 ))
    ga = GridLayout(fig[1, 1])
    
    let
        axs1 = plotrenewalequationsamples_r0!(
            ga, simulation1a_dataset["cases_counterfactual"], sim1afit1, 1;
            betafunctions=[ beta1a_a, beta1a_b ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=L"$\mathcal{R}_0$",
            rhoclip=3,
        )
        axs2 = plotrenewalequationsamples_causaleffect!(
            ga, 
            simulation1a_dataset["cases"], 
            simulation1a_dataset["cases_counterfactual"], 
            simulation1a_dataset["Ns"], 
            sim1afit1, 
            2;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle=L"Cumulative \\ incidence$$",
        )
        
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
     
        colgap!(ga, 1, 5)  
        for r ∈ [ 1, 3 ] rowgap!(ga, r, 5) end
    end

    fig
end

safesave(plotsdir("subsetsim1aplot.pdf"), subsetsim1aplot)

subsetsim1aplotsuppl = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 350 ))
    ga = GridLayout(fig[1, 1])
    
    let
        axs1 = plotrenewalequationsamples_w!(
            ga, 
            simulation1a_dataset["cases"], 
            W_sim1a, sim1fit0, 
            fitws(
                simulation1a_dataset["cases"], 
                simulation1a_dataset["Ns"], 
                sim1afit1
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=L"$w_{jt}$",
        )
        axs2 = plotrenewalequationsamples_cases!(
            ga, 
            simulation1a_dataset["cases"], 
            simulation1a_dataset["Ns"], 
            sim1afit1, 
            2;
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=L"Without \\ intervetion$$",
        )
        axs3 = plotrenewalequationsamples_cases!(
            ga, 
            simulation1a_dataset["cases"], 
            simulation1a_dataset["Ns"],
            sim1afit1,
            3;
            hidex=false,
            markersize=2, fittedparameter=:y_matrix_det_vec,
            xtitle="Time, days",
            ytitle=L"With \\ intervetion$$",
        )

    
        linkaxes!(axs2..., axs3...)
    
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
     
        colgap!(ga, 1, 5)  
        for r ∈ [ 1, 4 ] rowgap!(ga, r, 5) end
    end

    fig
end

safesave(plotsdir("subsetsim1aplotsuppl.pdf"), subsetsim1aplotsuppl)

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
sim2fit1kv = keyvalues(sim2chain1, sim2fit1)
println(sim2fit1kv)

subsetsim2plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 590 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[2, 1])
    
    let
        axs1 = plotrenewalequationsamples_r0!(
            ga, simulation2dataset["cases"], sim2fit0, 1;
            betafunctions=[ beta2a, beta2b, beta2c ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=L"$\mathcal{R}_0$",
            rhoclip=3,
        )
        axs2 = plotrenewalequationsamples_causaleffect!(
            ga, 
            simulation2dataset["cases"], 
            simulation2dataset["cases_counterfactual"], 
            simulation2dataset["Ns"], 
            sim2fit0, 
            2;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle=L"Cumulative \\ incidence$$",
        )
        
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2", "Group 3" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
     
        colgap!(ga, 1, 5)  
        for r ∈ [ 1, 3 ] rowgap!(ga, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_r0!(
            gb, simulation2dataset["cases"], sim2fit1, 1;
            betafunctions=[ beta2a, beta2b, beta2c ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=L"$\mathcal{R}_0$",
            rhoclip=3,
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            gb, 
            simulation2dataset["cases"], 
            simulation2dataset["cases_counterfactual"], 
            simulation2dataset["Ns"],
            sim2fit1,
            2;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle=L"Cumulative \\ incidence$$",
        )
        
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2", "Group 3" ])
            Label(
                gb[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
    
        colgap!(gb, 1, 5)  
        for r ∈ [ 1, 3 ] rowgap!(gb, r, 5) end
    end

    labelplots!([ "A", "B" ], [ ga, gb ])
    rowgap!(fig.layout, 1, 5) 

    fig
end

safesave(plotsdir("subsetsim2plot.pdf"), subsetsim2plot)

subsetsim2plotsuppl = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 600 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[2, 1])
    
    let
        axs1 = plotrenewalequationsamples_w!(
            ga, 
            simulation2dataset["cases"], 
            W_sim2, sim2fit0, 
            fitws(
                simulation2dataset["cases"], 
                simulation2dataset["Ns"], 
                sim2fit0
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=L"$w_{jt}$",
        )
        axs2 = plotrenewalequationsamples_cases!(
            ga, 
            simulation2dataset["cases"], 
            simulation2dataset["Ns"], 
            sim2fit0, 
            2;
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=L"Without \\ intervetion$$",
        )
        axs3 = plotrenewalequationsamples_cases!(
            ga, 
            simulation2dataset["cases"], 
            simulation2dataset["Ns"],
            sim2fit0,
            3;
            hidex=false,
            markersize=2, fittedparameter=:y_matrix_det_vec,
            xtitle="Time, days",
            ytitle=L"With \\ intervetion$$",
        )

        linkaxes!(axs2..., axs3...)
    
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2", "Group 3" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
     
        colgap!(ga, 1, 5)  
        for r ∈ [ 1, 4 ] rowgap!(ga, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_w!(
            gb, 
            simulation2dataset["cases"], 
            W_sim2, 
            sim2fit1, 
            fitws(
                simulation2dataset["cases"], 
                simulation2dataset["Ns"], 
                sim2fit1
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=L"$w_{jt}$",
        )
        axs2 = plotrenewalequationsamples_cases!(
            gb, 
            simulation2dataset["cases"], 
            simulation2dataset["Ns"], 
            sim2fit1, 
            2;
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=L"Without \\ intervetion$$",
        )
        axs3 = plotrenewalequationsamples_cases!(
            gb, 
            simulation2dataset["cases"], 
            simulation2dataset["Ns"],
            sim2fit1,
            3;
            hidex=false,
            markersize=2, fittedparameter=:y_matrix_det_vec,
            xtitle="Time, days",
            ytitle=L"With \\ intervetion$$",
        )

        linkaxes!(axs2..., axs3...)
    
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2", "Group 3" ])
            Label(
                gb[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
    
        colgap!(gb, 1, 5)  
        for r ∈ [ 1, 4 ] rowgap!(gb, r, 5) end    end

    labelplots!([ "A", "B" ], [ ga, gb ])
    rowgap!(fig.layout, 1, 5) 

    fig
end

safesave(plotsdir("subsetsim2plotsuppl.pdf"), subsetsim2plotsuppl)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 3 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim3chain0 = loadanalysisdictsasdf("sim3model0", 8, maxrounds, 300)
plotchains(sim3chain0)
sim3fit0 = samplerenewalequation_2sets(
    fseir, sim3chain0, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases_counterfactual"][1:20, :], 
    Ns=simulation3dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)

sim3chain0leadlag = loadanalysisdictsasdf("sim3model0leadlag", 8, maxrounds, 305)
plotchains(sim3chain0leadlag)
sim3fit0leadlag = samplerenewalequation_2sets(
    fseir, sim3chain0leadlag, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases_counterfactual"][1:20, :], 
    Ns=simulation3dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
    secondaryinterventions=[
        InterventionsMatrix([ nothing, 36 ], 100),
        InterventionsMatrix([ nothing, 64 ], 100)
    ],
)

subsetsim3plot0 = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 587, 250 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    
    let
        axs1 = plotrenewalequationsamples_r0!(
            ga, simulation3dataset["cases_counterfactual"], sim3fit0, 1;
            betafunctions=[ beta3a, beta3bcounterfactual ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=L"$\mathcal{R}_0$",
            rhoclip=3,
        )
        axs2 = plotrenewalequationsamples_causaleffect!(
            ga, 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit0, 
            2;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle=L"Cumulative \\ incidence$$",
        )
        
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
     
        colgap!(ga, 1, 5)  
        for r ∈ [ 1, 3 ] rowgap!(ga, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_r0!(
            gb, simulation3dataset["cases_counterfactual"], sim3fit0leadlag, 1;
            betafunctions=[ beta3a, beta3bcounterfactual ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=L"$\mathcal{R}_0$",
            rhoclip=3,
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            gb, 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit0leadlag,
            2;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle=L"Cumulative \\ incidence$$",
        )
        
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                gb[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
    
        colgap!(gb, 1, 5)  
        for r ∈ [ 1, 3 ] rowgap!(gb, r, 5) end
    end

    labelplots!([ "A", "B" ], [ ga, gb ])

    fig
end

safesave(plotsdir("subsetsim3plot0.pdf"), subsetsim3plot0)

subsetsim3plot0suppl = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 587, 350 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    
    let
        axs1 = plotrenewalequationsamples_w!(
            ga, 
            simulation3dataset["cases_counterfactual"], 
            W_sim3_0, sim3fit0, 
            fitws(
                simulation3dataset["cases_counterfactual"], 
                simulation3dataset["Ns"], 
                sim3fit0
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=L"$w_{jt}$",
        )
        axs2 = plotrenewalequationsamples_cases!(
            ga, 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit0, 
            2;
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=L"Without \\ intervetion$$",
        )
        axs3 = plotrenewalequationsamples_cases!(
            ga, 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"],
            sim3fit0,
            3;
            hidex=false,
            markersize=2, fittedparameter=:y_matrix_det_vec,
            xtitle="Time, days",
            ytitle=L"With \\ intervetion$$",
        )

    
        linkaxes!(axs2..., axs3...)
    
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
     
        colgap!(ga, 1, 5)  
        for r ∈ [ 1, 4 ] rowgap!(ga, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_w!(
            gb, 
            simulation3dataset["cases_counterfactual"], 
            W_sim3_0, 
            sim3fit0leadlag, 
            fitws(
                simulation3dataset["cases_counterfactual"], 
                simulation3dataset["Ns"], 
                sim3fit0leadlag
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=L"$w_{jt}$",
        )
        axs2 = plotrenewalequationsamples_cases!(
            gb, 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit0leadlag, 
            2;
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=L"Without \\ intervetion$$",
        )
        axs3 = plotrenewalequationsamples_cases!(
            gb, 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"],
            sim3fit0leadlag,
            3;
            hidex=false,
            markersize=2, fittedparameter=:y_matrix_det_vec,
            xtitle="Time, days",
            ytitle=L"With \\ intervetion$$",
        )

        linkaxes!(axs2..., axs3...)
    
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                gb[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
    
        colgap!(gb, 1, 5)  
        for r ∈ [ 1, 4 ] rowgap!(gb, r, 5) end    end

    labelplots!([ "A", "B" ], [ ga, gb ])

    fig
end

safesave(plotsdir("subsetsim3plot0suppl.pdf"), subsetsim3plot0suppl)

## with intervention

sim3chain1 = loadanalysisdictsasdf("sim3model1", 8, maxrounds, 310)
plotchains(sim3chain1)
sim3fit1 = samplerenewalequation_2sets(
    fseir, sim3chain1, simulation3dataset["interventions"]; 
    initialvalues=simulation3dataset["cases"][1:20, :], 
    Ns=simulation3dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)

sim3fit1kv = keyvalues(sim3chain1, sim3fit1)
println(sim3fit1kv)

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

subsetsim3plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 587, 250 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    
    let
        axs1 = plotrenewalequationsamples_r0!(
            ga, simulation3dataset["cases"], sim3fit1, 1;
            betafunctions=[ beta3a, beta3b ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=L"$\mathcal{R}_0$",
            rhoclip=4.1,
        )
        axs2 = plotrenewalequationsamples_causaleffect!(
            ga, 
            simulation3dataset["cases"], 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit1, 
            2;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle=L"Cumulative \\ incidence$$",
        )
        
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
     
        colgap!(ga, 1, 5)  
        for r ∈ [ 1, 3 ] rowgap!(ga, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_r0!(
            gb, simulation3dataset["cases"], sim3fit2, 1;
            betafunctions=[ beta3a, beta3b ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=L"$\mathcal{R}_0$",
            rhoclip=4.1,
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            gb, 
            simulation3dataset["cases"], 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit2,
            2;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle=L"Cumulative \\ incidence$$",
        )
        
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                gb[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
    
        colgap!(gb, 1, 5)  
        for r ∈ [ 1, 3 ] rowgap!(gb, r, 5) end
    end

    labelplots!([ "A", "B" ], [ ga, gb ])

    fig
end

safesave(plotsdir("subsetsim3plot.pdf"), subsetsim3plot)

subsetsim3plotsuppl = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 587, 350 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    
    let
        axs1 = plotrenewalequationsamples_w!(
            ga, 
            simulation3dataset["cases"], 
            W_sim3, sim3fit1, 
            fitws(
                simulation3dataset["cases"], 
                simulation3dataset["Ns"], 
                sim3fit1
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=L"$w_{jt}$",
        )
        axs2 = plotrenewalequationsamples_cases!(
            ga, 
            simulation3dataset["cases"], 
            simulation3dataset["Ns"], 
            sim3fit1, 
            2;
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=L"Without \\ intervetion$$",
        )
        axs3 = plotrenewalequationsamples_cases!(
            ga, 
            simulation3dataset["cases"], 
            simulation3dataset["Ns"],
            sim3fit1,
            3;
            hidex=false,
            markersize=2, fittedparameter=:y_matrix_det_vec,
            xtitle="Time, days",
            ytitle=L"With \\ intervetion$$",
        )

    
        linkaxes!(axs2..., axs3...)
    
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
     
        colgap!(ga, 1, 5)  
        for r ∈ [ 1, 4 ] rowgap!(ga, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_w!(
            gb, 
            simulation3dataset["cases"], 
            W_sim3, 
            sim3fit2, 
            fitws(
                simulation3dataset["cases"], 
                simulation3dataset["Ns"], 
                sim3fit0leadlag
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=L"$w_{jt}$",
        )
        axs2 = plotrenewalequationsamples_cases!(
            gb, 
            simulation3dataset["cases"], 
            simulation3dataset["Ns"], 
            sim3fit2, 
            2;
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=L"Without \\ intervetion$$",
        )
        axs3 = plotrenewalequationsamples_cases!(
            gb, 
            simulation3dataset["cases"], 
            simulation3dataset["Ns"],
            sim3fit2,
            3;
            hidex=false,
            markersize=2, fittedparameter=:y_matrix_det_vec,
            xtitle="Time, days",
            ytitle=L"With \\ intervetion$$",
        )

        linkaxes!(axs2..., axs3...)
    
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                gb[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
    
        colgap!(gb, 1, 5)  
        for r ∈ [ 1, 4 ] rowgap!(gb, r, 5) end    
    end

    labelplots!([ "A", "B" ], [ ga, gb ])

    fig
end

safesave(plotsdir("subsetsim3plotsuppl.pdf"), subsetsim3plotsuppl)


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
            betafunctions=[ beta1a_a, beta1a_b ], 
            betafunctions_counterfactual=[ beta1a_a, beta1a_bcounterfactual ],
            infectiousduration=2.5,
            rhoclip = 45,
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
            rhoclip = 4,
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
            rhoclip = 4,
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
            rhoclip = 4,
            simcolour=COLOURVECTOR[3],
            columntitles=[ "Group 1", "Group 2" ],
            columntitlefontsize=10,
            markersize=2,
            xtitle="Time, days",
        )

        colsize!(gu, 1, Auto(2.5))
        colsize!(gu, 2, Auto(3.5))
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
datafit1plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 587, 411 ))
    ga = GridLayout(fig[1, 1])
    plotrenewalequationsamples!(
        ga, allcovidcases, W_allcoviddata, selectpops, datafit1;
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
    
    fig
end

safesave(plotsdir("datafit1plot.pdf"), datafit1plot)

subsetdatafit1plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 587, 411 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_r0!(
        ga, allcovidcases, datafit1, 1;
        locationinds=[ 1:5; [ 9 ] ],
        plotcounterfactuals=true, 
        xticks=xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"~ \\ ~ \\ $\mathcal{R}_0$",
    )
    axs2 = plotrenewalequationsamples_cases!(
        ga, allcovidcases, selectpops, datafit1, 2;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
        fittedcolour=( COLOURVECTOR[2], 0.75 ), 
        xticks=xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"~ \\ ~ \\ Without \\ intervetion$$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, allcovidcases, selectpops, datafit1, 3;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"~ \\ ~ \\ With \\ intervetion$$",
    )
    axs4 = plotrenewalequationsamples_causaleffect!(
        ga, allcovidcases, nothing, selectpops, datafit1, 4;
        cumulativedifference=true,
        fittedparameter=:y_matrix_det_vec,
        counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
        locationinds=[ 1:5; [ 9 ] ],
        xticklabelrotation=-π/4,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
        xtitle="Date, 2020–2021",
        ytitle=L"~ \\ ~ \\ Cumulative \\ difference$$",
    )

    linkaxes!(axs2..., axs3...)

    for (i, ℓ) ∈ enumerate([ 
        "Halton", 
        "Knowsley", 
        "Liverpool", 
        "Sefton", 
        "St Helens",  
        "Wirral" 
    ])
        Label(
        ga[0, i], ℓ; 
        fontsize=10, halign=:left, tellwidth=false
    )
    end

    braceaxis = Axis(ga[2:4, 0])
    bracket!(braceaxis, 0, 0, 0, 10; color=:black)
    formataxis!(braceaxis; hidex=true, hidexticks=true, hidey=true, hideyticks=true)
    hidespines!(braceaxis, :l, :r, :t, :b)
    Label(
        ga[2:4, -1], L"Recorded infections per $100\,000$"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    colgap!(ga, 1, -5)  
    colgap!(ga, 2, 5)  
    for r ∈ [ 1, 5 ] rowgap!(ga, r, 5) end

    fig
end

safesave(plotsdir("subsetdatafit1plot.pdf"), subsetdatafit1plot)


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
datafit2plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 587, 411 ))
    ga = GridLayout(fig[1, 1])
    plotrenewalequationsamples!(
        ga, pil1covidcases, W_pil1coviddata, selectpops, datafit2;
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
    
    fig 
end

safesave(plotsdir("datafit2plot.pdf"), datafit2plot)

datafit2kv = keyvalues(datachain2, datafit2)
println(datafit2kv)

subsetdatafit2plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 587, 411 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_r0!(
        ga, pil1covidcases, datafit2, 1;
        locationinds=[ 1:5; [ 9 ] ],
        plotcounterfactuals=true, 
        xticks=xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"~ \\ ~ \\ $\mathcal{R}_0$",
    )
    axs2 = plotrenewalequationsamples_cases!(
        ga, pil1covidcases, selectpops, datafit2, 2;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
        fittedcolour=( COLOURVECTOR[2], 0.75 ), 
        xticks=xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"~ \\ ~ \\ Without \\ intervetion$$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, pil1covidcases, selectpops, datafit2, 3;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"~ \\ ~ \\ With \\ intervetion$$",
    )
    axs4 = plotrenewalequationsamples_causaleffect!(
        ga, pil1covidcases, nothing, selectpops, datafit2, 4;
        cumulativedifference=true,
        fittedparameter=:y_matrix_det_vec,
        counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
        locationinds=[ 1:5; [ 9 ] ],
        xticklabelrotation=-π/4,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
        xtitle="Date, 2020–2021",
        ytitle=L"~ \\ ~ \\ Cumulative \\ difference$$",
    )

    linkaxes!(axs2..., axs3...)

    for (i, ℓ) ∈ enumerate([ 
        "Halton", 
        "Knowsley", 
        "Liverpool", 
        "Sefton", 
        "St Helens",  
        "Wirral" 
    ])
        Label(
        ga[0, i], ℓ; 
        fontsize=10, halign=:left, tellwidth=false
    )
    end

    braceaxis = Axis(ga[2:4, 0])
    bracket!(braceaxis, 0, 0, 0, 10; color=:black)
    formataxis!(braceaxis; hidex=true, hidexticks=true, hidey=true, hideyticks=true)
    hidespines!(braceaxis, :l, :r, :t, :b)
    Label(
        ga[2:4, -1], L"Pillar 1 diagnoses per $100\,000$"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    colgap!(ga, 1, -5)  
    colgap!(ga, 2, 5)  
    for r ∈ [ 1, 5 ] rowgap!(ga, r, 5) end

    fig
end

safesave(plotsdir("subsetdatafit2plot.pdf"), subsetdatafit2plot)


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
datafit3plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 587, 411 ))
    ga = GridLayout(fig[1, 1])
    plotrenewalequationsamples!(
        ga, allcovidcases, W_allcoviddata, selectpops, datafit3;
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

    fig
end

safesave(plotsdir("datafit3plot.pdf"), datafit3plot)

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
=#

## Analysis 5 
# Pillar 1 test results with lead and lag

datachain5 = loadanalysisdictsasdf("datamodel5", 8, maxrounds, 1050)
plotchains(datachain5)
datafit5 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain5, masstesting; 
    initialvalues=pil1covidcases[1:124, :], 
    Ns=selectpops,
    timeknots=[ collect(1.0:28:216); [216] ],
    secondaryinterventions=[ 
        InterventionsMatrix([ 172, 172, 146, 172, 172, 217, 217, 217, 172 ], 216), 
        InterventionsMatrix([ 200, 200, 174, 200, 200, 217, 217, 217, 200 ], 216), 
    ],
)
datafit5plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 587, 411 ))
    ga = GridLayout(fig[1, 1])
    plotrenewalequationsamples!(
        ga, pil1covidcases, W_pil1coviddata, selectpops, datafit5;
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

    fig 
end

safesave(plotsdir("datafit5plot.pdf"), datafit5plot)


datafit5kv = keyvalues(datachain5, datafit5)
println(datafit5kv)

chainplot1 = let 
    @unpack colnames, plotnames_ind = processplotchains(datachain1)
    fig = with_theme(theme_latexfonts()) do 
        fig = Figure(; size=( 500, 700 ))
        ga = GridLayout(fig[1, 1])
        gb = GridLayout(fig[1, 2])
    
        plotchains!(ga, datachain1; colnames=colnames, plotnames_ind=plotnames_ind[1:13])
        plotchains!(gb, datachain1; colnames=colnames, plotnames_ind=plotnames_ind[14:end])
    
        fig 
    end
    fig
end

safesave(plotsdir("chainplot1.pdf"), chainplot1)

chainplot2 = let 
    @unpack colnames, plotnames_ind = processplotchains(datachain2)
    fig = with_theme(theme_latexfonts()) do 
        fig = Figure(; size=( 500, 700 ))
        ga = GridLayout(fig[1, 1])
        gb = GridLayout(fig[1, 2])
    
        plotchains!(ga, datachain2; colnames=colnames, plotnames_ind=plotnames_ind[1:13])
        plotchains!(gb, datachain2; colnames=colnames, plotnames_ind=plotnames_ind[14:end])
    
        fig 
    end
    fig
end

safesave(plotsdir("chainplot2.pdf"), chainplot2)

chainplot3 = let 
    @unpack colnames, plotnames_ind = processplotchains(datachain3)
    fig = with_theme(theme_latexfonts()) do 
        fig = Figure(; size=( 500, 700 ))
        ga = GridLayout(fig[1, 1])
        gb = GridLayout(fig[1, 2])
    
        plotchains!(ga, datachain3; colnames=colnames, plotnames_ind=plotnames_ind[1:14])
        plotchains!(gb, datachain3; colnames=colnames, plotnames_ind=plotnames_ind[15:end])
    
        fig 
    end
    fig
end

safesave(plotsdir("chainplot3.pdf"), chainplot3)





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



