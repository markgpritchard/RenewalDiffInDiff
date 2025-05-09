
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
            fig[1, 0], L"serial interval, days, $g(x)$"; 
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

datagenerationintervalplot = let
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 250 ))
        ax = Axis(fig[1, 1]; yticks=[ 0, 0.05, 0.1, 0.15 ])
        scatter!(ax, 0:20, COVIDSERIALINTERVAL[1:21]; color=:black)
        formataxis!(ax; trimspines=true, hidespines=( :t, :r ), setpoint=0.15)
    
        Label(
            fig[1, 0], L"serial interval, $g(x)$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(fig[2, 1], L"time, days, $x$"; fontsize=11.84, tellwidth=false)
    
        colgap!(fig.layout, 1, 5)
        rowgap!(fig.layout, 1, 5)
        fig
    end
    fig
end

safesave(plotsdir("datagenerationintervalplot.pdf"), datagenerationintervalplot)


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
println("$(last(_modelquantiles(
    sim1fit0discrete, :y_matrix_det_vec, :y_matrix_det_vec_counterfactual, 2; 
    cumulativedifference=true
))) additional cases")
(_modelquantiles(
    sim1fit0discrete, :y_matrix_det_vec, :y_matrix_det_vec_counterfactual, 2;

        cumulativedifference=true  
)
) .* 100_000 / simulation1dataset["Ns"][2]

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
println("$(last(_modelquantiles(
    sim1fit1discrete, :y_matrix_det_vec, :y_matrix_det_vec_counterfactual, 2; 
    cumulativedifference=true
))) additional cases")
(sum(simulation1dataset["cases"][:, 2]) - sum(simulation1dataset["cases_counterfactual"][:, 2])) * 100_000 / simulation1dataset["Ns"][2]
(_modelquantiles(
           sim1fit1discrete, :y_matrix_det_vec, :y_matrix_det_vec_counterfactual, 2;

               cumulativedifference=true  
       )
       ) .* 100_000 / simulation1dataset["Ns"][2]

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
            counterfactualcases=simulation1dataset["cases_counterfactual"], 
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle="No\nintervention",
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
            counterfactualcases=simulation1dataset["cases_counterfactual"], 
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

    labelplots!([ "A", "B" ], [ ga, gb ]; cols=[ 0, 1 ],)

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
(_modelquantiles(
    sim1fit0, :y_matrix_det_vec, :y_matrix_det_vec_counterfactual, 2;

        cumulativedifference=true  
)
) .* 100_000 / simulation1dataset["Ns"][2]

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
(_modelquantiles(
    sim1fit1, :y_matrix_det_vec, :y_matrix_det_vec_counterfactual, 2;

        cumulativedifference=true  
)
) .* 100_000 / simulation1dataset["Ns"][2]

subsetsim1plotcsuppl = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    g0 = GridLayout(fig[1, 0])
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
            hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
        )
        axs2 = plotrenewalequationsamples_r0!(
            ga, simulation1dataset["cases_counterfactual"], sim1fit0, 2;
            betafunctions=[ beta1a, beta1bcounterfactual ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            rhoclip=3,
            yticks=[ 1.5, 2.5 ],
            ytitle=L"$\mathcal{R}_0$",
        )
        axs3 = plotrenewalequationsamples_cases!(
            ga, 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["Ns"], 
            sim1fit0, 
            3;
            counterfactualcases=simulation1dataset["cases_counterfactual"], 
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle="No\nintervention",
        )
        axs4 = plotrenewalequationsamples_cases!(
            ga, 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["Ns"],
            sim1fit0,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle="Intervention",
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            ga, simulation1dataset["cases_counterfactual"], 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["Ns"], 
            sim1fit0, 
            5;
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
            sim1fit0, 
            fitws(
                simulation1dataset["cases"], 
                simulation1dataset["Ns"], 
                sim1fit1
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=nothing,
            yticks=[ -1, 0, 1, 2 ], 
        )
        setvalue!(axs1[1], -1)
        axs2 = plotrenewalequationsamples_r0!(
            gb, simulation1dataset["cases"], sim1fit1, 2;
            betafunctions=[ beta1a, beta1b ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=nothing,
            rhoclip=3,
            yticks=[ 1, 2, 3 ],
        )
        setvalue!(axs2[1], 1)
        axs3 = plotrenewalequationsamples_cases!(
            gb, 
            simulation1dataset["cases"], 
            simulation1dataset["Ns"], 
            sim1fit1, 
            3;
            counterfactualcases=simulation1dataset["cases_counterfactual"], 
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=nothing,
        )
        axs4 = plotrenewalequationsamples_cases!(
            gb, 
            simulation1dataset["cases"], 
            simulation1dataset["Ns"],
            sim1fit1,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle=nothing,
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            gb, simulation1dataset["cases"], 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["Ns"], 
            sim1fit1, 
            5;
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

    labelplots!([ "A", "B" ], [ ga, gb ]; cols=[ 0, 1 ],)

    colgap!(fig.layout, 1, -5)
    colsize!(fig.layout, 0, Auto(0.05))
    colsize!(fig.layout, 2, Auto(0.75))

    fig
end

safesave(plotsdir("subsetsim1plotcsuppl.pdf"), subsetsim1plotcsuppl)

subsetsim1plotc = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 180 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    
    let
        axs1 = plotrenewalequationsamples_causaleffect!(
            ga, simulation1dataset["cases_counterfactual"], 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["Ns"], 
            sim1fit0, 
            1;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle="Cumulative\ndifference",
        )
        formataxis!(
            axs1[1]; 
            hidespines=( :r, :t ), trimspines=true,
        )
        formataxis!(
            axs1[2]; 
            hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t ), trimspines=true,
        )
    
        interventionax = Axis(ga[1, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )

        linkxaxes!(axs1[2], interventionax)

        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
   
        for r ∈ 1:2 rowgap!(ga, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_causaleffect!(
            gb, simulation1dataset["cases"], 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["Ns"], 
            sim1fit1, 
            1;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle=nothing,
        )

        formataxis!(
            axs1[1]; 
            hidespines=( :r, :t ), trimspines=true,
        )
        formataxis!(
            axs1[2]; 
            hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t ), trimspines=true,
        )

        interventionax = Axis(gb[1, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )

        linkxaxes!(axs1[2], interventionax)

        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                gb[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
    
        for r ∈ 1:2 rowgap!(gb, r, 5) end
    end

    colgap!(ga, 1, 5) 

    labelplots!([ "A", "B" ], [ ga, gb ]; cols=[ 0, 1 ],)

    colsize!(fig.layout, 2, Auto(0.75))

    fig
end

safesave(plotsdir("subsetsim1plotc.pdf"), subsetsim1plotc)


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

subsetsim2plotsuppl = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 800 ))
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
            hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
        )
        axs2 = plotrenewalequationsamples_r0!(
            ga, simulation2dataset["cases"], sim2fit0, 2;
            betafunctions=[ beta2a, beta2b, beta2c ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            rhoclip=3,
            yticks=[ 1, 2, 3 ],
            ytitle=L"$\mathcal{R}_0$",
        )
        setvalue!(axs2[1], 1)
        axs3 = plotrenewalequationsamples_cases!(
            ga, 
            simulation2dataset["cases"], 
            simulation2dataset["Ns"], 
            sim2fit0, 
            3;
            counterfactualcases=simulation2dataset["cases_counterfactual"], 
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle="No\nintervention",
        )
        axs4 = plotrenewalequationsamples_cases!(
            ga, 
            simulation2dataset["cases"], 
            simulation2dataset["Ns"],
            sim2fit0,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle="Intervention",
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            ga, simulation2dataset["cases"], 
            simulation2dataset["cases_counterfactual"], 
            simulation2dataset["Ns"], 
            sim2fit0, 
            5;
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
                for i ∈ 2:3
                    formataxis!(
                        axs[i]; 
                        hidey=true, hideyticks=true, 
                        hidespines=( :l, :r, :t ), trimspines=true,
                    )
                end
            else
                formataxis!(
                    axs[1]; 
                    hidex=true, hidexticks=true, 
                    hidespines=( :r, :t, :b ), trimspines=true,
                )
                for i ∈ 2:3
                    formataxis!(
                        axs[i]; 
                        hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                        hidespines=( :l, :r, :t, :b ), trimspines=true,
                    )
                end
            end
        end
    
        interventionax = Axis(ga[1:5, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
        interventionax2 = Axis(ga[1:5, 3])
        vlines!(interventionax2, 30; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax2; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )

        linkaxes!(axs3..., axs4...)
        linkxaxes!(axs1[2], axs2[2], axs3[2], axs4[2], axs5[2], interventionax, interventionax2)

        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2", "Group 3" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
        Label(
            ga[3:5, -1], L"Diagnoses, per $100\,000$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )

        colgap!(ga, 1, 5)  
        colgap!(ga, 2, 5)  

        for r ∈ [ 1, 6 ] rowgap!(ga, r, 5) end
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
            ytitle=L"$\ln\mathcal{R}_e$",
            #yticks=[ -1, 0, 1, 2 ], 
        )
        setvalue!(axs1[1], -1)
        axs2 = plotrenewalequationsamples_r0!(
            gb, simulation2dataset["cases"], sim2fit1, 2;
            betafunctions=[ beta2a, beta2b, beta2c ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=L"$\mathcal{R}_0$",
            rhoclip=3,
            yticks=[ 1, 2, 3 ],
        )
        setvalue!(axs2[1], 1)
        axs3 = plotrenewalequationsamples_cases!(
            gb, 
            simulation2dataset["cases"], 
            simulation2dataset["Ns"], 
            sim2fit1, 
            3;
            counterfactualcases=simulation2dataset["cases_counterfactual"], 
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle="No\nintervention",
        )
        axs4 = plotrenewalequationsamples_cases!(
            gb, 
            simulation2dataset["cases"], 
            simulation2dataset["Ns"],
            sim2fit1,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle="Intervention",
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            gb, simulation2dataset["cases"], 
            simulation2dataset["cases_counterfactual"], 
            simulation2dataset["Ns"], 
            sim2fit1, 
            5;
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
                for i ∈ 2:3
                    formataxis!(
                        axs[i]; 
                        hidey=true, hideyticks=true, 
                        hidespines=( :l, :r, :t ), trimspines=true,
                    )
                end
            else
                formataxis!(
                    axs[1]; 
                    hidex=true, hidexticks=true, 
                    hidespines=( :r, :t, :b ), trimspines=true,
                )
                for i ∈ 2:3
                    formataxis!(
                        axs[i]; 
                        hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                        hidespines=( :l, :r, :t, :b ), trimspines=true,
                    )
                end
            end
        end
    
        interventionax = Axis(gb[1:5, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )
        interventionax2 = Axis(gb[1:5, 3])
        vlines!(interventionax2, 30; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax2; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )

        linkaxes!(axs3..., axs4...)
        linkxaxes!(axs1[2], axs2[2], axs3[2], axs4[2], axs5[2], interventionax, interventionax2)

        Label(
            gb[3:5, -1], L"Diagnoses, per $100\,000$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )

        colgap!(gb, 1, 5)  
        colgap!(gb, 2, 5)  
        rowgap!(gb, 5, 5) 
    end

    labelplots!([ "A", "B" ], [ ga, gb ]; cols=-1, rows=[0, 1])

    fig
end

safesave(plotsdir("subsetsim2plotsuppl.pdf"), subsetsim2plotsuppl)


subsetsim2plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 275 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[2, 1])
    
    let
        axs1 = plotrenewalequationsamples_causaleffect!(
            ga, simulation2dataset["cases"], 
            simulation2dataset["cases_counterfactual"], 
            simulation2dataset["Ns"], 
            sim2fit0, 
            1;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle="Cumulative\ndifference",
        )
        formataxis!(
            axs1[1]; 
            hidespines=( :r, :t ), trimspines=true,
        )
        for i ∈ 2:3
            formataxis!(
                axs1[i]; 
                hidey=true, hideyticks=true, 
                hidespines=( :l, :r, :t ), trimspines=true,
            )
        end
    
        interventionax = Axis(ga[1, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
        interventionax2 = Axis(ga[1, 3])
        vlines!(interventionax2, 30; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax2; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )

        linkxaxes!(axs1[2], interventionax, interventionax2)

        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2", "Group 3" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end

        colgap!(ga, 1, 5)  
        colgap!(ga, 2, 5)  

        for r ∈ [ 1, 2 ] rowgap!(ga, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_causaleffect!(
            gb, simulation2dataset["cases"], 
            simulation2dataset["cases_counterfactual"], 
            simulation2dataset["Ns"], 
            sim2fit1, 
            1;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle="Cumulative\ndifference",
        )
        formataxis!(
            axs1[1]; 
            hidespines=( :r, :t ), trimspines=true,
        )
        for i ∈ 2:3
            formataxis!(
                axs1[i]; 
                hidey=true, hideyticks=true, 
                hidespines=( :l, :r, :t ), trimspines=true,
            )
        end
    
        interventionax = Axis(gb[1, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )
        interventionax2 = Axis(gb[1, 3])
        vlines!(interventionax2, 30; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax2; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )

        linkxaxes!(axs1[2], interventionax, interventionax2)

        colgap!(gb, 1, 5)  
        rowgap!(gb, 1, 5) 

        rowsize!(fig.layout, 2, Auto(0.95))
    end

    labelplots!([ "A", "B" ], [ ga, gb ]; cols=0, rows=[0, 1])

    fig
end

safesave(plotsdir("subsetsim2plot.pdf"), subsetsim2plot)


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
sim3fit0kv = keyvalues(sim3chain0, sim3fit0)
println(sim3fit0kv)

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
sim3fit0kvll = keyvalues(sim3chain0leadlag, sim3fit0leadlag)
println(sim3fit0kvll)
exp.(quantile(getproperty(sim3chain0leadlag, "logsecondarydelta1.logsecondarydelta"), [ 0.05, 0.5, 0.95 ]))
exp.(quantile(getproperty(sim3chain0leadlag, "logsecondarydelta2.logsecondarydelta"), [ 0.05, 0.5, 0.95 ]))

subsetsim3plot0suppl = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    g0 = GridLayout(fig[1, 0])
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
            hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
        )
        axs2 = plotrenewalequationsamples_r0!(
            ga, simulation3dataset["cases_counterfactual"], sim3fit0, 2;
            betafunctions=[ beta3a, beta3bcounterfactual ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            rhoclip=3,
            #yticks=[ 1.5, 2.5 ],
            ytitle=L"$\mathcal{R}_0$",
        )
        axs3 = plotrenewalequationsamples_cases!(
            ga, 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit0, 
            3;
            counterfactualcases=simulation3dataset["cases_counterfactual"], 
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle="No\nintervention",
        )
        axs4 = plotrenewalequationsamples_cases!(
            ga, 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"],
            sim3fit0,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle="Intervention",
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            ga, simulation3dataset["cases_counterfactual"], 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit0, 
            5;
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
            ytitle=nothing,
            #yticks=[ -1, 0, 1, 2 ], 
        )
        setvalue!(axs1[1], -1)
        axs2 = plotrenewalequationsamples_r0!(
            gb, simulation3dataset["cases_counterfactual"], sim3fit0leadlag, 2;
            betafunctions=[ beta3a, beta3bcounterfactual ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=nothing,
            rhoclip=3,
            #yticks=[ 1, 2, 3 ],
        )
        setvalue!(axs2[1], 1)
        axs3 = plotrenewalequationsamples_cases!(
            gb, 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit0leadlag, 
            3;
            counterfactualcases=simulation3dataset["cases_counterfactual"], 
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=nothing,
        )
        axs4 = plotrenewalequationsamples_cases!(
            gb, 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"],
            sim3fit0leadlag,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle=nothing,
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            gb, simulation3dataset["cases_counterfactual"], 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit0leadlag, 
            5;
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

    labelplots!([ "A", "B" ], [ ga, gb ]; cols=[ 0, 1 ],)

    colgap!(fig.layout, 1, -5)
    colsize!(fig.layout, 0, Auto(0.05))
    colsize!(fig.layout, 2, Auto(0.75))

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
sim3fit2kv = keyvalues(sim3chain2, sim3fit2)
println(sim3fit2kv)

subsetsim3plotsuppl = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    g0 = GridLayout(fig[1, 0])
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
            hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
        )
        axs2 = plotrenewalequationsamples_r0!(
            ga, simulation3dataset["cases"], sim3fit0, 2;
            betafunctions=[ beta3a, beta3b ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            rhoclip=3,
            yticks=[ 1, 2, 3 ],
            ytitle=L"$\mathcal{R}_0$",
        )
        setvalue!(axs2[1], 3)
        axs3 = plotrenewalequationsamples_cases!(
            ga, 
            simulation3dataset["cases"], 
            simulation3dataset["Ns"], 
            sim3fit1, 
            3;
            counterfactualcases=simulation3dataset["cases_counterfactual"], 
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle="No\nintervention",
        )
        axs4 = plotrenewalequationsamples_cases!(
            ga, 
            simulation3dataset["cases"], 
            simulation3dataset["Ns"],
            sim3fit1,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle="Intervention",
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            ga, simulation3dataset["cases"], 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit1, 
            5;
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
            simulation3dataset["cases"], 
            W_sim3, 
            sim3fit2, 
            fitws(
                simulation3dataset["cases"], 
                simulation3dataset["Ns"], 
                sim3fit2
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=nothing,
            #yticks=[ -1, 0, 1, 2 ], 
        )
        setvalue!(axs1[1], -1)
        axs2 = plotrenewalequationsamples_r0!(
            gb, simulation3dataset["cases"], sim3fit2, 2;
            betafunctions=[ beta3a, beta3b ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=nothing,
            rhoclip=3.6,
            yticks=[ 1, 2, 3 ],
        )
        setvalue!(axs2[1], 1)
        axs3 = plotrenewalequationsamples_cases!(
            gb, 
            simulation3dataset["cases"], 
            simulation3dataset["Ns"], 
            sim3fit2, 
            3;
            counterfactualcases=simulation3dataset["cases_counterfactual"], 
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=nothing,
        )
        axs4 = plotrenewalequationsamples_cases!(
            gb, 
            simulation3dataset["cases"], 
            simulation3dataset["Ns"],
            sim3fit2,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle=nothing,
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            gb, simulation3dataset["cases"], 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit2, 
            5;
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

    labelplots!([ "A", "B" ], [ ga, gb ]; cols=[ 0, 1 ],)

    colgap!(fig.layout, 1, -5)
    colsize!(fig.layout, 0, Auto(0.05))
    colsize!(fig.layout, 2, Auto(0.75))

    fig
end

safesave(plotsdir("subsetsim3plotsuppl.pdf"), subsetsim3plotsuppl)

subsetsim3plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 300 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    gc = GridLayout(fig[2, 1])
    gd = GridLayout(fig[2, 2])

    let
        axs1 = plotrenewalequationsamples_causaleffect!(
            ga, simulation3dataset["cases_counterfactual"], 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit0, 
            1;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle="Cumulative\ndifference",
        )

        formataxis!(
            axs1[1]; 
            hidespines=( :r, :t ), trimspines=true,
        )
        formataxis!(
            axs1[2]; 
            hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t ), trimspines=true,
        )
    
        interventionax = Axis(ga[1, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )

        linkxaxes!(axs1[2], interventionax)

        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
   
        for r ∈ 1:2 rowgap!(ga, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_causaleffect!(
            gb, simulation3dataset["cases_counterfactual"], 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit0leadlag, 
            1;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle=nothing,
        )
        formataxis!(
            axs1[1]; 
            hidespines=( :r, :t ), trimspines=true,
        )
        formataxis!(
            axs1[2]; 
            hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t ), trimspines=true,
        )
    
        interventionax = Axis(gb[1, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )

        linkxaxes!(axs1[2], interventionax)

        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                gb[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
    
        for r ∈ 1:2 rowgap!(gb, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_causaleffect!(
            gc, simulation3dataset["cases"], 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit1, 
            1;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle="Cumulative\ndifference",
        )
        formataxis!(
            axs1[1]; 
            hidespines=( :r, :t ), trimspines=true,
        )
        formataxis!(
            axs1[2]; 
            hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t ), trimspines=true,
        )
    
        interventionax = Axis(gc[1, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )

        linkxaxes!(axs1[2], interventionax)

        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                gc[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
   
        for r ∈ 1:2 rowgap!(gc, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_causaleffect!(
            gd, simulation3dataset["cases"], 
            simulation3dataset["cases_counterfactual"], 
            simulation3dataset["Ns"], 
            sim3fit2, 
            1;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle=nothing,
        )
        formataxis!(
            axs1[1]; 
            hidespines=( :r, :t ), trimspines=true,
        )
        formataxis!(
            axs1[2]; 
            hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t ), trimspines=true,
        )

        interventionax = Axis(gd[1, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )

        linkxaxes!(axs1[2], interventionax)
        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                gd[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
        for r ∈ 1:2 rowgap!(gd, r, 5) end
    end

    colgap!(ga, 1, 5) 
    colsize!(fig.layout, 2, Auto(0.75))

    labelplots!([ "A", "B", "C", "D" ], [ ga, gb, gc, gd ]; cols=[ 0, 1, 0, 1 ],)

    fig
end

safesave(plotsdir("subsetsim3plot.pdf"), subsetsim3plot)


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
sim4fit0kv = keyvalues(sim4chain0, sim4fit0)
println(sim4fit0kv)

sim4chain1 = loadanalysisdictsasdf("sim4model1", 8, maxrounds, 410)
plotchains(sim4chain1)
sim4fit1 = samplerenewalequation_2sets(
    fseir, sim4chain1, simulation4dataset["interventions"]; 
    initialvalues=simulation4dataset["cases"][1:20, :], 
    Ns=simulation4dataset["Ns"], 
    timeknots=[ [ 1 ]; collect(11:89/4:100) ],
)
sim4fit1kv = keyvalues(sim4chain1, sim4fit1)
println(sim4fit1kv)

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
sim4fit2kv = keyvalues(sim4chain2, sim4fit2)
println(sim4fit2kv)

subsetsim4plotsuppl0 = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    g0 = GridLayout(fig[1, 0])
    ga = GridLayout(fig[1, 1])
    
    let
        axs1 = plotrenewalequationsamples_w!(
            ga, 
            simulation4dataset["cases_counterfactual"], 
            W_sim4_0, sim4fit0, 
            fitws(
                simulation4dataset["cases_counterfactual"], 
                simulation4dataset["Ns"], 
                sim4fit0
            ), 
            1;
            markersize=2,
            hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
        )
        axs2 = plotrenewalequationsamples_r0!(
            ga, simulation4dataset["cases_counterfactual"], sim4fit0, 2;
            betafunctions=[ beta4a, beta4bcounterfactual ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            rhoclip=3,
            #yticks=[ 1.5, 2.5 ],
            ytitle=L"$\mathcal{R}_0$",
        )
        axs3 = plotrenewalequationsamples_cases!(
            ga, 
            simulation4dataset["cases_counterfactual"], 
            simulation4dataset["Ns"], 
            sim4fit0, 
            3;
            counterfactualcases=simulation4dataset["cases_counterfactual"], 
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle="No\nintervention",
        )
        axs4 = plotrenewalequationsamples_cases!(
            ga, 
            simulation4dataset["cases_counterfactual"], 
            simulation4dataset["Ns"],
            sim4fit0,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle="Intervention",
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            ga, simulation4dataset["cases_counterfactual"], 
            simulation4dataset["cases_counterfactual"], 
            simulation4dataset["Ns"], 
            sim4fit0, 
            5;
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

    Label(
        g0[3:5, 0], L"Diagnoses, per $100\,000$"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    colgap!(ga, 1, 5) 

    colgap!(fig.layout, 1, -5)
    colsize!(fig.layout, 0, Auto(0.05))

    fig
end

safesave(plotsdir("subsetsim4plotsuppl0.pdf"), subsetsim4plotsuppl0)

subsetsim4plotsuppl = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    g0 = GridLayout(fig[1, 0])
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    
    let
        axs1 = plotrenewalequationsamples_w!(
            ga, 
            simulation4dataset["cases"], 
            W_sim4, sim4fit1, 
            fitws(
                simulation4dataset["cases"], 
                simulation4dataset["Ns"], 
                sim4fit1
            ), 
            1;
            markersize=2,
            hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
        )
        axs2 = plotrenewalequationsamples_r0!(
            ga, simulation4dataset["cases"], sim4fit1, 2;
            betafunctions=[ beta4a, beta4b ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            rhoclip=3,
            #yticks=[ 1.5, 2.5 ],
            ytitle=L"$\mathcal{R}_0$",
        )
        axs3 = plotrenewalequationsamples_cases!(
            ga, 
            simulation4dataset["cases"], 
            simulation4dataset["Ns"], 
            sim4fit1, 
            3;
            counterfactualcases=simulation4dataset["cases_counterfactual"], 
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle="No\nintervention",
        )
        axs4 = plotrenewalequationsamples_cases!(
            ga, 
            simulation4dataset["cases"], 
            simulation4dataset["Ns"],
            sim4fit1,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle="Intervention",
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            ga, simulation4dataset["cases"], 
            simulation4dataset["cases_counterfactual"], 
            simulation4dataset["Ns"], 
            sim4fit1, 
            5;
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
            simulation4dataset["cases"], 
            W_sim4, 
            sim4fit2, 
            fitws(
                simulation4dataset["cases"], 
                simulation4dataset["Ns"], 
                sim4fit2
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=nothing,
            #yticks=[ -1, 0, 1, 2 ], 
        )
        setvalue!(axs1[1], -1)
        axs2 = plotrenewalequationsamples_r0!(
            gb, simulation4dataset["cases"], sim4fit2, 2;
            betafunctions=[ beta4a, beta4b ], infectiousduration=2.5,
            plotcounterfactuals=true, 
            ytitle=nothing,
            rhoclip=3,
            #yticks=[ 1, 2, 3 ],
        )
        setvalue!(axs2[1], 1)
        axs3 = plotrenewalequationsamples_cases!(
            gb, 
            simulation4dataset["cases"], 
            simulation4dataset["Ns"], 
            sim4fit2, 
            3;
            markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
            fittedcolour=( COLOURVECTOR[2], 0.75 ), 
            ytitle=nothing,
        )
        axs4 = plotrenewalequationsamples_cases!(
            gb, 
            simulation4dataset["cases"], 
            simulation4dataset["Ns"],
            sim4fit2,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle=nothing,
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            gb, simulation4dataset["cases"], 
            simulation4dataset["cases_counterfactual"], 
            simulation4dataset["Ns"], 
            sim4fit2, 
            5;
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

    labelplots!([ "A", "B" ], [ ga, gb ]; cols=[ 0, 1 ],)

    colgap!(fig.layout, 1, -5)
    colsize!(fig.layout, 0, Auto(0.05))
    colsize!(fig.layout, 2, Auto(0.75))

    fig
end

safesave(plotsdir("subsetsim4plotsuppl.pdf"), subsetsim4plotsuppl)

subsetsim4plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[2, 1])
    gc = GridLayout(fig[3, 1])

    let
        axs1 = plotrenewalequationsamples_causaleffect!(
            ga, simulation4dataset["cases_counterfactual"], 
            simulation4dataset["cases_counterfactual"], 
            simulation4dataset["Ns"], 
            sim4fit0, 
            1;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle="Cumulative\ndifference",
        )

        formataxis!(
            axs1[1]; 
            hidespines=( :r, :t ), trimspines=true,
        )
        formataxis!(
            axs1[2]; 
            hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t ), trimspines=true,
        )

        interventionax = Axis(ga[1, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )

        linkxaxes!(axs1[2], interventionax)

        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                ga[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
   
        for r ∈ 1:2 rowgap!(ga, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_causaleffect!(
            gb, simulation4dataset["cases"], 
            simulation4dataset["cases_counterfactual"], 
            simulation4dataset["Ns"], 
            sim4fit1, 
            1;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle="Cumulative\ndifference",
        )
        formataxis!(
            axs1[1]; 
            hidespines=( :r, :t ), trimspines=true,
        )
        formataxis!(
            axs1[2]; 
            hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t ), trimspines=true,
        )
    
        interventionax = Axis(gb[1, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )

        linkxaxes!(axs1[2], interventionax)

        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                gb[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
   
        for r ∈ 1:2 rowgap!(gb, r, 5) end
    end
    
    let
        axs1 = plotrenewalequationsamples_causaleffect!(
            gc, simulation4dataset["cases"], 
            simulation4dataset["cases_counterfactual"], 
            simulation4dataset["Ns"], 
            sim4fit2, 
            1;
            cumulativedifference=true,
            fittedparameter=:y_matrix_det_vec,
            counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
            ytickformat=(vs -> [ "$(round(Int, v))" for v ∈ vs ]),
            xtitle="Time, days",
            ytitle="Cumulative\ndifference",
        )

        formataxis!(
            axs1[1]; 
            hidespines=( :r, :t ), trimspines=true,
        )
        formataxis!(
            axs1[2]; 
            hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t ), trimspines=true,
        )
           
        interventionax = Axis(gc[1, 2])
        vlines!(interventionax, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        formataxis!(
            interventionax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )

        linkxaxes!(axs1[2], interventionax)

        for (i, ℓ) ∈ enumerate([ "Group 1", "Group 2" ])
            Label(
                gc[0, i], ℓ; 
                fontsize=10, halign=:left, tellwidth=false
            )
        end
    
        for r ∈ 1:2 rowgap!(gc, r, 5) end
    end

    colgap!(ga, 1, 5) 

    labelplots!([ "A", "B", "C" ], [ ga, gb, gc ]; cols=0,)

    fig
end

safesave(plotsdir("subsetsim4plot.pdf"), subsetsim4plot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Covid data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Priors 

### Mass testing

#### All data

Random.seed!(6)
datapriors1 = sample(datamodel1, Prior(), MCMCThreads(), 4000, 4)
datapriors1df = DataFrame(datapriors1)
plotchains(datapriors1df)
datapriorsfit1 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datapriors1df, masstesting; 
    initialvalues=allcovidcases[1:56, :], 
    Ns=selectpops,
    timeknots=[ collect(1.0:28:216); [216] ],
)

for i ∈ 16000:-1:1
    if isnan(maximum(datapriorsfit1[:y_matrix_det_vec][i])) || isnan(maximum(datapriorsfit1[:y_matrix_det_vec_counterfactual][i]))
        popat!(datapriorsfit1[:y_matrix_det_vec], i)
        popat!(datapriorsfit1[:y_matrix_det_vec_counterfactual], i)
    end
end

datapriorsfit1plot =  with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 500 ))
    ga = GridLayout(fig[1, 1:2])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        allcovidcases, 
        W_allcoviddata, datapriorsfit1, 
        fitws(
            allcovidcases, 
            selectpops, 
            datapriorsfit1
        ), 
        1;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_cases!(
        ga, allcovidcases, selectpops, datapriorsfit1, 2;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec,
        hidex=false,
        xticklabelrotation=-π/4,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        xtitle="Date, 2020–2021",
        ytitle=L"Incidence \\ per $100\,000$",
    )

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

    gb = GridLayout(fig[2, 1])
    axs3 = plotrenewalequationsamples_w!(
        gb, 
        allcovidcases, 
        W_allcoviddata, datapriorsfit1, 
        fitws(
            allcovidcases, 
            selectpops, 
            datapriorsfit1
        ), 
        1;
        locationinds=6:8,
        markersize=2,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        gb, allcovidcases, selectpops, datapriorsfit1, 2;
        locationinds=6:8,
        markersize=2, fittedparameter=:y_matrix_det_vec,
        hidex=false,
        xticklabelrotation=-π/4,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        xtitle="Date, 2020–2021",
        ytitle=L"Incidence \\ per $100\,000$",
    )

    for (i, ℓ) ∈ enumerate([ 
        "Warrington", 
        "W. Lancashire",  
        "Wigan" 
    ])
        Label(
        gb[0, i], ℓ; 
        fontsize=10, halign=:left, tellwidth=false
    )
    end

    for gl ∈ [ ga, gb ]
        colgap!(gl, 1, 5)  
        for r ∈ [ 1, 3 ] rowgap!(gl, r, 5) end

        if gl === ga 
            _imax = 6 
        else 
            _imax = 3 
        end
        for i ∈ 1:_imax 
            iax = Axis(gl[1:2, i]; xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ))
            formataxis!(
                iax; 
                hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                hidespines=( :l, :r, :t, :b)
            )
            iax.xgridstyle = ( :dot, :dense ) 
            iax.xgridwidth = 1
            iax.xgridvisible = true
            linkxaxes!(iax, axs1[i])
        end
    end
    
    for axs ∈ [ axs1, axs2, axs3, axs4 ]
        if axs === axs2 || axs === axs4
            formataxis!(
                axs[1]; 
                hidespines=( :r, :t ), trimspines=true,
            )
            if axs === axs2 
                _imax = 6 
            else 
                _imax = 3 
            end
            for i ∈ 2:_imax
                formataxis!(
                    axs[i]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            end
        else
            formataxis!(
                axs[1]; 
                hidex=true, hidexticks=true, 
                hidespines=( :r, :t, :b ), trimspines=true,
            )
            if axs === axs1 
                _imax = 6 
            else 
                _imax = 3 
            end
            for i ∈ 2:_imax 
                formataxis!(
                    axs[i]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    end

    colsize!(fig.layout, 1, Auto(1.3))

    fig
end

safesave(plotsdir("datapriorsfit1plot.pdf"), datapriorsfit1plot)

#### Pillar 1 data

Random.seed!(28)
datapriors2 = sample(datamodel2, Prior(), MCMCThreads(), 4000, 4)
datapriors2df = DataFrame(datapriors2)
plotchains(datapriors2df)
datapriorsfit2 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datapriors2df, masstesting; 
    initialvalues=pil1covidcases[1:56, :], 
    Ns=selectpops,
    timeknots=[ collect(1.0:28:216); [216] ],
)

for i ∈ 16000:-1:1
    if isnan(maximum(datapriorsfit2[:y_matrix_det_vec][i])) || isnan(maximum(datapriorsfit2[:y_matrix_det_vec_counterfactual][i]))
        popat!(datapriorsfit2[:y_matrix_det_vec], i)
        popat!(datapriorsfit2[:y_matrix_det_vec_counterfactual], i)
    end
end

datapriorsfit2plot =  with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 500 ))
    ga = GridLayout(fig[1, 1:2])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        pil1covidcases, 
        W_pil1coviddata, datapriorsfit2, 
        fitws(
            pil1covidcases, 
            selectpops, 
            datapriorsfit2
        ), 
        1;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_cases!(
        ga, pil1covidcases, selectpops, datapriorsfit2, 2;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec,
        hidex=false,
        xticklabelrotation=-π/4,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        xtitle="Date, 2020–2021",
        ytitle=L"Incidence \\ per $100\,000$",
    )

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

    gb = GridLayout(fig[2, 1])
    axs3 = plotrenewalequationsamples_w!(
        gb, 
        pil1covidcases, 
        W_pil1coviddata, datapriorsfit2, 
        fitws(
            pil1covidcases, 
            selectpops, 
            datapriorsfit2
        ), 
        1;
        locationinds=6:8,
        markersize=2,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        gb, pil1covidcases, selectpops, datapriorsfit2, 2;
        locationinds=6:8,
        markersize=2, fittedparameter=:y_matrix_det_vec,
        hidex=false,
        xticklabelrotation=-π/4,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        xtitle="Date, 2020–2021",
        ytitle=L"Incidence \\ per $100\,000$",
    )

    for (i, ℓ) ∈ enumerate([ 
        "Warrington", 
        "W. Lancashire",  
        "Wigan" 
    ])
        Label(
        gb[0, i], ℓ; 
        fontsize=10, halign=:left, tellwidth=false
    )
    end

    for gl ∈ [ ga, gb ]
        colgap!(gl, 1, 5)  
        for r ∈ [ 1, 3 ] rowgap!(gl, r, 5) end

        if gl === ga 
            _imax = 6 
        else 
            _imax = 3 
        end
        for i ∈ 1:_imax 
            iax = Axis(gl[1:2, i]; xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ))
            formataxis!(
                iax; 
                hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                hidespines=( :l, :r, :t, :b)
            )
            iax.xgridstyle = ( :dot, :dense ) 
            iax.xgridwidth = 1
            iax.xgridvisible = true
            linkxaxes!(iax, axs1[i])
        end
    end
    
    for axs ∈ [ axs1, axs2, axs3, axs4 ]
        if axs === axs2 || axs === axs4
            formataxis!(
                axs[1]; 
                hidespines=( :r, :t ), trimspines=true,
            )
            if axs === axs2 
                _imax = 6 
            else 
                _imax = 3 
            end
            for i ∈ 2:_imax
                formataxis!(
                    axs[i]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            end
        else
            formataxis!(
                axs[1]; 
                hidex=true, hidexticks=true, 
                hidespines=( :r, :t, :b ), trimspines=true,
            )
            if axs === axs1 
                _imax = 6 
            else 
                _imax = 3 
            end
            for i ∈ 2:_imax 
                formataxis!(
                    axs[i]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    end

    colsize!(fig.layout, 1, Auto(1.3))

    fig
end

safesave(plotsdir("datapriorsfit2plot.pdf"), datapriorsfit2plot)


### Mask use / requirements

Random.seed!(496)
maskingdatapriors1 = sample(maskingdatamodel2, Prior(), MCMCThreads(), 4000, 4)
maskingdatapriors1df = DataFrame(maskingdatapriors1)
plotchains(maskingdatapriors1df)
maskingdatapriorsfit1 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, maskingdatapriors1df, facialcoveringsrequired; 
    initialvalues=maskcovidcases[1:91, :], 
    Ns=POPULATION2020,
    timeknots=[ 1.0; collect(56.0:28:224); 257 ],
)

for i ∈ 16000:-1:1
    if isnan(maximum(maskingdatapriorsfit1[:y_matrix_det_vec][i])) || isnan(maximum(maskingdatapriorsfit1[:y_matrix_det_vec_counterfactual][i]))
        popat!(maskingdatapriorsfit1[:y_matrix_det_vec], i)
        popat!(maskingdatapriorsfit1[:y_matrix_det_vec_counterfactual], i)
    end
end

maskingdatapriorsfit1plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 250 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        maskcovidcases, 
        W_maskcoviddata, maskingdatapriorsfit1, 
        fitws(
            maskcovidcases, 
            POPULATION2020, 
            maskingdatapriorsfit1
        ), 
        1;
        markersize=2,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases, POPULATION2020, maskingdatapriorsfit1, 2;
        markersize=2, fittedparameter=:y_matrix_det_vec,
        hidex=false,
        xticklabelrotation=-π/4,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        xtitle="Date, 2020",
        ytitle=L"Incidence \\ per $100\,000$",
    )

    for (i, ℓ) ∈ enumerate([ 
        "England", 
        "Northern Ireland", 
        "Scotland",  
        "Wales" 
    ])
        Label(
        ga[0, i], ℓ; 
        fontsize=10, halign=:left, tellwidth=false
    )
    end

    colgap!(ga, 1, 5)  
    for r ∈ [ 1, 3 ] rowgap!(ga, r, 5) end
    for axs ∈ [ axs1, axs2 ]
        if axs === axs2 
            formataxis!(
                axs[1]; 
                hidespines=( :r, :t ), trimspines=true,
            )
            for i ∈ 2:4
                formataxis!(
                    axs[i]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            end
        else
            formataxis!(
                axs[1]; 
                hidex=true, hidexticks=true, 
                hidespines=( :r, :t, :b ), trimspines=true,
            )
            for i ∈ 2:4 
                formataxis!(
                    axs[i]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    end

    for i ∈ 1:4 
        iax = Axis(ga[1:2, i]; xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ))
        formataxis!(
            iax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )
        iax.xgridstyle = ( :dot, :dense ) 
        iax.xgridwidth = 1
        iax.xgridvisible = true
        linkxaxes!(iax, axs1[i])
    end

    fig
end

safesave(plotsdir("maskingdatapriorsfit1plot.pdf"), maskingdatapriorsfit1plot)


## Fitted data

datachain1 = loadanalysisdictsasdf("datamodel1", 8, maxrounds, 1010)

datachain1plot = with_theme(theme_latexfonts()) do  
    _names1 = [
        "log density",
        L"$\mu_{\zeta}$",
        L"$\sigma_{\zeta}^2$",
        L"$\ln\zeta_{1}$",
        L"$\ln\zeta_{2}$",
        L"$\ln\zeta_{3}$",
        L"$\ln\zeta_{4}$",
        L"$\ln\zeta_{5}$",
        L"$\ln\zeta_{6}$",
        L"$\ln\zeta_{7}$",
        L"$\ln\zeta_{8}$",
        L"$\ln\zeta_{9}$",
        L"$\mu_{\eta}$",
    ]
    _names2 = [
        L"$\sigma_{\eta}^2$",
        L"$\ln\eta(1)$",
        L"$\ln\eta(3)$",
        L"$\ln\eta(4)$",
        L"$\ln\eta(5)$",
        L"$\ln\eta(6)$",
        L"$\ln\eta(7)$",
        L"$\ln\eta(8)$",
        L"$\ln\eta(9)$",
        L"$\ln\tau_{\mathrm{ATT}}$",
        L"$\sigma^2$",
        L"$\theta$",
    ]

    fig = Figure(; size=( 500, 800 ))
    axs1 = [ 
        Axis(fig[i, 2*j-1], xticks=WilkinsonTicks(3), yticks=WilkinsonTicks(3)) 
        for i ∈ 1:13, j ∈ 1:2 
    ]
    axs2 = [ 
        Axis(fig[i, 2*j+3], xticks=WilkinsonTicks(3), yticks=WilkinsonTicks(3)) 
        for i ∈ 1:12, j ∈ 1:2 
    ]

    for i ∈ 4:-1:1
        _tdf = filter(:chain => x -> x == i, datachain1)
        lines!(
            axs1[1, 1], _tdf.iteration, _tdf.log_density; 
            color=COLOURVECTOR[i], linewidth=1,
        )
        density!(
            axs1[1, 2], _tdf.log_density; 
            color=( :white, 0 ), strokecolor=COLOURVECTOR[i], strokewidth=1,
        )
        for (j, v) ∈ enumerate(names(_tdf)[3:14])
            lines!(
                axs1[1+j, 1], _tdf.iteration, getproperty(_tdf, v); 
                color=COLOURVECTOR[i], linewidth=1,
            )
            density!(
                axs1[1+j, 2], getproperty(_tdf, v); 
                color=( :white, 0 ), strokecolor=COLOURVECTOR[i], strokewidth=1,
            )
        end
        for (j, v) ∈ enumerate(names(_tdf)[15:26])
            lines!(
                axs2[j, 1], _tdf.iteration, getproperty(_tdf, v); 
                color=COLOURVECTOR[i], linewidth=1,
            )
            density!(
                axs2[j, 2], getproperty(_tdf, v); 
                color=( :white, 0 ), strokecolor=COLOURVECTOR[i], strokewidth=1,
            )
        end
    end

    for i ∈ 1:13, j ∈ 1:2
        formataxis!(axs1[i, j]; hidespines=( :r, :t ), trimspines=true,)
        #Label(fig[i, 1], "Iteration"; fontsize=11.84, tellwidth=false)
        Label(fig[i, 2], "Density"; fontsize=11.84, rotation=π/2, tellheight=false,)
        Label(fig[i, 0], _names1[i]; fontsize=11.84, rotation=π/2, tellheight=false,)
        i == 13 && continue
        formataxis!(axs2[i, j]; hidespines=( :r, :t ), trimspines=true,)
        #Label(fig[i, 5], "Iteration"; fontsize=11.84, tellwidth=false)
        Label(fig[i, 6], "Density"; fontsize=11.84, rotation=π/2, tellheight=false,)
        Label(fig[i, 4], _names2[i]; fontsize=11.84, rotation=π/2, tellheight=false,)
    end

    Label(fig[14, 1], "Iteration"; fontsize=11.84, tellwidth=false,)
    Label(
        fig[13, 5], "Iteration"; 
        fontsize=11.84, tellwidth=false, tellheight=false, valign=:top,
    )
    for c ∈ [ 1, 3, 5, 7 ] colgap!(fig.layout, c, 5) end
    rowgap!(fig.layout, 13, 5)

    fig
end

safesave(plotsdir("datachain1plot.pdf"), datachain1plot)

datafit1 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain1, masstesting; 
    initialvalues=allcovidcases[1:56, :], 
    Ns=selectpops,
    timeknots=[ collect(1.0:28:216); [ 216 ] ],
)
datafit1kv = keyvalues(datachain1, datafit1)
println(datafit1kv)

subsetdatafit1plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        allcovidcases, 
        W_allcoviddata, datafit1, 
        fitws(
            allcovidcases, 
            selectpops, 
            datafit1
        ), 
        1;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_r0!(
        ga, allcovidcases, datafit1, 2;
        locationinds=[ 1:5; [ 9 ] ],
        plotcounterfactuals=true, 
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"$\mathcal{R}_0$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, allcovidcases, selectpops, datafit1, 3;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
        fittedcolour=( COLOURVECTOR[2], 0.75 ), 
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"Without \\ intervention$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, allcovidcases, selectpops, datafit1, 4;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"With \\ intervention$$",
    )
    axs5 = plotrenewalequationsamples_causaleffect!(
        ga, allcovidcases, nothing, selectpops, datafit1, 5;
        cumulativedifference=true,
        fittedparameter=:y_matrix_det_vec,
        counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
        locationinds=[ 1:5; [ 9 ] ],
        xticklabelrotation=-π/4,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
        xtitle="Date, 2020–2021",
        ytitle=L"Cumulative \\ difference$$",
    )

    linkaxes!(axs3..., axs4...)

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

    colgap!(ga, 1, 5)  
    for r ∈ [ 1, 6 ] rowgap!(ga, r, 5) end
    for axs ∈ [ axs1, axs2, axs3, axs4, axs5 ]
        if axs === axs5 
            formataxis!(
                axs[1]; 
                hidespines=( :r, :t ), trimspines=true,
            )
            for i ∈ 2:6
                formataxis!(
                    axs[i]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            end
        else
            formataxis!(
                axs[1]; 
                hidex=true, hidexticks=true, 
                hidespines=( :r, :t, :b ), trimspines=true,
            )
            for i ∈ 2:6 
                formataxis!(
                    axs[i]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    end

    for i ∈ 1:6 
        iax = Axis(ga[1:5, i]; xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ))
        if i == 3
            vlines!(iax, 160; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        else
            vlines!(iax, 186; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        end
        formataxis!(
            iax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )
        iax.xgridstyle=( :dot, :dense ) 
        iax.xgridwidth = 1
        iax.xgridvisible = true
        linkxaxes!(iax, axs1[i])
    end

    fig
end

safesave(plotsdir("subsetdatafit1plot.pdf"), subsetdatafit1plot)

subsetdatafit1plotb = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        allcovidcases, 
        W_allcoviddata, datafit1, 
        fitws(
            allcovidcases, 
            selectpops, 
            datafit1
        ), 
        1;
        locationinds=6:8,
        markersize=2,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_r0!(
        ga, allcovidcases, datafit1, 2;
        locationinds=6:8,
        plotcounterfactuals=true, 
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"$\mathcal{R}_0$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, allcovidcases, selectpops, datafit1, 3;
        locationinds=6:8,
        markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
        fittedcolour=( COLOURVECTOR[2], 0.75 ), 
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"Without \\ intervention$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, allcovidcases, selectpops, datafit1, 4;
        locationinds=6:8,
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"With \\ intervention$$",
    )
    axs5 = plotrenewalequationsamples_causaleffect!(
        ga, allcovidcases, nothing, selectpops, datafit1, 5;
        cumulativedifference=true,
        fittedparameter=:y_matrix_det_vec,
        counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
        locationinds=6:8,
        xticklabelrotation=-π/4,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
        xtitle="Date, 2020–2021",
        ytitle=L"Cumulative \\ difference$$",
    )

    linkaxes!(axs3..., axs4...)

    for (i, ℓ) ∈ enumerate([ 
        "Warrington", 
        "West Lancashire",  
        "Wigan" 
    ])
        Label(
        ga[0, i], ℓ; 
        fontsize=10, halign=:left, tellwidth=false
    )
    end

    colgap!(ga, 1, 5)  
    for r ∈ [ 1, 6 ] rowgap!(ga, r, 5) end
    for axs ∈ [ axs1, axs2, axs3, axs4, axs5 ]
        if axs === axs5 
            formataxis!(
                axs[1]; 
                hidespines=( :r, :t ), trimspines=true,
            )
            for i ∈ 2:3
                formataxis!(
                    axs[i]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            end
        else
            formataxis!(
                axs[1]; 
                hidex=true, hidexticks=true, 
                hidespines=( :r, :t, :b ), trimspines=true,
            )
            for i ∈ 2:3
                formataxis!(
                    axs[i]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    end

    for i ∈ 1:3
        iax = Axis(ga[1:5, i]; xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ))
        formataxis!(
            iax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )
        iax.xgridstyle=( :dot, :dense ) 
        iax.xgridwidth = 1
        iax.xgridvisible = true
        linkxaxes!(iax, axs1[i])
    end

    fig
end

safesave(plotsdir("subsetdatafit1plotb.pdf"), subsetdatafit1plotb)

datafit1kv = keyvalues(datachain1, datafit1)
println(datafit1kv)


## Analysis 2 
# Pillar 1 test results 

datachain2 = loadanalysisdictsasdf("datamodel2", 8, maxrounds, 1020)

pillar1chainsplot = with_theme(theme_latexfonts()) do  
    _names1 = [
        "log density",
        L"$\mu_{\zeta}$",
        L"$\sigma_{\zeta}^2$",
        L"$\ln\zeta_{1}$",
        L"$\ln\zeta_{2}$",
        L"$\ln\zeta_{3}$",
        L"$\ln\zeta_{4}$",
        L"$\ln\zeta_{5}$",
        L"$\ln\zeta_{6}$",
        L"$\ln\zeta_{7}$",
        L"$\ln\zeta_{8}$",
        L"$\ln\zeta_{9}$",
        L"$\mu_{\eta}$",
    ]
    _names2 = [
        L"$\sigma_{\eta}^2$",
        L"$\ln\eta(1)$",
        L"$\ln\eta(3)$",
        L"$\ln\eta(4)$",
        L"$\ln\eta(5)$",
        L"$\ln\eta(6)$",
        L"$\ln\eta(7)$",
        L"$\ln\eta(8)$",
        L"$\ln\eta(9)$",
        L"$\ln\tau_{\mathrm{ATT}}$",
        L"$\sigma^2$",
        L"$\theta$",
    ]

    fig = Figure(; size=( 500, 800 ))
    axs1 = [ 
        Axis(fig[i, 2*j-1], xticks=WilkinsonTicks(3), yticks=WilkinsonTicks(3)) 
        for i ∈ 1:13, j ∈ 1:2 
    ]
    axs2 = [ 
        Axis(fig[i, 2*j+3], xticks=WilkinsonTicks(3), yticks=WilkinsonTicks(3)) 
        for i ∈ 1:12, j ∈ 1:2 
    ]

    for i ∈ 4:-1:1
        _tdf = filter(:chain => x -> x == i, datachain2)
        lines!(
            axs1[1, 1], _tdf.iteration, _tdf.log_density; 
            color=COLOURVECTOR[i], linewidth=1,
        )
        density!(
            axs1[1, 2], _tdf.log_density; 
            color=( :white, 0 ), strokecolor=COLOURVECTOR[i], strokewidth=1,
        )
        for (j, v) ∈ enumerate(names(_tdf)[3:14])
            lines!(
                axs1[1+j, 1], _tdf.iteration, getproperty(_tdf, v); 
                color=COLOURVECTOR[i], linewidth=1,
            )
            density!(
                axs1[1+j, 2], getproperty(_tdf, v); 
                color=( :white, 0 ), strokecolor=COLOURVECTOR[i], strokewidth=1,
            )
        end
        for (j, v) ∈ enumerate(names(_tdf)[15:26])
            lines!(
                axs2[j, 1], _tdf.iteration, getproperty(_tdf, v); 
                color=COLOURVECTOR[i], linewidth=1,
            )
            density!(
                axs2[j, 2], getproperty(_tdf, v); 
                color=( :white, 0 ), strokecolor=COLOURVECTOR[i], strokewidth=1,
            )
        end
    end

    for i ∈ 1:13, j ∈ 1:2
        formataxis!(axs1[i, j]; hidespines=( :r, :t ), trimspines=true,)
        #Label(fig[i, 1], "Iteration"; fontsize=11.84, tellwidth=false)
        Label(fig[i, 2], "Density"; fontsize=11.84, rotation=π/2, tellheight=false,)
        Label(fig[i, 0], _names1[i]; fontsize=11.84, rotation=π/2, tellheight=false,)
        i == 13 && continue
        formataxis!(axs2[i, j]; hidespines=( :r, :t ), trimspines=true,)
        #Label(fig[i, 5], "Iteration"; fontsize=11.84, tellwidth=false)
        Label(fig[i, 6], "Density"; fontsize=11.84, rotation=π/2, tellheight=false,)
        Label(fig[i, 4], _names2[i]; fontsize=11.84, rotation=π/2, tellheight=false,)
    end

    Label(fig[14, 1], "Iteration"; fontsize=11.84, tellwidth=false,)
    Label(
        fig[13, 5], "Iteration"; 
        fontsize=11.84, tellwidth=false, tellheight=false, valign=:top,
    )
    for c ∈ [ 1, 3, 5, 7 ] colgap!(fig.layout, c, 5) end
    rowgap!(fig.layout, 13, 5)


    fig
end

safesave(plotsdir("pillar1chainsplot.pdf"), pillar1chainsplot)

datafit2 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain2, masstesting; 
    initialvalues=pil1covidcases[1:56, :], 
    Ns=selectpops,
    timeknots=[ collect(1.0:28:216); [216] ],
)

datafit2plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        pil1covidcases, 
        W_pil1coviddata, datafit2, 
        fitws(
            pil1covidcases, 
            selectpops, 
            datafit2
        ), 
        1;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_r0!(
        ga, pil1covidcases, datafit2, 2;
        locationinds=[ 1:5; [ 9 ] ],
        plotcounterfactuals=true, 
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"$\mathcal{R}_0$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, pil1covidcases, selectpops, datafit2, 3;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
        fittedcolour=( COLOURVECTOR[2], 0.75 ), 
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"Without \\ intervention$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, pil1covidcases, selectpops, datafit2, 4;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"With \\ intervention$$",
    )
    axs5 = plotrenewalequationsamples_causaleffect!(
        ga, pil1covidcases, nothing, selectpops, datafit2, 5;
        cumulativedifference=true,
        fittedparameter=:y_matrix_det_vec,
        counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
        locationinds=[ 1:5; [ 9 ] ],
        xticklabelrotation=-π/4,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
        xtitle="Date, 2020–2021",
        ytitle=L"Cumulative \\ difference$$",
    )

    linkaxes!(axs3..., axs4...)

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

    colgap!(ga, 1, 5)  
    for r ∈ [ 1, 6 ] rowgap!(ga, r, 5) end
    for axs ∈ [ axs1, axs2, axs3, axs4, axs5 ]
        if axs === axs5 
            formataxis!(
                axs[1]; 
                hidespines=( :r, :t ), trimspines=true,
            )
            for i ∈ 2:6
                formataxis!(
                    axs[i]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            end
        else
            formataxis!(
                axs[1]; 
                hidex=true, hidexticks=true, 
                hidespines=( :r, :t, :b ), trimspines=true,
            )
            for i ∈ 2:6 
                formataxis!(
                    axs[i]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    end

    for i ∈ 1:6 
        iax = Axis(ga[1:5, i]; xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ))
        if i == 3
            vlines!(iax, 160; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        else
            vlines!(iax, 186; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        end
        formataxis!(
            iax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )
        iax.xgridstyle=( :dot, :dense ) 
        iax.xgridwidth = 1
        iax.xgridvisible = true
        linkxaxes!(iax, axs1[i])
    end

    fig
end

safesave(plotsdir("datafit2plot.pdf"), datafit2plot)

datafit2plotb = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        pil1covidcases, 
        W_pil1coviddata, datafit2, 
        fitws(
            pil1covidcases, 
            selectpops, 
            datafit2
        ), 
        1;
        locationinds=6:8,
        markersize=2,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_r0!(
        ga, pil1covidcases, datafit2, 2;
        locationinds=6:8,
        plotcounterfactuals=true, 
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"$\mathcal{R}_0$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, pil1covidcases, selectpops, datafit2, 3;
        locationinds=6:8,
        markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
        fittedcolour=( COLOURVECTOR[2], 0.75 ), 
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"Without \\ intervention$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, pil1covidcases, selectpops, datafit2, 4;
        locationinds=6:8,
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"With \\ intervention$$",
    )
    axs5 = plotrenewalequationsamples_causaleffect!(
        ga, pil1covidcases, nothing, selectpops, datafit2, 5;
        cumulativedifference=true,
        fittedparameter=:y_matrix_det_vec,
        counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
        locationinds=6:8,
        xticklabelrotation=-π/4,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
        xtitle="Date, 2020–2021",
        ytitle=L"Cumulative \\ difference$$",
    )

    linkaxes!(axs3..., axs4...)

    for (i, ℓ) ∈ enumerate([ 
        "Warrington", 
        "West Lancashire",  
        "Wigan" 
    ])
        Label(
        ga[0, i], ℓ; 
        fontsize=10, halign=:left, tellwidth=false
    )
    end

    colgap!(ga, 1, 5)  
    for r ∈ [ 1, 6 ] rowgap!(ga, r, 5) end
    for axs ∈ [ axs1, axs2, axs3, axs4, axs5 ]
        if axs === axs5 
            formataxis!(
                axs[1]; 
                hidespines=( :r, :t ), trimspines=true,
            )
            for i ∈ 2:3
                formataxis!(
                    axs[i]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            end
        else
            formataxis!(
                axs[1]; 
                hidex=true, hidexticks=true, 
                hidespines=( :r, :t, :b ), trimspines=true,
            )
            for i ∈ 2:3
                formataxis!(
                    axs[i]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    end

    for i ∈ 1:3
        iax = Axis(ga[1:5, i]; xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ))
        formataxis!(
            iax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )
        iax.xgridstyle=( :dot, :dense ) 
        iax.xgridwidth = 1
        iax.xgridvisible = true
        linkxaxes!(iax, axs1[i])
    end

    fig
end

safesave(plotsdir("datafit2plotb.pdf"), datafit2plotb)

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
        ytitle=L"~ \\ ~ \\ Without \\ intervention$$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, pil1covidcases, selectpops, datafit2, 3;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"~ \\ ~ \\ With \\ intervention$$",
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
    initialvalues=allcovidcases[1:56, :], 
    Ns=selectpops,
    timeknots=[ collect(1.0:28:216); [216] ],
    secondaryinterventions=[ 
        InterventionsMatrix([ 172, 172, 146, 172, 172, 217, 217, 217, 172 ], 216), 
        InterventionsMatrix([ 200, 200, 174, 200, 200, 217, 217, 217, 200 ], 216), 
    ],
)

datafit3kv = keyvalues(datachain3, datafit3)
println(datafit3kv)


datafit3plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        allcovidcases, 
        W_allcoviddata, datafit3, 
        fitws(
            allcovidcases, 
            selectpops, 
            datafit3
        ), 
        1;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_r0!(
        ga, allcovidcases, datafit3, 2;
        locationinds=[ 1:5; [ 9 ] ],
        plotcounterfactuals=true, 
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"$\mathcal{R}_0$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, allcovidcases, selectpops, datafit3, 3;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
        fittedcolour=( COLOURVECTOR[2], 0.75 ), 
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"Without \\ intervention$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, allcovidcases, selectpops, datafit3, 4;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"With \\ intervention$$",
    )
    axs5 = plotrenewalequationsamples_causaleffect!(
        ga, allcovidcases, nothing, selectpops, datafit3, 5;
        cumulativedifference=true,
        fittedparameter=:y_matrix_det_vec,
        counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
        locationinds=[ 1:5; [ 9 ] ],
        xticklabelrotation=-π/4,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ),
        xtitle="Date, 2020–2021",
        ytitle=L"Cumulative \\ difference$$",
    )

    linkaxes!(axs3..., axs4...)

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

    colgap!(ga, 1, 5)  
    for r ∈ [ 1, 6 ] rowgap!(ga, r, 5) end
    for axs ∈ [ axs1, axs2, axs3, axs4, axs5 ]
        if axs === axs5 
            formataxis!(
                axs[1]; 
                hidespines=( :r, :t ), trimspines=true,
            )
            for i ∈ 2:6
                formataxis!(
                    axs[i]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            end
        else
            formataxis!(
                axs[1]; 
                hidex=true, hidexticks=true, 
                hidespines=( :r, :t, :b ), trimspines=true,
            )
            for i ∈ 2:6 
                formataxis!(
                    axs[i]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    end

    for i ∈ 1:6 
        iax = Axis(ga[1:5, i]; xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ))
        if i == 3
            vlines!(iax, 160; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        else
            vlines!(iax, 186; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        end
        formataxis!(
            iax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )
        iax.xgridstyle=( :dot, :dense ) 
        iax.xgridwidth = 1
        iax.xgridvisible = true
        linkxaxes!(iax, axs1[i])
    end

    fig
end

safesave(plotsdir("datafit3plot.pdf"), datafit3plot)


## Analysis 5 
# Pillar 1 test results with lead and lag

datachain5 = loadanalysisdictsasdf("datamodel5", 8, maxrounds, 1050)
plotchains(datachain5)
datafit5 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datachain5, masstesting; 
    initialvalues=pil1covidcases[1:56, :], 
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

maskingdatachain1chainsplot = with_theme(theme_latexfonts()) do  
    _names1 = [
        "log density",
        L"$\mu_{\zeta}$",
        L"$\sigma_{\zeta}^2$",
        L"$\ln\zeta_{1}$",
        L"$\ln\zeta_{2}$",
        L"$\ln\zeta_{3}$",
        L"$\ln\zeta_{4}$",
        L"$\mu_{\eta}$",
        L"$\sigma_{\eta}^2$",
    ]
    _names2 = [
        L"$\ln\eta(1)$",
        L"$\ln\eta(3)$",
        L"$\ln\eta(4)$",
        L"$\ln\eta(5)$",
        L"$\ln\eta(6)$",
        L"$\ln\eta(7)$",
        L"$\ln\tau_{\mathrm{ATT}}$",
        L"$\sigma^2$",
        L"$\theta$",
    ]

    fig = Figure(; size=( 500, 800 ))
    axs1 = [ 
        Axis(fig[i, 2*j-1], xticks=WilkinsonTicks(3), yticks=WilkinsonTicks(3)) 
        for i ∈ 1:9, j ∈ 1:2 
    ]
    axs2 = [ 
        Axis(fig[i, 2*j+3], xticks=WilkinsonTicks(3), yticks=WilkinsonTicks(3)) 
        for i ∈ 1:9, j ∈ 1:2 
    ]

    for i ∈ 4:-1:1
        _tdf = filter(:chain => x -> x == i, maskingdatachain1)
        lines!(
            axs1[1, 1], _tdf.iteration, _tdf.log_density; 
            color=COLOURVECTOR[i], linewidth=1,
        )
        density!(
            axs1[1, 2], _tdf.log_density; 
            color=( :white, 0 ), strokecolor=COLOURVECTOR[i], strokewidth=1,
        )
        for (j, v) ∈ enumerate(names(_tdf)[3:10])
            lines!(
                axs1[1+j, 1], _tdf.iteration, getproperty(_tdf, v); 
                color=COLOURVECTOR[i], linewidth=1,
            )
            density!(
                axs1[1+j, 2], getproperty(_tdf, v); 
                color=( :white, 0 ), strokecolor=COLOURVECTOR[i], strokewidth=1,
            )
        end
        for (j, v) ∈ enumerate(names(_tdf)[11:19])
            lines!(
                axs2[j, 1], _tdf.iteration, getproperty(_tdf, v); 
                color=COLOURVECTOR[i], linewidth=1,
            )
            density!(
                axs2[j, 2], getproperty(_tdf, v); 
                color=( :white, 0 ), strokecolor=COLOURVECTOR[i], strokewidth=1,
            )
        end
    end

    for i ∈ 1:9, j ∈ 1:2
        formataxis!(axs1[i, j]; hidespines=( :r, :t ), trimspines=true,)
        #Label(fig[i, 1], "Iteration"; fontsize=11.84, tellwidth=false)
        Label(fig[i, 2], "Density"; fontsize=11.84, rotation=π/2, tellheight=false,)
        Label(fig[i, 0], _names1[i]; fontsize=11.84, rotation=π/2, tellheight=false,)
        i == 13 && continue
        formataxis!(axs2[i, j]; hidespines=( :r, :t ), trimspines=true,)
        #Label(fig[i, 5], "Iteration"; fontsize=11.84, tellwidth=false)
        Label(fig[i, 6], "Density"; fontsize=11.84, rotation=π/2, tellheight=false,)
        Label(fig[i, 4], _names2[i]; fontsize=11.84, rotation=π/2, tellheight=false,)
    end

    Label(fig[10, 1], "Iteration"; fontsize=11.84, tellwidth=false,)
    Label(fig[10, 5], "Iteration"; fontsize=11.84, tellwidth=false, )
    for c ∈ [ 1, 3, 5, 7 ] colgap!(fig.layout, c, 5) end
    rowgap!(fig.layout, 9, 5)

    fig
end

safesave(plotsdir("maskingdatachain1chainsplot.pdf"), maskingdatachain1chainsplot)

maskingdatafit1 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, maskingdatachain1, facialcoveringsrecommended[1:191, :]; 
    initialvalues=maskcovidcases[1:91, :], 
    Ns=POPULATION2020,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=[ 1.0; collect(56.0:28:191); 191 ],
)

maskdata1kv = keyvalues(maskingdatachain1, maskingdatafit1)
print(maskdata1kv)

subsetmaskdatafit1plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        maskcovidcases[1:191, :], 
        W_maskcoviddata[1:191, :], maskingdatafit1, 
        fitws(
            maskcovidcases[1:191, :], 
            POPULATION2020, 
            maskingdatafit1
        ), 
        1;
        markersize=2,
        xticks=( [ 1, 92, 183 ], [ "Jan.", "April", "July" ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_r0!(
        ga, maskcovidcases[1:191, :], maskingdatafit1, 2;
        plotcounterfactuals=true, 
        xticks=( [ 1, 92, 183 ], [ "Jan.", "April", "July" ] ), 
        ytitle=L"$\mathcal{R}_0$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases[1:191, :], POPULATION2020, maskingdatafit1, 3;
        markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
        fittedcolour=( COLOURVECTOR[2], 0.75 ), 
        xticks=( [ 1, 92, 183 ], [ "Jan.", "April", "July" ] ), 
        ytitle=L"Without \\ intervention$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases[1:191, :], POPULATION2020, maskingdatafit1, 4;
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 92, 183 ], [ "Jan.", "April", "July" ] ), 
        ytitle=L"With \\ intervention$$",
    )
    axs5 = plotrenewalequationsamples_causaleffect!(
        ga, maskcovidcases[1:191, :], nothing, POPULATION2020, maskingdatafit1, 5;
        cumulativedifference=true,
        fittedparameter=:y_matrix_det_vec,
        counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
        xticklabelrotation=-π/4,
        xticks=( [ 1, 92, 183 ], [ "Jan.", "April", "July" ] ), 
        xtitle="Date, 2020",
        ytitle=L"Cumulative \\ difference$$",
    )

    linkaxes!(axs3..., axs4...)

    for (i, ℓ) ∈ enumerate([ 
        "England", 
        "Northern Ireland", 
        "Scotland",  
        "Wales" 
    ])
        Label(
        ga[0, i], ℓ; 
        fontsize=10, halign=:left, tellwidth=false
    )
    end

    colgap!(ga, 1, 5)  
    for r ∈ [ 1, 6 ] rowgap!(ga, r, 5) end
    for axs ∈ [ axs1, axs2, axs3, axs4, axs5 ]
        if axs === axs5 
            formataxis!(
                axs[1]; 
                hidespines=( :r, :t ), trimspines=true,
            )
            for i ∈ 2:4
                formataxis!(
                    axs[i]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            end
        else
            formataxis!(
                axs[1]; 
                hidex=true, hidexticks=true, 
                hidespines=( :r, :t, :b ), trimspines=true,
            )
            for i ∈ 2:4 
                formataxis!(
                    axs[i]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    end

    for i ∈ 1:4 
        iax = Axis(ga[1:5, i]; xticks=( [ 1, 92, 183 ], [ "Jan.", "April", "July" ] ),)
        if i == 1
            vlines!(iax, 133; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        elseif i == 3
            vlines!(iax, 119; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        elseif i == 4
            vlines!(iax, 161; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        end
        formataxis!(
            iax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )
        iax.xgridstyle=( :dot, :dense ) 
        iax.xgridwidth = 1
        iax.xgridvisible = true
        linkxaxes!(iax, axs1[i])
    end

    fig
end

safesave(plotsdir("subsetmaskdatafit1plot.pdf"), subsetmaskdatafit1plot)


## Analysis 2 
# Effect of mask requirements. No other considerations of confounding 

maskingdatachain2 = loadanalysisdictsasdf("maskingdatamodel2", 8, maxrounds, 1120)
plotchains(maskingdatachain2)
maskingdatafit2 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, maskingdatachain2, facialcoveringsrequired; 
    initialvalues=maskcovidcases[1:91, :], 
    Ns=POPULATION2020,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=[ 1.0; collect(56.0:28:224); 257 ],
)

maskdata2kv = keyvalues(maskingdatachain2, maskingdatafit2)
print(maskdata2kv)

subsetmaskdatafit2plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        maskcovidcases, 
        W_maskcoviddata, maskingdatafit2, 
        fitws(
            maskcovidcases, 
            POPULATION2020, 
            maskingdatafit2
        ), 
        1;
        markersize=2,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_r0!(
        ga, maskcovidcases, maskingdatafit2, 2;
        plotcounterfactuals=true, 
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"$\mathcal{R}_0$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases, POPULATION2020, maskingdatafit2, 3;
        markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
        fittedcolour=( COLOURVECTOR[2], 0.75 ), 
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"Without \\ intervention$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases, POPULATION2020, maskingdatafit2, 4;
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"With \\ intervention$$",
    )
    axs5 = plotrenewalequationsamples_causaleffect!(
        ga, maskcovidcases, nothing, POPULATION2020, maskingdatafit2, 5;
        cumulativedifference=true,
        fittedparameter=:y_matrix_det_vec,
        counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
        xticklabelrotation=-π/4,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        xtitle="Date, 2020",
        ytitle=L"Cumulative \\ difference$$",
    )

    linkaxes!(axs3..., axs4...)

    for (i, ℓ) ∈ enumerate([ 
        "England", 
        "Northern Ireland", 
        "Scotland",  
        "Wales" 
    ])
        Label(
        ga[0, i], ℓ; 
        fontsize=10, halign=:left, tellwidth=false
    )
    end

    colgap!(ga, 1, 5)  
    for r ∈ [ 1, 6 ] rowgap!(ga, r, 5) end
    for axs ∈ [ axs1, axs2, axs3, axs4, axs5 ]
        if axs === axs5 
            formataxis!(
                axs[1]; 
                hidespines=( :r, :t ), trimspines=true,
            )
            for i ∈ 2:4
                formataxis!(
                    axs[i]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            end
        else
            formataxis!(
                axs[1]; 
                hidex=true, hidexticks=true, 
                hidespines=( :r, :t, :b ), trimspines=true,
            )
            for i ∈ 2:4 
                formataxis!(
                    axs[i]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    end

    for i ∈ 1:4 
        iax = Axis(ga[1:5, i]; xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ),)
        if i == 1
            vlines!(iax, 167; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        elseif i == 2
            vlines!(iax, 192; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        elseif i == 3
            vlines!(iax, 174; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        end
        formataxis!(
            iax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )
        iax.xgridstyle=( :dot, :dense ) 
        iax.xgridwidth = 1
        iax.xgridvisible = true
        linkxaxes!(iax, axs1[i])
    end

    fig
end

safesave(plotsdir("subsetmaskdatafit2plot.pdf"), subsetmaskdatafit2plot)



## Analysis 3 
# Effect of mask requirements with secondary interventions of end of stay-at-home and some
# businesses reopening

maskingdatachain3 = loadanalysisdictsasdf("maskingdatamodel3", 8, maxrounds, 1130)
plotchains(maskingdatachain3)
maskingdatafit3 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, maskingdatachain3, facialcoveringsrequired; 
    initialvalues=maskcovidcases[1:91, :], 
    Ns=POPULATION2020,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=[ 1.0; collect(56.0:28:224); 257 ],
    secondaryinterventions=[ endstayathometimes, somebusinessreopen ],
)

subsetmaskdatafit3plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        maskcovidcases, 
        W_maskcoviddata, maskingdatafit3, 
        fitws(
            maskcovidcases, 
            POPULATION2020, 
            maskingdatafit3
        ), 
        1;
        markersize=2,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_r0!(
        ga, maskcovidcases, maskingdatafit3, 2;
        plotcounterfactuals=true, 
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"$\mathcal{R}_0$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases, POPULATION2020, maskingdatafit3, 3;
        markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
        fittedcolour=( COLOURVECTOR[2], 0.75 ), 
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"Without \\ intervention$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases, POPULATION2020, maskingdatafit3, 4;
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"With \\ intervention$$",
    )
    axs5 = plotrenewalequationsamples_causaleffect!(
        ga, maskcovidcases, nothing, POPULATION2020, maskingdatafit3, 5;
        cumulativedifference=true,
        fittedparameter=:y_matrix_det_vec,
        counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
        xticklabelrotation=-π/4,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        xtitle="Date, 2020",
        ytitle=L"Cumulative \\ difference$$",
    )

    linkaxes!(axs3..., axs4...)

    for (i, ℓ) ∈ enumerate([ 
        "England", 
        "Northern Ireland", 
        "Scotland",  
        "Wales" 
    ])
        Label(
        ga[0, i], ℓ; 
        fontsize=10, halign=:left, tellwidth=false
    )
    end

    colgap!(ga, 1, 5)  
    for r ∈ [ 1, 6 ] rowgap!(ga, r, 5) end
    for axs ∈ [ axs1, axs2, axs3, axs4, axs5 ]
        if axs === axs5 
            formataxis!(
                axs[1]; 
                hidespines=( :r, :t ), trimspines=true,
            )
            for i ∈ 2:4
                formataxis!(
                    axs[i]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            end
        else
            formataxis!(
                axs[1]; 
                hidex=true, hidexticks=true, 
                hidespines=( :r, :t, :b ), trimspines=true,
            )
            for i ∈ 2:4 
                formataxis!(
                    axs[i]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    end

    for i ∈ 1:4 
        iax = Axis(ga[1:5, i]; xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ),)
        if i == 1
            vlines!(iax, 167; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        elseif i == 2
            vlines!(iax, 192; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        elseif i == 3
            vlines!(iax, 174; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        end
        formataxis!(
            iax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )
        iax.xgridstyle=( :dot, :dense ) 
        iax.xgridwidth = 1
        iax.xgridvisible = true
        linkxaxes!(iax, axs1[i])
    end

    fig
end

safesave(plotsdir("subsetmaskdatafit3plot.pdf"), subsetmaskdatafit3plot)



## Analysis 4 
# Add lead and lag times  

maskingdatachain4 = loadanalysisdictsasdf("maskingdatamodel4", 8, maxrounds, 1140)
plotchains(maskingdatachain4)
maskingdatafit4 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, maskingdatachain4, facialcoveringsrequired; 
    initialvalues=maskcovidcases[1:91, :], 
    Ns=POPULATION2020,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=[ 1.0; collect(56.0:28:224); 257 ],
    secondaryinterventions=[ 
        [ endstayathometimes, somebusinessreopen ]; secondaryinterventions_data 
    ],
)

subsetmaskdatafit4plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        maskcovidcases, 
        W_maskcoviddata, maskingdatafit4, 
        fitws(
            maskcovidcases, 
            POPULATION2020, 
            maskingdatafit4
        ), 
        1;
        markersize=2,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_r0!(
        ga, maskcovidcases, maskingdatafit4, 2;
        plotcounterfactuals=true, 
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"$\mathcal{R}_0$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases, POPULATION2020, maskingdatafit4, 3;
        markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
        fittedcolour=( COLOURVECTOR[2], 0.75 ), 
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"Without \\ intervention$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases, POPULATION2020, maskingdatafit4, 4;
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"With \\ intervention$$",
    )
    axs5 = plotrenewalequationsamples_causaleffect!(
        ga, maskcovidcases, nothing, POPULATION2020, maskingdatafit4, 5;
        cumulativedifference=true,
        fittedparameter=:y_matrix_det_vec,
        counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
        xticklabelrotation=-π/4,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        xtitle="Date, 2020",
        ytitle=L"Cumulative \\ difference$$",
    )

    linkaxes!(axs3..., axs4...)

    for (i, ℓ) ∈ enumerate([ 
        "England", 
        "Northern Ireland", 
        "Scotland",  
        "Wales" 
    ])
        Label(
        ga[0, i], ℓ; 
        fontsize=10, halign=:left, tellwidth=false
    )
    end

    colgap!(ga, 1, 5)  
    for r ∈ [ 1, 6 ] rowgap!(ga, r, 5) end
    for axs ∈ [ axs1, axs2, axs3, axs4, axs5 ]
        if axs === axs5 
            formataxis!(
                axs[1]; 
                hidespines=( :r, :t ), trimspines=true,
            )
            for i ∈ 2:4
                formataxis!(
                    axs[i]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            end
        else
            formataxis!(
                axs[1]; 
                hidex=true, hidexticks=true, 
                hidespines=( :r, :t, :b ), trimspines=true,
            )
            for i ∈ 2:4 
                formataxis!(
                    axs[i]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    end

    for i ∈ 1:4 
        iax = Axis(ga[1:5, i]; xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ),)
        if i == 1
            vlines!(iax, 167; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        elseif i == 2
            vlines!(iax, 192; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        elseif i == 3
            vlines!(iax, 174; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        end
        formataxis!(
            iax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )
        iax.xgridstyle=( :dot, :dense ) 
        iax.xgridwidth = 1
        iax.xgridvisible = true
        linkxaxes!(iax, axs1[i])
    end

    fig
end

safesave(plotsdir("subsetmaskdatafit4plot.pdf"), subsetmaskdatafit4plot)



## Analysis 5 

maskingdatachain5 = loadanalysisdictsasdf("maskingdatamodel5", 8, maxrounds, 1150)
plotchains(maskingdatachain5)
maskingdatafit5 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, maskingdatachain5, facialcoveringsrequired; 
    initialvalues=maskcovidcases[1:91, :], 
    Ns=POPULATION2020,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=[ 1.0; collect(56.0:28:224); 257 ],
    secondaryinterventions=[ endstayathometimes, somebusinessreopen, facialcoveringsrecommended ],
)

maskdata5kv = keyvalues(maskingdatachain5, maskingdatafit5)
print(maskdata5kv)

quantile(exp.(getproperty(maskingdatachain5, "logsecondarydelta1.logsecondarydelta")), [ 0.05, 0.5, 0.95 ])

quantile(exp.(getproperty(maskingdatachain5, "logsecondarydelta2.logsecondarydelta")), [ 0.05, 0.5, 0.95 ])

quantile(exp.(getproperty(maskingdatachain5, "logsecondarydelta3.logsecondarydelta")), [ 0.05, 0.5, 0.95 ])


subsetmaskdatafit5plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 450 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        maskcovidcases, 
        W_maskcoviddata, maskingdatafit5, 
        fitws(
            maskcovidcases, 
            POPULATION2020, 
            maskingdatafit5
        ), 
        1;
        markersize=2,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_r0!(
        ga, maskcovidcases, maskingdatafit5, 2;
        plotcounterfactuals=true, 
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"$\mathcal{R}_0$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases, POPULATION2020, maskingdatafit5, 3;
        markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
        fittedcolour=( COLOURVECTOR[2], 0.75 ), 
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"Without \\ intervention$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases, POPULATION2020, maskingdatafit5, 4;
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"With \\ intervention$$",
    )
    axs5 = plotrenewalequationsamples_causaleffect!(
        ga, maskcovidcases, nothing, POPULATION2020, maskingdatafit5, 5;
        cumulativedifference=true,
        fittedparameter=:y_matrix_det_vec,
        counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
        xticklabelrotation=-π/4,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        xtitle="Date, 2020",
        ytitle=L"Cumulative \\ difference$$",
    )

    linkaxes!(axs3..., axs4...)

    for (i, ℓ) ∈ enumerate([ 
        "England", 
        "Northern Ireland", 
        "Scotland",  
        "Wales" 
    ])
        Label(
        ga[0, i], ℓ; 
        fontsize=10, halign=:left, tellwidth=false
    )
    end

    colgap!(ga, 1, 5)  
    for r ∈ [ 1, 6 ] rowgap!(ga, r, 5) end
    for axs ∈ [ axs1, axs2, axs3, axs4, axs5 ]
        if axs === axs5 
            formataxis!(
                axs[1]; 
                hidespines=( :r, :t ), trimspines=true,
            )
            for i ∈ 2:4
                formataxis!(
                    axs[i]; 
                    hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t ), trimspines=true,
                )
            end
        else
            formataxis!(
                axs[1]; 
                hidex=true, hidexticks=true, 
                hidespines=( :r, :t, :b ), trimspines=true,
            )
            for i ∈ 2:4 
                formataxis!(
                    axs[i]; 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                    hidespines=( :l, :r, :t, :b ), trimspines=true,
                )
            end
        end
    end

    for i ∈ 1:4 
        iax = Axis(ga[1:5, i]; xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ),)
        if i == 1
            vlines!(iax, 167; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        elseif i == 2
            vlines!(iax, 192; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        elseif i == 3
            vlines!(iax, 174; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        end
        formataxis!(
            iax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b)
        )
        iax.xgridstyle=( :dot, :dense ) 
        iax.xgridwidth = 1
        iax.xgridvisible = true
        linkxaxes!(iax, axs1[i])
    end

    fig
end

safesave(plotsdir("subsetmaskdatafit5plot.pdf"), subsetmaskdatafit5plot)



