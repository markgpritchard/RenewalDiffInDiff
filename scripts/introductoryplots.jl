
using DrWatson
@quickactivate :RenewalDiffInDiff
include(srcdir("plottingfunctions.jl"))

using CairoMakie, DifferentialEquations 

function sirmodel!(du, u, p, t)
    du[1] = -0.4 * u[1] * u[2] / sum(u)
    du[2] = 0.4 * u[1] * u[2] / sum(u) - 0.2 * u[2]
    du[3] = 0.2 * u[2]
    du[4] = 0.4 * u[1] * u[2] / sum(u)  # cumulative infections
end

paralleltrendsfig = let 
    u0_1 = [ 989.0, 6.0, 5.0, 0 ]
    u0_2 = [ 997.0, 2.0, 1.0, 0 ]

    prob1 = ODEProblem(sirmodel!, u0_1, ( 0.0, 50.0 ))
    sol1 = solve(prob1; saveat = 1)
    df1 = DataFrame(sol1)
    incidence1 = [ 
        i == 1 ? 
            df1.value4[i] : 
            df1.value4[i] - df1.value4[i-1] 
        for i ∈ axes(df1, 1) 
    ]

    prob2 = ODEProblem(sirmodel!, u0_2, ( 0.0, 50.0 ))
    sol2 = solve(prob2; saveat = 1)
    df2 = DataFrame(sol2)
    incidence2 = [ 
        i == 1 ? 
            df2.value4[i] : 
            df2.value4[i] - df2.value4[i-1] 
        for i ∈ axes(df2, 1) 
    ]

    sol1cf_constant25 = df1.value2[10] - df2.value2[10]
    sol1cf = df2.value2[10:end] .+ sol1cf_constant25

    incidencecfconst25 = incidence1[10] - incidence2[10]
    incidencecf = incidence2[10:end] .+ incidencecfconst25

    fig = with_theme(theme_latexfonts()) do
        bandcolour = ( COLOURVECTOR[3], 0.3 )

        fig = Figure(; size=( 200, 350 ))
        ga = GridLayout(fig[1, 1])
        axs = [ Axis(ga[i, 1]; xticks=[ 0, 25, 50 ]) for i ∈ 1:2 ]

        band!(
            axs[1], df1.timestamp[10:end], df1.value2[10:end], sol1cf; 
            color=bandcolour
        )
        lines!(
            axs[1], df1.timestamp[10:end], sol1cf; 
            color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,
        )
        lines!(axs[1], df1.timestamp, df1.value2; color=COLOURVECTOR[1], linewidth=1)
        lines!(axs[1], df2.timestamp, df2.value2; color=:black, linewidth=1,)
        arrows!(axs[1], df1.timestamp[10:10], df1.value2[10:10] .+ 40, [ 0 ], [ -20 ])

        band!(
            axs[2], df1.timestamp[10:end], incidence1[10:end], incidencecf; 
            color=bandcolour
        )
        lines!(
            axs[2], df1.timestamp[10:end], incidencecf; 
            color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,
        )
        lines!(
            axs[2], df1.timestamp[2:end], incidence1[2:end]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[2], df2.timestamp[2:end], incidence2[2:end]; 
            color=:black, linewidth=1,
        )
        arrows!(axs[2], df1.timestamp[10:10], incidence1[10:10] .+ 8, [ 0 ], [ -4.2 ])
    
        formataxis!(
            axs[1]; 
            hidex=true, hidexticks=true, trimspines=true, hidespines=( :r, :t, :b )
        )
        formataxis!(axs[2]; trimspines=true, hidespines=( :r, :t ))

        linkxaxes!(axs...)

        Label(
            ga[1, 0], "prevalence"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            ga[2, 0], "incidence"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(ga[3, 1], "time, days"; fontsize=11.84, tellwidth=false)
        colgap!(ga, 1, 5)
        rowgap!(ga, 2, 5)
    end
    fig 
end

safesave(plotsdir("paralleltrendsfig.pdf"), paralleltrendsfig)
