
using DrWatson
@quickactivate :RenewalDiffInDiff
include(srcdir("plottingfunctions.jl"))

using CairoMakie, DifferentialEquations 


paralleltrendsfig = let 
    ts = 1.0:1.0:100.0

    xs1 = [ 
        0.3 - 0.1 * (t - 50) / 25 + 0.7 * (t / 100)^1.4 + cos(2π * t / 25) * 3 / (t + 100) 
        for t ∈ ts 
    ]
    xs2 = [ 
        t < 50 ? 
            xs1[t] + 0.05 :
            (1 - (t - 50) / 250) * xs1[t] + 0.05
        for t ∈ round.(Int, ts)
    ]
    xs2cf = xs1 .+ 0.05
        
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 300 ))
        ax = Axis(fig[1, 1])
        ax2 = Axis(fig[1, 2])
        
        lines!(ax, ts, xs1; color=:black)
        lines!(ax, ts, xs2; color=COLOURVECTOR[1])
        lines!(ax, ts[50:100], xs2cf[50:100]; color=:red, linestyle=:dash)
        band!(ax, ts[50:100], xs2cf[50:100], xs2[50:100]; color=( :red, 0.3 ))
            
        arrows!(ax, [ 50 ], xs2[50:50] .+ 0.075, [ 0 ], [ -0.05 ])
    
        text!(
            ax, 50, xs2[50] + 0.075; 
            text="Group 1\nIntervention", align=( :center, :bottom ), fontsize=10
        )
        text!(
            ax2, 0, xs1[100]; 
            text="Group 0\nobserved", align = ( :left, :top ), fontsize=10
        )
        text!(
            ax2, 0, xs2[100]; 
            text="Group 1\nobserved", align = ( :left, :center ), fontsize=10
        )
        text!(
            ax2, 0, xs2cf[100] + 0.005; 
            text="Group 1\ncounterfactual", align = ( :left, :top ), fontsize=10
        )
    
        formataxis!(ax; hidex=true, hidey=true)    
        formataxis!(ax2; hidex=true, hidexticks=true, hidey=true, hideyticks=true)
        hidespines!(ax2, :l, :r, :t, :b )    
        
        Label(fig[1, 0], "y"; fontsize=11.84, tellheight=false)
        Label(fig[2, 1], "time"; fontsize=11.84, tellwidth=false)
    
        colgap!(fig.layout, 1, 5)
        colgap!(fig.layout, 2, 0)
        rowgap!(fig.layout, 1, 5)
        colsize!(fig.layout, 2, Auto(0.2))
        linkaxes!(ax, ax2)
        
        fig   
    end
    fig 
end

safesave(plotsdir("paralleltrendsfig.pdf"), paralleltrendsfig)

function sirmodel!(du, u, p, t)
    du[1] = -0.4 * u[1] * u[2] / sum(u)
    du[2] = 0.4 * u[1] * u[2] / sum(u) - 0.2 * u[2]
    du[3] = 0.2 * u[2]
end

sirparraleltrendsfig = let 
    u0_1 = [ 989, 6, 5 ]
    u0_2 = [ 997, 2, 1 ]

    prob1 = ODEProblem(sirmodel!, u0_1, ( 0.0, 50.0 ))
    sol1 = solve(prob1; saveat = 1)
    df1 = DataFrame(sol1)

    prob2 = ODEProblem(sirmodel!, u0_2, ( 0.0, 50.0 ))
    sol2 = solve(prob2; saveat = 1)
    df2 = DataFrame(sol2)

    sol1cf_constant25 = df1.value2[10] - df2.value2[10]
    sol1cf = df2.value2[10:end] .+ sol1cf_constant25

    sol1cf_logconstant25 = log(df1.value2[10]) - log(df2.value2[10])
    logsol1cf = log.(df2.value2[10:end]) .+ sol1cf_logconstant25

    sol1cfr0_constant25 = df1.value1[10] .* 2 ./ 1000 - df2.value1[10] .* 2 ./ 1000
    sol1cfr0 = df2.value1[10:end] .* 2 ./ 1000 .+ sol1cfr0_constant25

    fig = Figure(; size=( 500, 350 ))
    axs = [ Axis(fig[1, i]) for i ∈ [ 1, 3 ] ]
    axs2 = [ Axis(fig[2, i]) for i ∈ [ 1, 3 ] ]

    lines!(axs[1], df1.timestamp, df1.value2; color=:black)
    lines!(axs[1], df2.timestamp, df2.value2; color=:blue)
    lines!(axs[1], df1.timestamp[10:end], sol1cf; color=:red, linestyle=:dash)
    band!(axs[1], df1.timestamp[10:end], df1.value2[10:end], sol1cf; color=( :red, 0.3 ))
    arrows!(axs[1], df1.timestamp[10:10], df1.value2[10:10] .+ 40, [ 0 ], [ -20 ])

    formataxis!(axs[1]; hidex=true, hidexticks=true, hidespines=( :b, :t, :r ))

    lines!(axs2[1], df1.timestamp, log.(df1.value2); color=:black)
    lines!(axs2[1], df2.timestamp, log.(df2.value2); color=:blue)
    lines!(axs2[1], df1.timestamp[10:end], logsol1cf; color=:red, linestyle=:dash)
    band!(
        axs2[1], df1.timestamp[10:end], log.(df1.value2[10:end]), logsol1cf; 
        color=( :red, 0.3 )
    )
    formataxis!(axs2[1])
    arrows!(axs2[1], df1.timestamp[10:10], log.(df1.value2[10:10]) .+ 1.5, [ 0 ], [ -0.75 ])

    lines!(axs[2], df1.timestamp, df1.value1 .* 2 ./ 1000; color=:black)
    lines!(axs[2], df2.timestamp, df2.value1 .* 2 ./ 1000; color=:blue)
    lines!(axs[2], df1.timestamp[10:end], sol1cfr0; color=:red, linestyle=:dash)
    band!(
        axs[2], df1.timestamp[10:end], df1.value1[10:end] .* 2 ./ 1000, sol1cfr0; 
        color=( :red, 0.3 )
    )
    arrows!(
        axs[2], df1.timestamp[10:10], df2.value1[10:10] .* 2 ./ 1000 .+ 0.3, [ 0 ], [ -0.2 ]
    )

    formataxis!(axs[2]; hidex=true, hidexticks=true, hidespines=( :b, :t, :r ))

    lines!(axs2[2], df1.timestamp, ones(length(df2.timestamp)) .* 2; color=:black)
    lines!(axs2[2], df2.timestamp, ones(length(df2.timestamp)) .* 2; color=:blue)
    arrows!(axs2[2], df1.timestamp[10:10], [ 2.3 ], [ 0 ], [ -0.2 ])

    formataxis!(axs2[2]; )

    linkxaxes!(axs..., axs2...)
    linkyaxes!(axs[2], axs2[2])
    Label(fig[1, 0], "Number infectious"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(fig[2, 0], "Log number infectious"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(
        fig[1, 2], "Effective reproduction\nnumber"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    Label(
        fig[2, 2], "Basic reproduction\nnumber"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    Label(fig[3, 1:3], "time"; fontsize=11.84, tellwidth=false)

    for c ∈ [ 1, 3 ] colgap!(fig.layout, c, 5) end
    rowgap!(fig.layout, 2, 5)

    fig 
end

safesave(plotsdir("sirparraleltrendsfig.pdf"), sirparraleltrendsfig)

