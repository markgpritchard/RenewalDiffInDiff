
using DrWatson
@quickactivate :RenewalDiffInDiff
include(srcdir("plottingfunctions.jl"))

using CairoMakie, DifferentialEquations 


parraleltrendsfig = let 
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
        
    fig = Figure(; size=( 400, 250 ))
    axs = [ Axis(fig[1, i]) for i ∈ 1:2 ]
    
    for ax ∈ axs 
        lines!(ax, ts, xs1; color=:blue)
        lines!(ax, ts, xs2; color=:black)
        arrows!(ax, [ 50 ], xs2[50:50] .+ 0.075, [ 0 ], [ -0.05 ])
        text!(
            ax, 50, xs2[50] + 0.075; 
            text="Group 1\nIntervention", align=( :center, :bottom ), fontsize=10
        )
    end

    text!(axs[1], 105, xs1[100]; text="Group 0", align = ( :right, :bottom ), fontsize=10)
    text!(axs[1], 105, xs2[90] - 0.01; text="Group 1", align = ( :right, :top ), fontsize=10)
    formataxis!(axs[1]; hidex=true, hidey=true)

    lines!(axs[2], ts[50:100], xs2cf[50:100]; color=:red, linestyle=:dash)
    band!(axs[2], ts[50:100], xs2cf[50:100], xs2[50:100]; color=( :red, 0.3 ))
    text!(
        axs[2], 110, xs2cf[100] + 0.03; 
        text="Counterfactual", align = ( :right, :top ), fontsize=10
    )
    text!(
        axs[2], 110, xs2[90] - 0.01; 
        text="Observed", align = ( :right, :top ), fontsize=10
    )      
    formataxis!(axs[2]; hidex=true, hidey=true, hideyticks=true, hidespines=( :l, :t, :r )) 
    
    linkaxes!(axs...)
    Label(fig[1, 0], "y"; fontsize=11.84, tellheight=false)
    Label(fig[2, 1:2], "time"; fontsize=11.84, tellwidth=false)

    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)
    
    fig    
end

safesave(plotsdir("parraleltrendsfig.svg"), parraleltrendsfig)

function sirmodel!(du, u, p, t)
    du[1] = -0.4 * u[1] * u[2] / sum(u)
    du[2] = 0.4 * u[1] * u[2] / sum(u) - 0.2 * u[2]
    du[3] = 0.2 * u[2]
end

sirparraleltrendsfig, R0_sirparraleltrendsfig = let 
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

    fig = Figure(; size=( 400, 350 ))
    axs = [ Axis(fig[1, i]) for i ∈ 1:2 ]
    axs2 = [ Axis(fig[2, i]) for i ∈ 1:2 ]

    for ax ∈ axs 
        lines!(ax, df1.timestamp, df1.value2; color=:black)
        lines!(ax, df2.timestamp, df2.value2; color=:blue)
        formataxis!(ax; hidex=true, hidey=true)
    end

    vlines!(axs[2], df1.timestamp[10]; color=:black, linestyle=:dot)

    lines!(axs[2], df1.timestamp[10:end], sol1cf; color=:red, linestyle=:dash)
    band!(axs[2], df1.timestamp[10:end], df1.value2[10:end], sol1cf; color=( :red, 0.3 ))

    for ax ∈ axs2 
        lines!(ax, df1.timestamp, log.(df1.value2); color=:black)
        lines!(ax, df2.timestamp, log.(df2.value2); color=:blue)
        formataxis!(ax; hidex=true, hidey=true)
    end

    vlines!(axs2[2], df1.timestamp[10]; color=:black, linestyle=:dot)

    lines!(axs2[2], df1.timestamp[10:end], logsol1cf; color=:red, linestyle=:dash)
    band!(axs2[2], df1.timestamp[10:end], log.(df1.value2[10:end]), logsol1cf; color=( :red, 0.3 ))

    linkaxes!(axs...)
    linkaxes!(axs2...)
    Label(fig[1, 0], "Number infectious"; fontsize = 11.84, rotation=π/2, tellheight = false)
    Label(fig[2, 0], "Log number infectious"; fontsize = 11.84, rotation=π/2, tellheight = false)
    Label(fig[3, 1:2], "time"; fontsize = 11.84, tellwidth = false)

    fig2 = Figure(; size=( 400, 350 ))
    axs = [ Axis(fig2[1, i]) for i ∈ 1:2 ]
    axs2 = [ Axis(fig2[2, i]) for i ∈ 1:2 ]

    for ax ∈ axs 
        lines!(ax, df1.timestamp, df1.value1 .* 2 ./ 1000; color=:black)
        lines!(ax, df2.timestamp, df2.value1 .* 2 ./ 1000; color=:blue)
        formataxis!(ax; hidex=true, hidey=true)
    end

    vlines!(axs[2], df1.timestamp[10]; color=:black, linestyle=:dot)

    lines!(axs[2], df1.timestamp[10:end], sol1cfr0; color=:red, linestyle=:dash)
    band!(axs[2], df1.timestamp[10:end], df1.value1[10:end] .* 2 ./ 1000, sol1cfr0; color=( :red, 0.3 ))

    for ax ∈ axs2 
        lines!(ax, df1.timestamp, ones(length(df2.timestamp)) .* 2; color=:black)
        lines!(ax, df2.timestamp, ones(length(df2.timestamp)) .* 2; color=:blue)
        formataxis!(ax; hidex=true, hidey=true)
    end

    vlines!(axs2[2], df1.timestamp[10]; color=:black, linestyle=:dot)

   # lines!(axs2[2], df1.timestamp[10:end], logsol1cf; color=:red, linestyle=:dash)
   # band!(axs2[2], df1.timestamp[10:end], log.(df1.value2[10:end]), logsol1cf; color=( :red, 0.3 ))

    linkaxes!(axs..., axs2...)
    Label(fig2[1, 0], "Effective reproduction\nnumber"; fontsize = 11.84, rotation=π/2, tellheight = false)
    Label(fig2[2, 0], "Basic reproduction\nnumber"; fontsize = 11.84, rotation=π/2, tellheight = false)
    Label(fig2[3, 1:2], "time"; fontsize = 11.84, tellwidth = false)


    ( fig, fig2 )    
end

safesave(plotsdir("sirparraleltrendsfig.svg"), sirparraleltrendsfig)
safesave(plotsdir("R0_sirparraleltrendsfig.svg"), R0_sirparraleltrendsfig)

