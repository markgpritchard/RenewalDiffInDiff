
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

    bandcolour = ( COLOURVECTOR[3], 0.3 )

    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 400 ))
        ga = GridLayout(fig[1, 1])
        gb = GridLayout(fig[1, 2])
        gc = GridLayout(fig[2, 1:2])

        let 
            ax = Axis(ga[1, 1])
            ax2 = Axis(ga[1, 2])

            band!(ax, ts[50:100], xs2cf[50:100], xs2[50:100]; color=bandcolour)
            lines!(ax, ts, xs1; color=:black, linewidth=1,)
            lines!(ax, ts, xs2; color=COLOURVECTOR[1], linewidth=1,)
            lines!(
                ax, ts[50:100], xs2cf[50:100]; 
                color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,
            )
                
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
        
            formataxis!(
                ax; 
                hidex=true, hidey=true, trimspines=true, hidespines=( :r, :t )
            )    
            formataxis!(
                ax2; 
                hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                hidespines=( :l, :r, :t, :b )  
            )
            
            Label(ga[1, 0], "y"; fontsize=11.84, tellheight=false)
            Label(ga[2, 1], "time"; fontsize=11.84, tellwidth=false)
        
            colgap!(ga, 1, 5)
            colgap!(ga, 2, 0)
            rowgap!(ga, 1, 5)
            colsize!(ga, 2, Auto(0.32))
            linkaxes!(ax, ax2)
        end

        let 
            u0_1 = [ 989.0, 6.0, 5.0, 0 ]
            u0_2 = [ 997.0, 2.0, 1.0, 0 ]
        
            prob1 = ODEProblem(sirmodel!, u0_1, ( 0.0, 50.0 ))
            sol1 = solve(prob1; saveat = 1)
            df1 = DataFrame(sol1)
        
            #incidence1 = 0.4 .* df1.value1 .* df1.value2 ./ 1000
            incidence1 = [ i == 1 ? df1.value4[i] : df1.value4[i] - df1.value4[i-1] for i ∈ axes(df1, 1) ]

            prob2 = ODEProblem(sirmodel!, u0_2, ( 0.0, 50.0 ))
            sol2 = solve(prob2; saveat = 1)
            df2 = DataFrame(sol2)
        
            #incidence2 = 0.4 .* df2.value1 .* df2.value2 ./ 1000
            incidence2 = [ i == 1 ? df2.value4[i] : df2.value4[i] - df2.value4[i-1] for i ∈ axes(df2, 1) ]

            sol1cf_constant25 = df1.value2[10] - df2.value2[10]
            sol1cf = df2.value2[10:end] .+ sol1cf_constant25
        
            sol1cf_logconstant25 = log(df1.value2[10]) - log(df2.value2[10])
            logsol1cf = log.(df2.value2[10:end]) .+ sol1cf_logconstant25
        
            sol1cf_cumulativeconstant25 = df1.value4[10] - df2.value4[10]
            cumulativesol1cf = df2.value4[10:end] .+ sol1cf_cumulativeconstant25
        
            sol1cf_cumulativelogconstant25 = log(df1.value4[10]) - log(df2.value4[10])
            cumulativelogsol1cf = log.(df2.value4[10:end]) .+ sol1cf_cumulativelogconstant25
        
            sol1cfr0_constant25 = df1.value1[10] .* 2 ./ 1000 - df2.value1[10] .* 2 ./ 1000
            sol1cfr0 = df2.value1[10:end] .* 2 ./ 1000 .+ sol1cfr0_constant25
        
            incidencecfconst25 = incidence1[10] - incidence2[10]
            incidencecf = incidence2[10:end] .+ incidencecfconst25

            logincidencecfconst25 = log(incidence1[10]) - log(incidence2[10])
            logincidencecf = log.(incidence2[10:end]) .+ logincidencecfconst25

            let 
                axs = [ Axis(gb[i, 1]; xticks=[ 0, 25, 50 ]) for i ∈ 1:2 ]

                band!(
                    axs[1], df1.timestamp[10:end], df1.value2[10:end], sol1cf; 
                    color=bandcolour
                )
                lines!(axs[1], df1.timestamp[10:end], sol1cf; color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,)
                lines!(axs[1], df1.timestamp, df1.value2; color=:black, linewidth=1,)
                lines!(axs[1], df2.timestamp, df2.value2; color=COLOURVECTOR[1], linewidth=1,)
                arrows!(axs[1], df1.timestamp[10:10], df1.value2[10:10] .+ 40, [ 0 ], [ -20 ])

                band!(
                    axs[2], df1.timestamp[10:end], incidence1[10:end], incidencecf; 
                    color=bandcolour
                )
                lines!(axs[2], df1.timestamp[10:end], incidencecf; color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,)
                lines!(axs[2], df1.timestamp[2:end], incidence1[2:end]; color=:black, linewidth=1,)
                lines!(axs[2], df2.timestamp[2:end], incidence2[2:end]; color=COLOURVECTOR[1], linewidth=1,)
                arrows!(axs[2], df1.timestamp[10:10], incidence1[10:10] .+ 18, [ 0 ], [ -9 ])
            
                formataxis!(axs[1]; hidex=true, hidexticks=true, trimspines=true, hidespines=( :r, :t, :b ))
                formataxis!(axs[2]; trimspines=true, hidespines=( :r, :t ))

                linkxaxes!(axs...)

                Label(
                    gb[1, 0], "prevalence"; 
                    fontsize=11.84, rotation=π/2, tellheight=false
                )
                Label(
                    gb[2, 0], "incidence"; 
                    fontsize=11.84, rotation=π/2, tellheight=false
                )
                Label(gb[3, 1], "time, days"; fontsize=11.84, tellwidth=false)
                colgap!(gb, 1, 5)
                rowgap!(gb, 2, 5)
            end

            let 
                axs = [ Axis(gc[1, 2*i]; xticks=[ 0, 25, 50 ]) for i ∈ 1:3 ]

                band!(
                    axs[1], df1.timestamp[10:end], log.(df1.value2[10:end]), logsol1cf; 
                    color=bandcolour
                )
                lines!(axs[1], df1.timestamp[10:end], logsol1cf; color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,)
                lines!(axs[1], df1.timestamp, log.(df1.value2); color=:black, linewidth=1,)
                lines!(axs[1], df2.timestamp, log.(df2.value2); color=COLOURVECTOR[1], linewidth=1,)
                arrows!(
                    axs[1], df1.timestamp[10:10], log.(df1.value2[10:10]) .+ 1.5, [ 0 ], [ -0.75 ]
                )

                band!(
                    axs[2], df1.timestamp[10:end], log.(incidence1[10:end]), logincidencecf; 
                    color=bandcolour
                )
                lines!(axs[2], df1.timestamp[10:end], logincidencecf; color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,)
                lines!(axs[2], df1.timestamp[2:end], log.(incidence1[2:end]); color=:black, linewidth=1,)
                lines!(axs[2], df2.timestamp[2:end], log.(incidence2[2:end]); color=COLOURVECTOR[1], linewidth=1,)
                arrows!(axs[2], df1.timestamp[10:10], log.(incidence1[10:10]) .+ 1.2, [ 0 ], [ -0.6 ])

                band!(
                    axs[3], df1.timestamp[10:end], log.(df1.value4[10:end]), cumulativelogsol1cf; 
                    color=bandcolour
                )
                lines!(axs[3], df1.timestamp[10:end], cumulativelogsol1cf; color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,)
                lines!(axs[3], df1.timestamp[2:end], log.(df1.value4[2:end]); color=:black, linewidth=1,)
                lines!(axs[3], df2.timestamp[2:end], log.(df2.value4[2:end]); color=COLOURVECTOR[1], linewidth=1,)
                arrows!(
                    axs[3], df1.timestamp[10:10], log.(df1.value4[10:10]) .+ 1.125, [ 0 ], [ -0.6 ]
                )

                for i ∈ 1:3 
                    formataxis!(axs[i]; trimspines=true, hidespines=( :r, :t ))
                    Label(gc[2, 2*i], "time, days"; fontsize=11.84, tellwidth=false)
                end
                linkxaxes!(axs...)

                Label(
                    gc[1, 1], "log prevalence"; 
                    fontsize=11.84, rotation=π/2, tellheight=false
                )
                Label(
                    gc[1, 3], "log incidence"; 
                    fontsize=11.84, rotation=π/2, tellheight=false
                )
                Label(
                    gc[1, 5], "log cumulative\nincidence"; 
                    fontsize=11.84, rotation=π/2, tellheight=false
                )

                for c ∈ [ 1, 3, 5 ] colgap!(gc, c, 5) end
                rowgap!(gc, 1, 5)
            end
        end
        
        
        labelplots!(
            [ "A", "B", "C" ], 
            [ ga, gb, gc ]; 
            cols=[ 0, 0, 1 ],
            rows=1,
            padding = ( 0, 0, 5, 0 )
        )

        colsize!(fig.layout, 2, Auto(0.5))
        rowsize!(fig.layout, 2, Auto(0.5))


        fig   
    end
    fig 
end

safesave(plotsdir("paralleltrendsfig.pdf"), paralleltrendsfig)



sirparraleltrendsfig = let 
    u0_1 = [ 989.0, 6.0, 5.0, 0 ]
    u0_2 = [ 997.0, 2.0, 1.0, 0 ]

    prob1 = ODEProblem(sirmodel!, u0_1, ( 0.0, 50.0 ))
    sol1 = solve(prob1; saveat = 1)
    df1 = DataFrame(sol1)

    incidence1 = 0.4 .* df1.value1 .* df1.value2 ./ 1000

    prob2 = ODEProblem(sirmodel!, u0_2, ( 0.0, 50.0 ))
    sol2 = solve(prob2; saveat = 1)
    df2 = DataFrame(sol2)

    incidence2 = 0.4 .* df2.value1 .* df2.value2 ./ 1000

    sol1cf_constant25 = df1.value2[10] - df2.value2[10]
    sol1cf = df2.value2[10:end] .+ sol1cf_constant25

    sol1cf_logconstant25 = log(df1.value2[10]) - log(df2.value2[10])
    logsol1cf = log.(df2.value2[10:end]) .+ sol1cf_logconstant25

    sol1cf_cumulativeconstant25 = df1.value4[10] - df2.value4[10]
    cumulativesol1cf = df2.value4[10:end] .+ sol1cf_cumulativeconstant25

    sol1cf_cumulativelogconstant25 = log(df1.value4[10]) - log(df2.value4[10])
    cumulativelogsol1cf = log.(df2.value4[10:end]) .+ sol1cf_cumulativelogconstant25

    sol1cfr0_constant25 = df1.value1[10] .* 2 ./ 1000 - df2.value1[10] .* 2 ./ 1000
    sol1cfr0 = df2.value1[10:end] .* 2 ./ 1000 .+ sol1cfr0_constant25

    incidencecfconst25 = incidence2[10] - incidence1[10]
    incidencecf = incidence2[10:end] .+ incidencecfconst25

    bandcolour = ( COLOURVECTOR[3], 0.3 )

    fig = with_theme(theme_latexfonts()) do 
        fig = Figure(; size=( 500, 250 ))
        axs = [ Axis(fig[1, i]) for i ∈ [ 1, 3, 5 ] ]
        axs2 = [ Axis(fig[2, i]) for i ∈ [ 1, 3, 5 ] ]
    
        band!(
            axs[1], df1.timestamp[10:end], df1.value2[10:end], sol1cf; 
            color=bandcolour
        )
        lines!(axs[1], df1.timestamp[10:end], sol1cf; color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,)
        lines!(axs[1], df1.timestamp, df1.value2; color=:black, linewidth=1,)
        lines!(axs[1], df2.timestamp, df2.value2; color=COLOURVECTOR[1], linewidth=1,)
        arrows!(axs[1], df1.timestamp[10:10], df1.value2[10:10] .+ 40, [ 0 ], [ -20 ])
    
        formataxis!(axs[1]; hidex=true, hidexticks=true, trimspines=true, hidespines=( :r, :t, :b ))
    
        band!(
            axs2[1], df1.timestamp[10:end], log.(df1.value2[10:end]), logsol1cf; 
            color=bandcolour
        )
        lines!(axs2[1], df1.timestamp[10:end], logsol1cf; color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,)
        lines!(axs2[1], df1.timestamp, log.(df1.value2); color=:black, linewidth=1,)
        lines!(axs2[1], df2.timestamp, log.(df2.value2); color=COLOURVECTOR[1], linewidth=1,)
        formataxis!(axs2[1]; trimspines=true, hidespines=( :r, :t ))
        arrows!(
            axs2[1], df1.timestamp[10:10], log.(df1.value2[10:10]) .+ 1.5, [ 0 ], [ -0.75 ]
        )
    
        band!(
            axs[2], df1.timestamp[10:end], df1.value4[10:end], cumulativesol1cf; 
            color=bandcolour
        )
        lines!(axs[2], df1.timestamp[10:end], cumulativesol1cf; color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,)
        lines!(axs[2], df1.timestamp, df1.value4; color=:black, linewidth=1,)
        lines!(axs[2], df2.timestamp, df2.value4; color=COLOURVECTOR[1], linewidth=1,)
        arrows!(axs[2], df1.timestamp[10:10], df1.value4[10:10] .+ 200, [ 0 ], [ -100 ])
    
        formataxis!(axs[2]; hidex=true, hidexticks=true, trimspines=true, hidespines=( :r, :t, :b ))
    
        band!(
            axs2[2], df1.timestamp[10:end], log.(df1.value4[10:end]), cumulativelogsol1cf; 
            color=bandcolour
        )
        lines!(axs2[2], df1.timestamp[10:end], cumulativelogsol1cf; color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,)
        lines!(axs2[2], df1.timestamp, log.(df1.value4); color=:black, linewidth=1,)
        lines!(axs2[2], df2.timestamp, log.(df2.value4); color=COLOURVECTOR[1], linewidth=1,)
        arrows!(
            axs2[2], df1.timestamp[10:10], log.(df1.value4[10:10]) .+ 1.125, [ 0 ], [ -0.6 ]
        )

        formataxis!(axs2[2]; trimspines=true, hidespines=( :r, :t ))
    
        band!(
            axs[3], df1.timestamp[10:end], df1.value1[10:end] .* 2 ./ 1000, sol1cfr0; 
            color=bandcolour
        )
        lines!(axs[3], df1.timestamp[10:end], sol1cfr0; color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,)
        lines!(axs[3], df1.timestamp, df1.value1 .* 2 ./ 1000; color=:black, linewidth=1,)
        lines!(axs[3], df2.timestamp, df2.value1 .* 2 ./ 1000; color=COLOURVECTOR[1], linewidth=1,)
        arrows!(
            axs[3], 
            df1.timestamp[10:10], 
            df2.value1[10:10] .* 2 ./ 1000 .+ 0.3, 
            [ 0 ], 
            [ -0.2 ]
        )
    
        formataxis!(axs[3]; hidex=true, hidexticks=true, trimspines=true, hidespines=( :r, :t, :b ))
    
        lines!(axs2[3], df1.timestamp, ones(length(df2.timestamp)) .* 2; color=:black, linewidth=1,)
        lines!(axs2[3], df2.timestamp, ones(length(df2.timestamp)) .* 2; color=COLOURVECTOR[1], linewidth=1,)
        arrows!(axs2[3], df1.timestamp[10:10], [ 2.3 ], [ 0 ], [ -0.2 ])
    
        formataxis!(axs2[3]; trimspines=true, hidespines=( :r, :t ))
    
        linkxaxes!(axs..., axs2...)
        linkyaxes!(axs[3], axs2[3])
        Label(
            fig[1, 0], "Number\ninfectious"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            fig[2, 0], "Log number\ninfectious"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            fig[1, 2], "Cumulative\ninfectious"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            fig[2, 2], "Log cumulative\ninfectious"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            fig[1, 4], "Effective\nreproduction ratio"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            fig[2, 4], "Basic reproduction\nratio"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
    
        for c ∈ [ 1, 3, 5 ] 
            Label(fig[3, c], "time"; fontsize=11.84, tellwidth=false)
            colgap!(fig.layout, c, 5) 
        end
        rowgap!(fig.layout, 2, 5)
    
        fig 
    end
    fig 
end

safesave(plotsdir("sirparraleltrendsfig.pdf"), sirparraleltrendsfig)















u0_1 = [ 989.0, 6.0, 5.0, 0 ]
            u0_2 = [ 997.0, 2.0, 1.0, 0 ]
        
            prob1 = ODEProblem(sirmodel!, u0_1, ( 0.0, 50.0 ))
            sol1 = solve(prob1; saveat = 1)
            df1 = DataFrame(sol1)
        
            incidence1 = 0.4 .* df1.value1 .* df1.value2 ./ 1000
        
            prob2 = ODEProblem(sirmodel!, u0_2, ( 0.0, 50.0 ))
            sol2 = solve(prob2; saveat = 1)
            df2 = DataFrame(sol2)
        
            incidence2 = 0.4 .* df2.value1 .* df2.value2 ./ 1000
        
            sol1cf_constant25 = df1.value2[10] - df2.value2[10]
            sol1cf = df2.value2[10:end] .+ sol1cf_constant25
        
            sol1cf_logconstant25 = log(df1.value2[10]) - log(df2.value2[10])
            logsol1cf = log.(df2.value2[10:end]) .+ sol1cf_logconstant25
        
            sol1cf_cumulativeconstant25 = df1.value4[10] - df2.value4[10]
            cumulativesol1cf = df2.value4[10:end] .+ sol1cf_cumulativeconstant25
        
            sol1cf_cumulativelogconstant25 = log(df1.value4[10]) - log(df2.value4[10])
            cumulativelogsol1cf = log.(df2.value4[10:end]) .+ sol1cf_cumulativelogconstant25
        
            sol1cfr0_constant25 = df1.value1[10] .* 2 ./ 1000 - df2.value1[10] .* 2 ./ 1000
            sol1cfr0 = df2.value1[10:end] .* 2 ./ 1000 .+ sol1cfr0_constant25
        
            incidencecfconst25 = incidence2[10] - incidence1[10]
            incidencecf = incidence2[10:end] .+ incidencecfconst25