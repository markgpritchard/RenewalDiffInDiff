
using DrWatson
@quickactivate :RenewalDiffInDiff
include("loadanalysisforplots.jl")
include(srcdir("plottingfunctions.jl"))
using CairoMakie, DifferentialEquations, StatsBase


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Introductory plot / figure 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

paralleltrendsfig = let 
    bandcolour = ( COLOURVECTOR[3], 0.3 )

    function _sirmodel!(du, u, p, t)
        du[1] = -0.4 * u[1] * u[2] / sum(u)
        du[2] = 0.4 * u[1] * u[2] / sum(u) - 0.2 * u[2]
        du[3] = 0.2 * u[2]
        du[4] = 0.4 * u[1] * u[2] / sum(u)  # cumulative infections
    end

    function _sirmodeldf(u0, tspan=( 0.0, 50.0 ); saveat=1)
        prob = ODEProblem(_sirmodel!, u0, tspan)
        sol = solve(prob; saveat)
        df = DataFrame(sol)
        return df 
    end

    function _modelincidence(cumi)
        return [ i == 1 ? cumi[i] : cumi[i] - cumi[i-1] for i ∈ eachindex(cumi) ]
    end

    u0_1 = [ 989.0, 6.0, 5.0, 0 ]
    u0_2 = [ 997.0, 2.0, 1.0, 0 ]

    df1 = _sirmodeldf(u0_1)
    incidence1 = _modelincidence(df1.value4)

    df2 = _sirmodeldf(u0_2)
    incidence2 = _modelincidence(df2.value4)

    sol1cf_constant25 = df1.value2[10] - df2.value2[10]
    sol1cf = df2.value2[10:end] .+ sol1cf_constant25

    incidencecfconst25 = incidence1[10] - incidence2[10]
    incidencecf = incidence2[10:end] .+ incidencecfconst25

    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 200, 350 ))
        ga = GridLayout(fig[1, 1])
        axs = [ Axis(ga[i, 1]; xticks=[ 0, 25, 50 ]) for i ∈ 1:2 ]

        band!(
            axs[1], df1.timestamp[10:end], df1.value2[10:end], sol1cf; 
            color=bandcolour,
        )
        lines!(
            axs[1], df1.timestamp[10:end], sol1cf; 
            color=COLOURVECTOR[3], linestyle=:dash, linewidth=1,
        )
        lines!(
            axs[1], df1.timestamp, df1.value2; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[1], df2.timestamp, df2.value2; 
            color=:black, linewidth=1,
        )
        band!(
            axs[2], df1.timestamp[10:end], incidence1[10:end], incidencecf; 
            color=bandcolour,
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
    
        formataxis!(
            axs[1]; 
            hidex=true, hidexticks=true, trimspines=true, hidespines=( :r, :t, :b ),
        )
        formataxis!(
            axs[2]; 
            trimspines=true, hidespines=( :r, :t ), setpoint=30,
        )

        linkxaxes!(axs...)

        Label(
            ga[1, 0], "prevalence, per 1000"; 
            fontsize=11.84, rotation=π/2, tellheight=false,
        )
        Label(
            ga[2, 0], "daily incidence, per 1000"; 
            fontsize=11.84, rotation=π/2, tellheight=false,
        )
        Label(
            ga[3, 1], "time, days"; 
            fontsize=11.84, tellwidth=false,
        )
        colgap!(ga, 1, 5)
        rowgap!(ga, 2, 5)

        fig
    end
    fig 
end
safesave(plotsdir("paralleltrendsfig.pdf"), paralleltrendsfig)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Results for simulation 1 / figure 2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim1_0_fitzfig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 500 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    attax1 = Axis(ga[3, 1:2])
    casesax1 = Axis(ga[5, 1:2])
    attax2 = Axis(gb[3, 1:2])
    casesax2 = Axis(gb[5, 1:2])

    let 
        axs = [ Axis(ga[1, i]) for i ∈ 1:2 ]

        for i ∈ 1:2 
            lines!(
                axs[i], 1:100, sim1nointerventionlogeffectivereproductionratios.med[:, i];
                color=COLOURVECTOR[1], linewidth=1,
            )
            band!(
                axs[i], 
                1:100, 
                sim1nointerventionlogeffectivereproductionratios.lci[:, i],
                sim1nointerventionlogeffectivereproductionratios.uci[:, i];
                color=( COLOURVECTOR[1], 0.5 ), 
            )
            noninfindex = findall(x -> -Inf < x < Inf, W_sim1_0[:, i])
            scatter!(
                axs[i], collect(1:100)[noninfindex], W_sim1_0[noninfindex, i];
                color=:black, markersize=2,
            )
            linkaxes!(axs...)
            formataxis!(
                axs[i]; 
                hidespines=( :r, :t ), hidey=(i != 1), hideyticks=(i != 1), trimspines=true,
            )
            if i == 2 hidespines!(axs[i], :l) end
        end

        lines!(
            attax1, [ -21, 21 ], [ 0, 0 ]; 
            color=COLOURVECTOR[2], linestyle=( :dot, :dense ), linewidth=1,
        )
        plotphi_att!(attax1, sim1leadlagnointerventiondf, -21:7:21; )

        lines!(
            casesax1, 
            1:100, 
            sim1nointerventioncasesdiff.med[:, 2] .* 10_000 ./ simulation1dataset["Ns"][2];
            color=COLOURVECTOR[1], linewidth=1,
        )
        band!(
            casesax1, 
            1:100, 
            sim1nointerventioncasesdiff.lci[:, 2] .* 10_000 ./ simulation1dataset["Ns"][2],
            sim1nointerventioncasesdiff.uci[:, 2] .* 10_000 ./ simulation1dataset["Ns"][2];
            color=( COLOURVECTOR[1], 0.5 ), 
        )
        hlines!(
            casesax1, 0; 
            color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
        )
        vlines!(casesax1, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)

        formataxis!(casesax1; hidespines=( :r, :t ), trimspines=true,)
    end

    let 
        axs = [ Axis(gb[1, i]) for i ∈ 1:2 ]

        for i ∈ 1:2 
            lines!(
                axs[i], 1:100, sim1logeffectivereproductionratios.med[:, i];
                color=COLOURVECTOR[1], linewidth=1,
            )
            band!(
                axs[i], 
                1:100, 
                sim1logeffectivereproductionratios.lci[:, i],
                sim1logeffectivereproductionratios.uci[:, i];
                color=( COLOURVECTOR[1], 0.5 ), 
            )
            noninfindex = findall(x -> -Inf < x < Inf, W_sim1[:, i])
            scatter!(
                axs[i], collect(1:100)[noninfindex], W_sim1[noninfindex, i];
                color=:black, markersize=2,
            )

            formataxis!(
                axs[i]; 
                hidespines=( :r, :t ), hidey=(i != 1), hideyticks=(i != 1), trimspines=true,
            )
            if i == 2 hidespines!(axs[i], :l) end
        end
        
        lines!(
            attax2, [ -21, 0, 0, 21 ], [ 0, 0, log(0.8), log(0.8) ]; 
            color=COLOURVECTOR[2], linestyle=( :dot, :dense ), linewidth=1,
        )
        plotphi_att!(attax2, sim1model1lagleaddf, -21:7:21; )

        lines!(
            casesax2, 
            1:100, 
            sim1casesdiff.med[:, 2] .* 10_000 ./ simulation1dataset["Ns"][2];
            color=COLOURVECTOR[1], linewidth=1,
        )
        band!(
            casesax2, 
            1:100, 
            sim1casesdiff.lci[:, 2] .* 10_000 ./ simulation1dataset["Ns"][2],
            sim1casesdiff.uci[:, 2] .* 10_000 ./ simulation1dataset["Ns"][2];
            color=( COLOURVECTOR[1], 0.5 ), 
        )
        hlines!(
            casesax2, 0; 
            color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
        )
        vlines!(casesax2, 50; color=:red, linestyle=( :dot, :dense ), linewidth=1,)

        formataxis!(casesax2; hidespines=( :r, :t ), trimspines=true,)
    end

    for gl ∈ [ ga, gb ]
        Label(
            gl[1, 0], L"log $\mathcal{R}_{e}$"; 
            fontsize=11.84, rotation=π/2, tellheight=false,
        )
        Label(
            gl[3, 0], L"$\varphi$"; 
            fontsize=11.84, rotation=π/2, tellheight=false,
        )
        Label(
            gl[5, 0], L"cumulative difference \\ in diagnoses, per $10\,000$"; 
            fontsize=11.84, rotation=π/2, tellheight=false,
        )
        Label(gl[2, 1:2], "Time, days"; fontsize=11.84, tellwidth=false,)
        Label(gl[4, 1:2], "Time, days from intervention"; fontsize=11.84, tellwidth=false,)
        Label(gl[6, 1:2], "Time, days"; fontsize=11.84, tellwidth=false,)
        colgap!(gl, 1, 5) 
        for r ∈ [ 1, 3, 5 ] rowgap!(gl, r, 5) end
    end
    
    linkaxes!(attax1, attax2)
    linkaxes!(casesax1, casesax2)

    labelplots!([ "A", "B" ], [ ga, gb ]; rows=1)
    fig
end
safesave(plotsdir("sim1_0_fitzfig.pdf"), sim1_0_fitzfig)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Causal effect from other simulations / figure 3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

simfig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 250, 450 ))
    ga = GridLayout(fig[1, 1])
    axs = [ Axis(ga[(2 * i - 1), j]; xticks=[ -21, 0, 21 ]) for i ∈ 1:3, j ∈ 1:2 ]

    lines!(
        axs[1, 1], [ -21, 21 ], [ 0, 0 ]; 
        color=COLOURVECTOR[2], linestyle=( :dot, :dense ), linewidth=1,
    )
    plotphi_att!(axs[1, 1], sim2leadlagconfoundernointerventiondf, -21:7:21; )
    lines!(
        axs[1, 2], [ -21, 0, 0, 21 ], [ 0, 0, log(0.8), log(0.8) ]; 
        color=COLOURVECTOR[2], linestyle=( :dot, :dense ), linewidth=1,
    )
    plotphi_att!(axs[1, 2], sim2model1leadlagdf, -21:7:21; )
    lines!(
        axs[2, 1], [ -21, 21 ], [ 0, 0 ]; 
        color=COLOURVECTOR[2], linestyle=( :dot, :dense ), linewidth=1,
    )
    plotphi_att!(axs[2, 1], sim3leadlagnointerventiondf, -21:7:21; )
    lines!(
        axs[2, 2], [ -21, 0, 0, 21 ], [ 0, 0, log(0.8), log(0.8) ]; 
        color=COLOURVECTOR[2], linestyle=( :dot, :dense ), linewidth=1,
    )
    plotphi_att!(axs[2, 2], sim3model1leadlagdf, -21:7:21; )
    lines!(
        axs[3, 1], [ -21, 21 ], [ 0, 0 ]; 
        color=COLOURVECTOR[2], linestyle=( :dot, :dense ), linewidth=1,
    )
    plotphi_att!(axs[3, 1], sim4leadlagnointerventiondf, -21:7:21; )
    lines!(
        axs[3, 2], [ -21, 0, 0, 21 ], [ 0, 0, log(0.8), log(0.8) ]; 
        color=COLOURVECTOR[2], linestyle=( :dot, :dense ), linewidth=1,
    )
    plotphi_att!(axs[3, 2], sim4model1leadlagdf, -21:7:21; )

    for row ∈ [ 1, 3, 5 ] 
        Label(
            ga[row, 0], L"$\varphi$"; 
            fontsize=11.84, rotation=π/2, tellheight=false,
        )
        Label(
            ga[row+1, 1:2], "Time, days from intervention"; 
            fontsize=11.84, tellwidth=false,
        )
    end

    linkaxes!(axs...)
    colgap!(ga, 1, 5) 
    for r ∈ [ 1, 3, 5 ] rowgap!(ga, r, 5) end 

    labelplots!([ "A", "B", "C" ], ga; rows=[ 1, 3, 5 ])
    fig
end
safesave(plotsdir("simfig.pdf"), simfig)

exp.(quantile(sim2leadlagconfoundernointerventiondf.phi_plus21, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(sim2model1leadlagdf.phi_plus21, [ 0.05, 0.5, 0.95 ]))

exp.(quantile(sim3leadlagnointerventiondf.phi_plus21, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(sim3model1leadlagdf.phi_plus21, [ 0.05, 0.5, 0.95 ]))

exp.(quantile(sim4leadlagnointerventiondf.phi_minus21, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(sim4model1leadlagdf.phi_minus21, [ 0.05, 0.5, 0.95 ]))

exp.(quantile(sim4leadlagnointerventiondf.phi_plus21, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(sim4model1leadlagdf.phi_plus21, [ 0.05, 0.5, 0.95 ]))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# US data / figure 4
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

usfig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 250, 450 ))
    gl = GridLayout(fig[1, 1])

    attaxs = [ Axis(gl[i, 1]) for i ∈ [ 1, 3 ] ]
    plotphi_att!(attaxs[1], usdataleadlagdf, -21:7:21; )
    plotphi_att!(attaxs[2], usdataconfoundersleadlagdf, -21:7:21; )

    for row ∈ [ 1, 3 ] 
        Label(
            gl[row, 0], L"$\varphi$"; 
            fontsize=11.84, rotation=π/2, tellheight=false,
        )
        Label(
            gl[row+1, 1], "Time, days from intervention"; 
            fontsize=11.84, tellwidth=false,
        )
    end

    linkaxes!(attaxs...)

    colgap!(gl, 1, 5)
    for r ∈ [ 1, 3 ] rowgap!(gl, r, 5) end

    labelplots!([ "A", "B" ], gl; cols=0, rows=[ 1, 3 ])
    fig
end
safesave(plotsdir("usfig.pdf"), usfig)

exp.(quantile(usdataleadlagdf.phi_plus21, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(usdataconfoundersleadlagdf.phi_plus21, [ 0.05, 0.5, 0.95 ]))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# UK data / figure 5
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ukdatafig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 400 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])

    attax1 = Axis(ga[3, 1:4])
    attax2 = Axis(gb[3, 1:4])

    let 
        axs = [ 
            Axis(
                ga[1, i]; 
                xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ),
                xticklabelrotation=-π/4,
            ) 
            for i ∈ 1:4 
        ]

        for i ∈ 1:4
            lines!(
                axs[i], 1:257, ukdatalogeffectivereproductionratios.med[:, i];
                color=COLOURVECTOR[1], linewidth=1,
            )
            band!(
                axs[i], 
                1:257, 
                ukdatalogeffectivereproductionratios.lci[:, i],
                ukdatalogeffectivereproductionratios.uci[:, i];
                color=( COLOURVECTOR[1], 0.5 ), 
            )
            noninfindex = findall(x -> -Inf < x < Inf, W_maskcoviddata[:, i])
            scatter!(
                axs[i], collect(1:257)[noninfindex], W_maskcoviddata[noninfindex, i];
                color=:black, markersize=2,
            )

            linkaxes!(axs...)
            formataxis!(
                axs[i]; 
                hidespines=( :r, :t ), hidey=(i != 1), hideyticks=(i != 1), trimspines=true,
            )
            if i != 1 hidespines!(axs[i], :l) end
        end
    end

    plotphi_att!(attax1, ukdataleadlagdf, -21:7:21; )

    let 
        axs = [ 
            Axis(
                gb[1, i]; 
                xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ),
                xticklabelrotation=-π/4,
            ) 
            for i ∈ 1:4 
        ]

        for i ∈ 1:4
            lines!(
                axs[i], 1:257, ukdataconfounderslogeffectivereproductionratios.med[:, i];
                color=COLOURVECTOR[1], linewidth=1,
            )
            band!(
                axs[i], 
                1:257, 
                ukdataconfounderslogeffectivereproductionratios.lci[:, i],
                ukdataconfounderslogeffectivereproductionratios.uci[:, i];
                color=( COLOURVECTOR[1], 0.5 ), 
            )
            noninfindex = findall(x -> -Inf < x < Inf, W_maskcoviddata[:, i])
            scatter!(
                axs[i], collect(1:257)[noninfindex], W_maskcoviddata[noninfindex, i];
                color=:black, markersize=2,
            )

            linkaxes!(axs...)
            formataxis!(
                axs[i]; 
                hidespines=( :r, :t ), hidey=(i != 1), hideyticks=(i != 1), trimspines=true,
            )
            if i != 1 hidespines!(axs[i], :l) end
        end
    end

    plotphi_att!(attax2, ukdataconfounders2leadlagdf, -21:7:21; )


    for gl ∈ [ ga, gb ]
        for (i, c) ∈ enumerate([ "England", "N. Ireland", "Scotland", "Wales"])
            Label(gl[0, i], c; fontsize=10, halign=:left, tellwidth=false)
        end

        Label(
            gl[1, 0], L"log $\mathcal{R}_{e}$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            gl[3, 0], L"$\varphi$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(gl[2, 1:4], "Date, 2020"; fontsize=11.84, tellwidth=false)
        Label(gl[4, 1:4], "Time, days from intervention"; fontsize=11.84, tellwidth=false)
        colgap!(gl, 1, 5) 
        for c ∈ 2:4 colgap!(gl, c, 7) end
        for r ∈ [ 1, 2, 4 ] rowgap!(gl, r, 5) end
        #rowsize!(gl, 5, Auto(1.3))
    end
    
    linkaxes!(attax1, attax2)
    labelplots!([ "A", "B" ], [ ga, gb ])
    fig
end
safesave(plotsdir("ukdatafig.pdf"), ukdatafig)

exp.(quantile(ukdataleadlagdf.phi_minus7, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(ukdataleadlagdf.phi_plus21, [ 0.05, 0.5, 0.95 ]))

exp.(quantile(ukdataconfounders2leadlagdf.phi_minus7, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(ukdataconfounders2leadlagdf.phi_plus21, [ 0.05, 0.5, 0.95 ]))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulated data / supplementary figures 1 and 2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

simulationsR0fig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 600 ))
    gls = [ GridLayout(fig[i, 0:3]) for i ∈ 1:4 ]

    let 
        axs = [ Axis(gls[1][1, i]) for i ∈ 1:3 ]
        lines!(
            axs[1], 1:100, [ beta1a(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[2], 1:100, [ beta1bcounterfactual(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[3], 1:100, [ beta1b(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        linkaxes!(axs...)
        for i ∈ 1:3 
            formataxis!(
                axs[i];
                trimspines=true, hidey=(i > 1), hideyticks=(i > 1), hidespines=( :t, :r ),
            )
            i == 1 && continue
            hidespines!(axs[i], :l)
        end
    end
    let 
        axs = [ Axis(gls[2][j, i]) for i ∈ 1:3, j ∈ 1:2 ]
        lines!(
            axs[1, 1], 1:100, [ beta2a(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[2, 1], 1:100, [ beta2bcounterfactual(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[3, 1], 1:100, [ beta2b(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[1, 2], 1:100, [ beta2a(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[2, 2], 1:100, [ beta2ccounterfactual(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[3, 2], 1:100, [ beta2c(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        linkaxes!(axs...)
        for i ∈ 1:3, j ∈ 1:2
            formataxis!(
                axs[i, j];
                trimspines=true, 
                hidex=(j == 1), hidexticks=(j == 1), 
                hidey=(i > 1), hideyticks=(i > 1), 
                hidespines=( :t, :r ),
            )
            if i > 1 hidespines!(axs[i, j], :l) end
            if j == 1 hidespines!(axs[i, j], :b) end
        end
    end
    let 
        axs = [ Axis(gls[3][1, i]) for i ∈ 1:3 ]
        lines!(
            axs[1], 1:100, [ beta3a(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[2], 1:100, [ beta3bcounterfactual(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[3], 1:100, [ beta3b(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        linkaxes!(axs...)
        for i ∈ 1:3 
            formataxis!(
                axs[i];
                trimspines=true, hidey=(i > 1), hideyticks=(i > 1), hidespines=( :t, :r ),
            )
            i == 1 && continue
            hidespines!(axs[i], :l)
        end
    end
    let 
        axs = [ Axis(gls[4][1, i]) for i ∈ 1:3 ]
        lines!(
            axs[1], 1:100, [ beta4a(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[2], 1:100, [ beta4bcounterfactual(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[3], 1:100, [ beta4b(t) / 0.4 for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        linkaxes!(axs...)
        for i ∈ 1:3 
            formataxis!(
                axs[i];
                trimspines=true, hidey=(i > 1), hideyticks=(i > 1), hidespines=( :t, :r ),
            )
            i == 1 && continue
            hidespines!(axs[i], :l)
        end
    end
    

    for (i, gl) ∈ enumerate(gls)
        i == 2 && continue
        if i != 5 
            Label(
                gl[1, 0], L"$\mathcal{R}_0$"; 
                fontsize=11.84, rotation=π/2, tellheight=false
            )
        else
            Label(
                gl[1, 0], "Proportion\ndiagnosed"; 
                fontsize=11.84, rotation=π/2, tellheight=false
            )
        end
        Label(gl[2, 1:3], "Time, days"; fontsize=11.84, tellwidth=false)
        colgap!(gl, 1, 5)
        rowgap!(gl, 1, 5)
    end
    Label(
        gls[2][1:2, 0], L"$\mathcal{R}_0$"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    Label(gls[2][3, 1:3], "Time, days"; fontsize=11.84, tellwidth=false)
    Label(
        fig[0, 1], "Control"; 
        fontsize=11.84, halign=:left, tellwidth=false
    )
    Label(
        fig[0, 2], "Ineffective intervention"; 
        fontsize=11.84, halign=:left, tellwidth=false
    )
    Label(
        fig[0, 3], "Effective intervention"; 
        fontsize=11.84, halign=:left, tellwidth=false
    )
    rowgap!(fig.layout, 1, 5)
    colgap!(gls[2], 1, 5)
    rowgap!(gls[2], 1, 7)
    rowgap!(gls[2], 2, 5)
    colsize!(fig.layout, 0, Auto(0.2))
    rowsize!(fig.layout, 2, Auto(1.5))

    labelplots!([ "A", "B", "C", "D", ], gls; rows=1)
    fig
end
safesave(plotsdir("simulationsR0fig.pdf"), simulationsR0fig)

simulationsthetafig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 200 ))
    ga = GridLayout(fig[1, 1])
    let 
        axs = [ Axis(ga[1, i]) for i ∈ 1:3 ]
        lines!(
            axs[1], 1:100, [ theta4a(t) for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[2], 1:100, [ theta4b(t) for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        lines!(
            axs[3], 1:100, [ theta4b(t) for t ∈ 1:100 ]; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        linkaxes!(axs...)
        for i ∈ 1:3 
            formataxis!(
                axs[i];
                trimspines=true, 
                hidey=(i > 1), hideyticks=(i > 1), hidespines=( :t, :r ),
                setpoint=(i == 1 ? 0.2 : 0.4),
            )
            i == 1 && continue
            hidespines!(axs[i], :l)
        end
    end

    Label(
        ga[1, 0], "Proportion diagnosed"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    Label(ga[2, 1:3], "Time, days"; fontsize=11.84, tellwidth=false)
    colgap!(ga, 1, 5)
    rowgap!(ga, 1, 5)
    Label(
        ga[0, 1], "Control"; 
        fontsize=11.84, halign=:left, tellwidth=false
    )
    Label(
        ga[0, 2], "Ineffective intervention"; 
        fontsize=11.84, halign=:left, tellwidth=false
    )
    Label(
        ga[0, 3], "Effective intervention"; 
        fontsize=11.84, halign=:left, tellwidth=false
    )
    for r ∈ 1:2 rowgap!(ga, r, 5) end
    colgap!(ga, 1, 5)
 
    fig
end
safesave(plotsdir("simulationsthetafig.pdf"), simulationsthetafig)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# US data trace plots / supplementary figure 3 and 4
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function traceplotlabels(text)
    if text == "logzeta_mean" 
        return "ln μ_ζ"
    elseif text == "zeta_sigma2"
        return "ln σ_ζ"
    elseif length(text) >= 7 && text[1:7] == "logzeta"
        _d = findfirst('.', text)
        return "ln ζ_$(text[10:(_d - 1)])"
    elseif text == "logeta_mean" 
        return "ln μ_η"
        elseif text == "eta_sigma2"
        return "ln σ_η"
    elseif length(text) >= 4 && text[1:3] == "eta"
        _d = findfirst('.', text)
        return "ln η_$(text[6:(_d - 1)])"
    elseif text == "logdelta"
        return "ln τ"
    elseif length(text) >= 17 && text[1:17] == "logsecondarydelta"
        _d = findfirst('.', text)
        return "ln τ_$(text[18:(_d - 1)])"
    elseif text == "log_density"
        return "density"
    elseif length(text) >= 5 && text[1:5] == "phi_m"
        return "ϕ_-$(text[10:end])"
    elseif length(text) >= 5 && text[1:5] == "phi_p"
        return "ϕ_$(text[9:end])"
    elseif text == "phi_0"
        return "ϕ_0"
    else
        return text
    end
end

ustraceplot = with_theme(theme_latexfonts()) do 
    ylbls = traceplotlabels.(names(usdataleadlagdf))
    fig = Figure(; size=( 600, 800 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    gc = GridLayout(fig[1, 3])
    gd = GridLayout(fig[1, 4])
    plotchains!(
        ga, usdataleadlagdf; 
        plotnames_ind=[ [ 72 ]; 3:21 ], ylabels=ylbls[[ [ 72 ]; 3:21 ]], ylabelrotation=π/3,
    )
    plotchains!(
        gb, usdataleadlagdf; 
        plotnames_ind=22:41, ylabels=ylbls[22:41], ylabelrotation=π/3,
    )
    plotchains!(
        gc, usdataleadlagdf; 
        plotnames_ind=42:61, ylabels=ylbls[42:61], ylabelrotation=π/3,
    )
    plotchains!(
        gd, usdataleadlagdf; 
        plotnames_ind=[ 62:71; 73:79 ], ylabels=ylbls[[ 62:71; 73:79 ]], ylabelrotation=π/3,
    )

    fig
end
safesave(plotsdir("ustraceplot.pdf"), ustraceplot)

ustraceplotconfounders = with_theme(theme_latexfonts()) do 
    ylbls = traceplotlabels.(names(usdataconfoundersleadlagdf))
    fig = Figure(; size=( 600, 800 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    gc = GridLayout(fig[1, 3])
    gd = GridLayout(fig[1, 4])
    plotchains!(
        ga, usdataconfoundersleadlagdf; 
        plotnames_ind=[ [ 77 ]; 3:22 ], ylabels=ylbls[[ [ 77 ]; 3:22 ]], ylabelrotation=π/3,
    )
    plotchains!(
        gb, usdataconfoundersleadlagdf; 
        plotnames_ind=23:43, ylabels=ylbls[23:43], ylabelrotation=π/3,
    )
    plotchains!(
        gc, usdataconfoundersleadlagdf; 
        plotnames_ind=44:64, ylabels=ylbls[44:64], ylabelrotation=π/3,
    )
    plotchains!(
        gd, usdataconfoundersleadlagdf; 
        plotnames_ind=[ 65:76; 78:84 ], ylabels=ylbls[[ 65:76; 78:84 ]], ylabelrotation=π/3,
    )

    fig
end
safesave(plotsdir("ustraceplotconfounders.pdf"), ustraceplotconfounders)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# US data trace plots / supplementary figure 5 and 6
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

uktraceplot = with_theme(theme_latexfonts()) do 
    ylbls = traceplotlabels.(names(ukdataleadlagdf))
    fig = Figure(; size=( 600, 800 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    gc = GridLayout(fig[1, 3])
    plotchains!(
        ga, ukdataleadlagdf; 
        plotnames_ind=[ [ 28 ]; 3:11 ], ylabels=ylbls[[ [ 28 ]; 3:11 ]], ylabelrotation=π/3,
    )
    plotchains!(
        gb, ukdataleadlagdf; 
        plotnames_ind=12:23, ylabels=ylbls[12:23], ylabelrotation=π/3,
    )
    plotchains!(
        gc, ukdataleadlagdf; 
        plotnames_ind=[ 24:27; 29:35 ], ylabels=ylbls[[ 24:27; 29:35 ]], ylabelrotation=π/3,
    )

    fig
end
safesave(plotsdir("uktraceplot.pdf"), uktraceplot)

uktraceplotconfounders = with_theme(theme_latexfonts()) do 
    ylbls = traceplotlabels.(names(ukdataconfounders2leadlagdf))
    colnames = traceplotlabels.(names(ukdataconfounders2leadlagdf))
    fig = Figure(; size=( 600, 800 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    gc = GridLayout(fig[1, 3])
    plotchains!(
        ga, ukdataconfounders2leadlagdf; 
        plotnames_ind=[ [ 31 ]; 3:12 ], ylabels=ylbls[[ [ 31 ]; 3:12 ]], ylabelrotation=π/3,
    )
    plotchains!(
        gb, ukdataconfounders2leadlagdf; 
        plotnames_ind=13:25, ylabels=ylbls[13:25], ylabelrotation=π/3,
    )
    plotchains!(
        gc, ukdataconfounders2leadlagdf; 
        plotnames_ind=[ 26:30; 32:38 ], ylabels=ylbls[[ 26:30; 32:38 ]], ylabelrotation=π/3,
    )

    fig
end
safesave(plotsdir("uktraceplotconfounders.pdf"), uktraceplotconfounders)
