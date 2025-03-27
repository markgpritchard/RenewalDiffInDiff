
using DrWatson
@quickactivate :RenewalDiffInDiff

using CairoMakie, StatsBase
include("loadanalysisforplots.jl")
include(srcdir("plottingfunctions.jl"))


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
# Simulations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim1_0_fitzfig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 500 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])

    attax1 = Axis(ga[3, 1:2])
    attax2 = Axis(gb[3, 1:2])
    casesax1 = Axis(ga[5, 1:2])
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

        plottau_att!(attax1, sim1leadlagnointerventiondf, -21:7:21; )

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

        plottau_att!(attax2, sim1model1lagleaddf, -21:7:21; )

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

        formataxis!(casesax2; hidespines=( :r, :t ), trimspines=true,)
    end

    for gl ∈ [ ga, gb ]
        Label(
            gl[1, 0], L"log $\mathcal{R}_{e}$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            gl[3, 0], L"$\varphi_{\mathrm{ATT}}$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            gl[5, 0], L"cumulative difference \\ in diagnoses, per $10\,000$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(gl[2, 1:2], "Time"; fontsize=11.84, tellwidth=false)
        Label(gl[4, 1:2], "Time from intervention"; fontsize=11.84, tellwidth=false)
        Label(gl[6, 1:2], "Time"; fontsize=11.84, tellwidth=false)
        colgap!(gl, 1, 5) 
        for r ∈ [ 1, 3, 5 ] rowgap!(gl, r, 5) end
    end
    
    linkaxes!(attax1, attax2)
    linkaxes!(casesax1, casesax2)

    labelplots!([ "A", "B" ], [ ga, gb ]; rows=1)
    fig
end

safesave(plotsdir("sim1_0_fitzfig.pdf"), sim1_0_fitzfig)

exp.(quantile(sim1leadlagnointerventiondf.tau_plus21, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(sim1model1lagleaddf.tau_plus21, [ 0.05, 0.5, 0.95 ]))

simfig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 250, 450 ))
    ga = GridLayout(fig[1, 1])
    axs = [ Axis(ga[(2 * i - 1), j]; xticks=[ -21, 0, 21 ]) for i ∈ 1:3, j ∈ 1:2 ]

    plottau_att!(axs[1, 1], sim2leadlagconfoundernointerventiondf, -21:7:21; )
    plottau_att!(axs[1, 2], sim2model1leadlagdf, -21:7:21; )
    plottau_att!(axs[2, 1], sim3leadlagnointerventiondf, -21:7:21; )
    plottau_att!(axs[2, 2], sim3model1leadlagdf, -21:7:21; )
    plottau_att!(axs[3, 1], sim4leadlagnointerventiondf, -21:7:21; )
    plottau_att!(axs[3, 2], sim4model1leadlagdf, -21:7:21; )

    for row ∈ [ 1, 3, 5 ] 
        Label(
            ga[row, 0], L"$\varphi_{\mathrm{ATT}}$"; 
            fontsize=11.84, rotation=π/2, tellheight=false,
        )
        Label(
            ga[row+1, 1:2], "Time, relative to intervention"; 
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

exp.(quantile(sim2leadlagconfoundernointerventiondf.tau_plus21, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(sim2model1leadlagdf.tau_plus21, [ 0.05, 0.5, 0.95 ]))

exp.(quantile(sim3leadlagnointerventiondf.tau_plus21, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(sim3model1leadlagdf.tau_plus21, [ 0.05, 0.5, 0.95 ]))

exp.(quantile(sim4leadlagnointerventiondf.tau_minus21, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(sim4model1leadlagdf.tau_minus21, [ 0.05, 0.5, 0.95 ]))

exp.(quantile(sim4leadlagnointerventiondf.tau_plus21, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(sim4model1leadlagdf.tau_plus21, [ 0.05, 0.5, 0.95 ]))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# UK data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ukdatafig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 500 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    ga1 = GridLayout(ga[5, 0:4])
    gb1 = GridLayout(gb[5, 0:4])

    attax1 = Axis(ga[3, 1:4])
    attax2 = Axis(gb[3, 1:4])

    casesaxs1 = [ 
        Axis(
            ga1[1, i]; 
            xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ),
            xticklabelrotation=-π/4,
        ) 
        for i ∈ 1:3
    ]
    casesaxs2 = [ 
        Axis(
            gb1[1, i]; 
            xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ),
            xticklabelrotation=-π/4,
        ) 
        for i ∈ 1:3
    ]

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

        plottau_att!(attax1, ukdataleadlagdf, -21:7:21; )

        for (i, cases) ∈ enumerate(
                [ 
                ukdatacasesdiff_167, ukdatacasesdiff_192, ukdatacasesdiff_174 
                ]
            )
            lines!(
                casesaxs1[i], 1:257, cases.med[:, i] .* 10_000 ./ POPULATION2020[i];
                color=COLOURVECTOR[1], linewidth=1,
            )
            band!(
                casesaxs1[i], 
                1:257, 
                cases.lci[:, i] .* 10_000 ./ POPULATION2020[i],
                cases.uci[:, i] .* 10_000 ./ POPULATION2020[i];
                color=( COLOURVECTOR[1], 0.5 ), 
            )
            hlines!(
                casesaxs1[i], 0; 
                color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
            )

            formataxis!(
                casesaxs1[i]; 
                hidespines=( :r, :t ), hidey=(i != 1), hideyticks=(i != 1), trimspines=true,
            )
            if i != 1 hidespines!(casesaxs1[i], :l) end
        end
    end

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

        plottau_att!(attax2, ukdataconfounders2leadlagdf, -21:7:21; )

        for (i, cases) ∈ enumerate(
            [ 
                ukdataconfounderscasesdiff_167, 
                ukdataconfounderscasesdiff_192, 
                ukdataconfounderscasesdiff_174 
            ]
        )
            lines!(
                casesaxs2[i], 1:257, cases.med[:, i] .* 10_000 ./ POPULATION2020[i];
                color=COLOURVECTOR[1], linewidth=1,
            )
            band!(
                casesaxs2[i], 
                1:257, 
                cases.lci[:, i] .* 10_000 ./ POPULATION2020[i],
                cases.uci[:, i] .* 10_000 ./ POPULATION2020[i];
                color=( COLOURVECTOR[1], 0.5 ), 
            )
            hlines!(
                casesaxs2[i], 0; 
                color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
            )

            formataxis!(
                casesaxs2[i]; 
                hidespines=( :r, :t ), hidey=(i != 1), hideyticks=(i != 1), trimspines=true,
            )
            if i != 1 hidespines!(casesaxs2[i], :l) end
        end
    end

    for gl ∈ [ ga, gb ]
        for (i, c) ∈ enumerate([ "England", "N. Ireland", "Scotland", "Wales"])
            Label(gl[0, i], c; fontsize=10, halign=:left, tellwidth=false)
        end

        Label(
            gl[1, 0], L"log $\mathcal{R}_{e}$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            gl[3, 0], L"$\varphi_{\mathrm{ATT}}$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(gl[2, 1:4], "Date, 2020"; fontsize=11.84, tellwidth=false)
        Label(gl[4, 1:4], "Time from intervention"; fontsize=11.84, tellwidth=false)
        colgap!(gl, 1, 5) 
        for c ∈ 2:4 colgap!(gl, c, 7) end
        for r ∈ [ 1, 2, 4 ] rowgap!(gl, r, 5) end
        rowsize!(gl, 5, Auto(1.3))
    end

    for gl ∈ [ ga1, gb1 ]
        for (i, c) ∈ enumerate([ "England", "N. Ireland", "Scotland" ])
            Label(gl[0, i], c; fontsize=10, halign=:left, tellwidth=false)
        end
        Label(
            gl[1, 0], L"cumulative difference \\ in diagnoses, per $10\,000$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(gl[2, 1:3], "Date, 2020"; fontsize=11.84, tellwidth=false)
        colgap!(gl, 1, 5) 
        for c ∈ 2:3 colgap!(gl, c, 7) end
        for r ∈ [ 1, 2 ] rowgap!(gl, r, 5) end 
    end
    
    linkaxes!(attax1, attax2)
    linkaxes!(casesaxs1..., casesaxs2...)
    labelplots!([ "A", "B" ], [ ga, gb ])
    fig
end

safesave(plotsdir("ukdatafig.pdf"), ukdatafig)

exp.(quantile(ukdataleadlagdf.tau_minus7, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(ukdataleadlagdf.tau_plus21, [ 0.05, 0.5, 0.95 ]))

exp.(quantile(ukdataconfounders2leadlagdf.tau_plus21, [ 0.05, 0.5, 0.95 ]))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# US data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

usfig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 400 ))
    gl = GridLayout(fig[1, 1])

    attaxs = [ Axis(gl[1, j]) for j ∈ [ 1, 3 ] ]
    incidentaxs = [ 
        Axis(
            gl[3, j]; 
            xticks=( 1:13, usstatedataleadlag.plotstates),
            xticklabelrotation=-π/2,
        ) 
        for j ∈ [ 1, 3 ] 
    ]

    plottau_att!(attaxs[1], usdataleadlagdf, -21:7:21; )
    plottau_att!(attaxs[2], usdataconfoundersleadlagdf, -21:7:21; )

    for (i, d) ∈ enumerate([ usstatedataleadlag, usstatedataleadlag_confounder ])
        scatter!(
            incidentaxs[i], 1:13, d.plotmedians; 
            color=COLOURVECTOR[1], marker=:x, markersize=5,
        )
        rangebars!(
            incidentaxs[i], 1:13, d.plotlcis, d.plotucis; 
            color=COLOURVECTOR[1], linewidth=1,
        )
        hlines!(
            incidentaxs[i], 0; 
            color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
        )
        formataxis!(incidentaxs[i]; hidespines=( :r, :t ), trimspines=true)
    end

    for col ∈ [ 1, 3 ] 
        Label(
            gl[1, col-1], L"$\varphi_{\mathrm{ATT}}$"; 
            fontsize=11.84, rotation=π/2, tellheight=false,
        )
        Label(
            gl[2, col], "Time, relative to intervention"; 
            fontsize=11.84, tellwidth=false,
        )
        Label(
            gl[3, col-1], L"cumulative difference \\ in diagnosis, per $10\,000$"; 
            fontsize=11.84, rotation=π/2, tellheight=false,
        )
        Label(
            gl[4, col], "State"; 
            fontsize=11.84, tellwidth=false,
        )
    end

    linkaxes!(attaxs...)
    linkaxes!(incidentaxs...)

    for cr ∈ [ 1, 3 ]
        colgap!(gl, cr, 5)
        rowgap!(gl, cr, 5)
    end

    labelplots!([ "A", "B" ], gl; cols=[ 0, 2 ], rows=1)
    fig
end

safesave(plotsdir("usfig.pdf"), usfig)

exp.(quantile(usdataleadlagdf.tau_plus21, [ 0.05, 0.5, 0.95 ]))
exp.(quantile(usdataconfoundersleadlagdf.tau_plus21, [ 0.05, 0.5, 0.95 ]))

