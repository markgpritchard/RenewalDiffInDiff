
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
    fig = Figure(; size=( 500, 350 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[2, 1])

    let 
        axs = [ Axis(ga[1, i]) for i ∈ 1:2 ]
        attax = Axis(ga[1, 4])
        casesax = Axis(ga[1, 6])

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

            formataxis!(
                axs[i]; 
                hidespines=( :r, :t ), hidey=(i != 1), hideyticks=(i != 1), trimspines=true,
            )
            if i == 2 hidespines!(axs[i], :l) end
        end

        plottau_att!(attax, sim1leadlagnointerventiondf, -21:7:21; )

        lines!(
            casesax, 1:100, sim1nointerventioncasesdiff.med[:, 2] ./ 10_000;
            color=COLOURVECTOR[1], linewidth=1,
        )
        band!(
            casesax, 
            1:100, 
            sim1nointerventioncasesdiff.lci[:, 2] ./ 10_000,
            sim1nointerventioncasesdiff.uci[:, 2] ./ 10_000;
            color=( COLOURVECTOR[1], 0.5 ), 
        )
        hlines!(
            casesax, 0; 
            color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
        )

        formataxis!(casesax; hidespines=( :r, :t ), trimspines=true,)
    end

    let 
        axs = [ Axis(gb[1, i]) for i ∈ 1:2 ]
        attax = Axis(gb[1, 4])
        casesax = Axis(gb[1, 6])

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

        plottau_att!(attax, sim1model1lagleaddf, -21:7:21; )

        lines!(
            casesax, 1:100, sim1casesdiff.med[:, 2] ./ 10_000;
            color=COLOURVECTOR[1], linewidth=1,
        )
        band!(
            casesax, 
            1:100, 
            sim1casesdiff.lci[:, 2] ./ 10_000,
            sim1casesdiff.uci[:, 2] ./ 10_000;
            color=( COLOURVECTOR[1], 0.5 ), 
        )
        hlines!(
            casesax, 0; 
            color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
        )

        formataxis!(casesax; hidespines=( :r, :t ), trimspines=true,)
    end

    for gl ∈ [ ga, gb ]
        Label(
            gl[1, 0], "log effective\nreproduction ratio"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            gl[1, 3], L"$\tau_{\mathrm{ATT}}$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            gl[1, 5], L"cumulative difference \\ in incidence, per $10\,000$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(gl[2, 1:2], "Time"; fontsize=11.84, tellwidth=false)
        Label(gl[2, 4], "Time, relative\nto intervention"; fontsize=11.84, tellwidth=false)
        Label(gl[2, 6], "Time"; fontsize=11.84, tellwidth=false)
        for c ∈ [ 1, 4, 6 ] colgap!(gl, c, 5) end
        rowgap!(gl, 1, 5)
    end

    labelplots!([ "A", "B" ], [ ga, gb ]; rows=1)
    fig
end

safesave(plotsdir("sim1_0_fitzfig.pdf"), sim1_0_fitzfig)

simfig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 250, 800 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[2, 1])
    gc = GridLayout(fig[3, 1])
    gd = GridLayout(fig[4, 1])
    ge = GridLayout(fig[5, 1])

    plottau_att!(ga, sim1leadlagnointerventiondf, sim1model1lagleaddf, -21:7:21; )
    plottau_att!(gb, sim2leadlagnointerventiondf, sim2modelleadlagdf, -21:7:21; )
    plottau_att!(gc, sim2leadlagconfoundernointerventiondf, sim2model1leadlagdf, -21:7:21; )
    plottau_att!(gd, sim3leadlagnointerventiondf, sim3model1leadlagdf, -21:7:21; )
    plottau_att!(ge, sim4leadlagnointerventiondf, sim4model1leadlagdf, -21:7:21; )

    labelplots!([ "A", "B", "C", "D", "E" ], [ ga, gb, gc, gd, ge ]; rows=1)
    fig
end

safesave(plotsdir("simfig.pdf"), simfig)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Covid data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ukdatafig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 800 ))
    ga = GridLayout(fig[1, 1])

    let 
        axs = [ 
            Axis(
                ga[j, i]; 
                xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] )
            ) 
            for i ∈ 1:2, j ∈ [ 1, 3 ] 
        ]
        plotrenewalequationsamples_w!(
            axs, 
            maskcovidcases, W_maskcoviddata, ukdataleadlagfittedoutput, fitws(
                maskcovidcases, 
                POPULATION2020, 
                ukdataleadlagfittedoutput
            );
            markersize=2,
        )

        vlines!(axs[1], 167; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        vlines!(axs[2], 192; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        vlines!(axs[3], 174; color=:red, linestyle=( :dot, :dense ), linewidth=1,)

        Label(
            ga[0, 1], "England"; 
            fontsize=10, halign=:left, tellwidth=false
        )
        formataxis!(
            axs[1]; 
            hidespines=( :r, :t, :b ), hidex=true, hidexticks=true, trimspines=true
        )
        Label(
            ga[0, 2], "Northern Ireland"; 
            fontsize=10, halign=:left, tellwidth=false
        )
        formataxis!(
            axs[2]; 
            hidespines=( :r, :t, :b, :l ), 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, trimspines=true
        )
        Label(
            ga[2, 1], "Scotland"; 
            fontsize=10, halign=:left, tellwidth=false
        )
        formataxis!(axs[3]; hidespines=( :r, :t, ), trimspines=true)
        Label(
            ga[2, 2], "Wales"; 
            fontsize=10, halign=:left, tellwidth=false
        )
        formataxis!(
            axs[4]; 
            hidespines=( :r, :t, :l ), hidey=true, hideyticks=true, trimspines=true
        )
        Label(
            ga[1:3, 0], L"$\ln\mathcal{R}_e$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(ga[4, 1:2], "Date, 2020"; fontsize=11.84, tellwidth=false)
    end

    let 
        axs = [ Axis(ga[j, 4]) for j ∈ [ 1:3, 5:7 ]]
        plottau_att!(axs[1], ukdataleadlagdf, -21:7:21; )
        plottau_att!(axs[2], ukdataconfounders2leadlagdf, -21:7:21; )

        for row in [ 1, 5 ]
            Label(
                ga[row:(row+3), 3], L"$\tau_{\mathrm{ATT}}$"; 
                fontsize=11.84, rotation=π/2, tellheight=false,
            )
            Label(
                ga[row+3, 4], "Time, relative to intervention"; 
                fontsize=11.84, tellwidth=false,
            )
        end
        linkaxes!(axs...)
    end

    ga1 = GridLayout(ga[0:4, 5])
    let 
        axs = [ 
            Axis(
                ga1[j, 1]; 
                xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] )
            ) 
            for j ∈ [ 1, 3, 5 ]
        ]
        for (i, d) ∈ enumerate([ 
                ukdataleadlagfittedoutput_167, 
                ukdataleadlagfittedoutput_192, 
                ukdataleadlagfittedoutput_174
            ])
            plotrenewalequationsamples_causaleffect!(
                axs[i], maskcovidcases, nothing, POPULATION2020, d, i;
                cumulativedifference=true,
                fittedparameter=:y_matrix_det_vec,
                counterfactualfittedparameter=:y_matrix_det_vec_counterfactual,
                xticklabelrotation=-π/4,
                xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
                xtitle="Date, 2020",
                ytitle=L"Cumulative \\ difference$$",
            )
        end
        vlines!(axs[1], 167; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        vlines!(axs[2], 192; color=:red, linestyle=( :dot, :dense ), linewidth=1,)
        vlines!(axs[3], 174; color=:red, linestyle=( :dot, :dense ), linewidth=1,)

        Label(
            ga1[0, 1], "England"; 
            fontsize=10, halign=:left, tellwidth=false
        )
        formataxis!(
            axs[1]; 
            hidespines=( :r, :t, :b ), hidex=true, hidexticks=true, trimspines=true
        )
        Label(
            ga1[2, 1], "Northern Ireland"; 
            fontsize=10, halign=:left, tellwidth=false
        )
        formataxis!(
            axs[2]; 
            hidespines=( :r, :t, :b ), hidex=true, hidexticks=true, trimspines=true
        )
        Label(
            ga1[4, 1], "Scotland"; 
            fontsize=10, halign=:left, tellwidth=false
        )
        formataxis!(axs[3]; hidespines=( :r, :t, ), trimspines=true)

        Label(
            ga1[1:5, 0], "Change in cumulative incidence"; 
            fontsize=11.84, rotation=π/2, tellheight=false,
        )
        Label(
            ga1[6, 1], "Date, 2020"; 
            fontsize=11.84, tellwidth=false,
        )

        colgap!(ga1, 1, 5)
        for r ∈ [ 1, 3, 5, 6 ] rowgap!(ga1, r, 5) end
    end

    for c ∈ [ 1, 4 ] colgap!(ga, c, 5) end
    for r ∈ [ 1, 3, 4 ] rowgap!(ga, r, 5) end
        
    labelplots!([ "A"  ], ga; rows=[ 0 ])
    fig
end

safesave(plotsdir("ukdatafig.pdf"), ukdatafig)

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
        ytitle=L"Without \\ intervetion$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, allcovidcases, selectpops, datafit1, 4;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"With \\ intervetion$$",
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
        ytitle=L"Without \\ intervetion$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, allcovidcases, selectpops, datafit1, 4;
        locationinds=6:8,
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"With \\ intervetion$$",
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
        ytitle=L"Without \\ intervetion$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, pil1covidcases, selectpops, datafit2, 4;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"With \\ intervetion$$",
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
        ytitle=L"Without \\ intervetion$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, pil1covidcases, selectpops, datafit2, 4;
        locationinds=6:8,
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"With \\ intervetion$$",
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
        ytitle=L"Without \\ intervetion$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, allcovidcases, selectpops, datafit3, 4;
        locationinds=[ 1:5; [ 9 ] ],
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 93, 215 ], [ "June", "Sept.", "Jan." ] ), 
        ytitle=L"With \\ intervetion$$",
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
        ytitle=L"Without \\ intervetion$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases[1:191, :], POPULATION2020, maskingdatafit1, 4;
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 92, 183 ], [ "Jan.", "April", "July" ] ), 
        ytitle=L"With \\ intervetion$$",
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
        ytitle=L"Without \\ intervetion$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases, POPULATION2020, maskingdatafit2, 4;
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"With \\ intervetion$$",
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
        ytitle=L"Without \\ intervetion$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases, POPULATION2020, maskingdatafit3, 4;
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"With \\ intervetion$$",
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
        ytitle=L"Without \\ intervetion$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases, POPULATION2020, maskingdatafit4, 4;
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"With \\ intervetion$$",
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
        ytitle=L"Without \\ intervetion$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, maskcovidcases, POPULATION2020, maskingdatafit5, 4;
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 92, 183, 245 ], [ "Jan.", "April", "July", "Sept." ] ), 
        ytitle=L"With \\ intervetion$$",
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



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with no effective intervention
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

simineffectiveplot = with_theme(theme_latexfonts()) do 
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
            sim1fit0, 
            fitws(
                simulation1dataset["cases"], 
                simulation1dataset["Ns"], 
                sim1fit0
            ), 
            1;
            markersize=2,
            hidex=true,
            ytitle=nothing,
            yticks=[ -1, 0, 1, 2 ], 
        )
        setvalue!(axs1[1], -1)
        axs2 = plotrenewalequationsamples_r0!(
            gb, simulation1dataset["cases"], sim1fit0, 2;
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
            sim1fit0, 
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
            sim1fit0,
            4;
            markersize=2, fittedparameter=:y_matrix_det_vec,
            ytitle=nothing,
        )
        axs5 = plotrenewalequationsamples_causaleffect!(
            gb, 
            simulation1dataset["cases"], 
            simulation1dataset["cases_counterfactual"], 
            simulation1dataset["Ns"], 
            sim1fit0, 
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

safesave(plotsdir("simineffectiveplot.pdf"), simineffectiveplot)



###########################################################################################

## Analysis 1 
# Effect of mask recommendations. No other considerations of confounding 

datamodeluschain1 = loadanalysisdictsasdf("datamodelus1", 12, maxrounds, 101)

datamodelusfit1 = samplerenewalequation_2sets(
    COVIDSERIALINTERVAL, datamodeluschain1, maskday; 
    initialvalues=incidence[1:70, :], 
    Ns=populations,
    #psi=0.4, timeknots=collect(1:303/10:304),
    timeknots=[ collect(1.0:28:113); [ 123 ] ],
)

datamodelus1kv = keyvalues(datamodeluschain1, datamodelusfit1)
print(datamodelus1kv)

datamodelus1plot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 5000, 450 ))
    ga = GridLayout(fig[1, 1])
    axs1 = plotrenewalequationsamples_w!(
        ga, 
        incidence, 
        W_uscoviddata, datamodelusfit1, 
        fitws(
            incidence, 
            populations, 
            datamodelusfit1
        ), 
        1;
        markersize=2,
        #xticks=( [ 1, 92, 183 ], [ "Jan.", "April", "July" ] ), 
        hidex=true, ytitle=L"$\ln\mathcal{R}_e$",
    )
    axs2 = plotrenewalequationsamples_r0!(
        ga, incidence, datamodelusfit1, 2;
        plotcounterfactuals=true, 
        #xticks=( [ 1, 92, 183 ], [ "Jan.", "April", "July" ] ), 
        ytitle=L"$\mathcal{R}_0$",
    )
    axs3 = plotrenewalequationsamples_cases!(
        ga, incidence, populations, datamodelusfit1, 3;
        markersize=2, fittedparameter=:y_matrix_det_vec_counterfactual,
        fittedcolour=( COLOURVECTOR[2], 0.75 ), 
        xticks=( [ 1, 92, 183 ], [ "Jan.", "April", "July" ] ), 
        ytitle=L"Without \\ intervetion$$",
    )
    axs4 = plotrenewalequationsamples_cases!(
        ga, incidence, populations, datamodelusfit1, 4;
        markersize=2, fittedparameter=:y_matrix_det_vec,
        xticks=( [ 1, 92, 183 ], [ "Jan.", "April", "July" ] ), 
        ytitle=L"With \\ intervetion$$",
    )
    axs5 = plotrenewalequationsamples_causaleffect!(
        ga, incidence, nothing, populations, datamodelusfit1, 5;
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
#=
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
=#
    fig
end

safesave(plotsdir("subsetmaskdatafit1plot.pdf"), subsetmaskdatafit1plot)



#
#
#


fn = RenewalDiffInDiff._findanalysisfilename("datamodelus1", 12, 12, 101 + 4)
#isnothing(fn) && continue
chain = load(fn)["chain"]
datamodeluschain1 = DataFrame(chain)
_tdf.chain = [ i for _ ∈ axes(_tdf, 1) ]
df = vcat(df, _tdf) 



#
#
#
#



sim1model0laglead

fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, 1:7, exp.(median.([ getproperty(sim1leadlagnointerventiondf, "tau$x") for x ∈ 1:7 ])))
rangebars!(ax, 1:7, exp.(quantile.([ getproperty(sim1leadlagnointerventiondf, "tau$x") for x ∈ 1:7 ], 0.05)), exp.(quantile.([ getproperty(sim1leadlagnointerventiondf, "tau$x") for x ∈ 1:7 ], 0.95)))


fig 




tq = tauquantiles(sim2model1leadlagdf, -21:7:21)

fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, 1:7, exp.(tq.med))
rangebars!(ax, 1:7, exp.(tq.lci), exp.(tq.uci))


fig 




using CubicSplines



asdf2 = logeffectivereproductionratios(
    ukdataleadlagdf,
    4,
    maskcovidcases,
    populations,
    [ 1.0; collect(56.0:28:224); 257 ],
    facialcoveringsrequired,
    lagleadinterventionsmatrix(facialcoveringsrequired_IM, -21:7:21);
)

asdf4 = quantilelogeffectivereproductionratios(
    ukdataleadlagdf,
    4,
    maskcovidcases,
    populations,
    [ 1.0; collect(56.0:28:224); 257 ],
    facialcoveringsrequired,
    lagleadinterventionsmatrix(facialcoveringsrequired_IM, -21:7:21);
)

lines(asdf.med[:, 1])


asdf = basicreproductionratios(
    ukdataleadlagdf,
    4,
    [ 1.0; collect(56.0:28:224); 257 ],
    facialcoveringsrequired,
    lagleadinterventionsmatrix(facialcoveringsrequired_IM, -21:7:21);
)


jkl = similar(asdf)



predictcases(x -> COVIDSERIALINTERVAL[x], asdf, 1, maskcovidcases)


asdf5 = predictcases(
    x -> x > 31 ? 0 : COVIDSERIALINTERVAL[x], 
    ukdataleadlagdf, 
    asdf, 
    maskcovidcases[1:167, :], 
    populations
)

asdf6 = predictcases(
    x -> x > 31 ? 0 : COVIDSERIALINTERVAL[x], 
    ukdataleadlagdf, 
    asdf, 
    maskcovidcases[1:167, :], 
    populations
)