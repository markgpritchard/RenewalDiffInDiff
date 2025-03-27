
using DrWatson
@quickactivate :RenewalDiffInDiff

using StatsBase
include("analysis.jl")  # will run analyses if they are not already available 

const maxrounds = 12  # the greatest number of rounds that the analysis may have run 

sim1leadlagnointerventiondf = loadanalysisdictsasdf("sim1model0laglead", 8, maxrounds, 105)
insertcumulativeeffects!(sim1leadlagnointerventiondf, -21:7:21)
sim1nointerventionlogeffectivereproductionratios = quantilelogeffectivereproductionratios(
    sim1leadlagnointerventiondf,
    2,
    simulation1dataset["cases_counterfactual"],
    simulation1dataset["Ns"],
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation1dataset["interventions"],
    lagleadinterventionsmatrix(simulation1dataset["interventions"], -21:7:21);
    logdelta=0, 
    logsecondarydelta1=0,
    logsecondarydelta2=0, 
    logsecondarydelta3=0, 
    logsecondarydelta4=0, 
    logsecondarydelta5=0, 
    logsecondarydelta6=0, 
)
sim1nointerventioncasesdiff = let 
    logr0a = logbasicreproductionratios(
        sim1leadlagnointerventiondf, 
        2, 
        [ [ 1 ]; collect(11:89/4:100) ], 
        simulation1dataset["interventions"], 
        lagleadinterventionsmatrix(simulation1dataset["interventions"], -21:7:21);
        logdelta=0, 
        logsecondarydelta1=0,
        logsecondarydelta2=0, 
        logsecondarydelta3=0, 
        logsecondarydelta4=0, 
        logsecondarydelta5=0, 
        logsecondarydelta6=0,
    )
    logr0b = logbasicreproductionratios(
        sim1leadlagnointerventiondf, 
        2, 
        [ [ 1 ]; collect(11:89/4:100) ], 
        simulation1dataset["interventions"], 
        lagleadinterventionsmatrix(simulation1dataset["interventions"], -21:7:21);
    )
    quantilepredictcumulativedifferenceincases(
        fseir, 
        sim1leadlagnointerventiondf, 
        logr0a, 
        logr0b,
        simulation1dataset["cases_counterfactual"][1:29, :], 
        simulation1dataset["Ns"]
    )
end
sim1model1lagleaddf = loadanalysisdictsasdf("sim1model1laglead", 8, maxrounds, 115)
insertcumulativeeffects!(sim1model1lagleaddf, -21:7:21)
sim1logeffectivereproductionratios = quantilelogeffectivereproductionratios(
    sim1model1lagleaddf,
    2,
    simulation1dataset["cases"],
    simulation1dataset["Ns"],
    [ [ 1 ]; collect(11:89/4:100) ],
    simulation1dataset["interventions"],
    lagleadinterventionsmatrix(simulation1dataset["interventions"], -21:7:21)
)
sim1casesdiff = let 
    logr0a = logbasicreproductionratios(
        sim1model1lagleaddf, 
        2, 
        [ [ 1 ]; collect(11:89/4:100) ], 
        simulation1dataset["interventions"], 
        lagleadinterventionsmatrix(simulation1dataset["interventions"], -21:7:21);
        logdelta=0, 
        logsecondarydelta1=0,
        logsecondarydelta2=0, 
        logsecondarydelta3=0, 
        logsecondarydelta4=0, 
        logsecondarydelta5=0, 
        logsecondarydelta6=0,
    )
    logr0b = logbasicreproductionratios(
        sim1model1lagleaddf, 
        2, 
        [ [ 1 ]; collect(11:89/4:100) ], 
        simulation1dataset["interventions"], 
        lagleadinterventionsmatrix(simulation1dataset["interventions"], -21:7:21);
    )
    quantilepredictcumulativedifferenceincases(
        fseir, 
        sim1model1lagleaddf, 
        logr0a, 
        logr0b,
        simulation1dataset["cases"][1:29, :], 
        simulation1dataset["Ns"]
    )
end
sim2leadlagnointerventiondf = loadanalysisdictsasdf("sim2model1_0laglead", 8, maxrounds, 255)
insertcumulativeeffects!(sim2leadlagnointerventiondf, -21:7:21)
sim2confoundernointerventiondf = loadanalysisdictsasdf("sim2model2_0", 8, maxrounds, 260)
sim2leadlagconfoundernointerventiondf = loadanalysisdictsasdf(
    "sim2model2_0laglead", 8, maxrounds, 265
)
insertcumulativeeffects!(sim2leadlagconfoundernointerventiondf, -21:7:21; deltaindex=2:7)
sim2modelleadlagdf = loadanalysisdictsasdf("sim2model0laglead", 8, maxrounds, 205)
insertcumulativeeffects!(sim2modelleadlagdf, -21:7:21)
sim2model1leadlagdf = loadanalysisdictsasdf("sim2model1laglead", 8, maxrounds, 215)
insertcumulativeeffects!(sim2model1leadlagdf, -21:7:21; deltaindex=2:7)
sim3leadlagnointerventiondf = loadanalysisdictsasdf("sim3model0leadlag", 8, maxrounds, 355)
insertcumulativeeffects!(sim3leadlagnointerventiondf, -21:7:21)
sim3model1leadlagdf = loadanalysisdictsasdf("sim3model2", 8, maxrounds, 315)
insertcumulativeeffects!(sim3model1leadlagdf, -21:7:21)
sim4leadlagnointerventiondf = loadanalysisdictsasdf("sim4model0leadlag", 8, maxrounds, 405)
insertcumulativeeffects!(sim4leadlagnointerventiondf, -21:7:21)
sim4model1leadlagdf = loadanalysisdictsasdf("sim4model2", 8, maxrounds, 415)
insertcumulativeeffects!(sim4model1leadlagdf, -21:7:21)
ukdataleadlagdf = loadanalysisdictsasdf("maskingdatamodel2leadlag", 8, maxrounds, 1125)
insertcumulativeeffects!(ukdataleadlagdf, -21:7:21)
ukdatalogeffectivereproductionratios = quantilelogeffectivereproductionratios(
    ukdataleadlagdf,
    4,
    maskcovidcases,
    POPULATION2020,
    [ 1.0; collect(56.0:28:224); 257 ],
    facialcoveringsrequired,
    lagleadinterventionsmatrix(facialcoveringsrequired_IM, -21:7:21)
)
let 
    logr0a = logbasicreproductionratios(
        ukdataleadlagdf, 
        4, 
        [ 1.0; collect(56.0:28:224); 257 ], 
        facialcoveringsrequired, 
        lagleadinterventionsmatrix(facialcoveringsrequired_IM, -21:7:21);
        logdelta=0, 
        logsecondarydelta1=0,
        logsecondarydelta2=0, 
        logsecondarydelta3=0, 
        logsecondarydelta4=0, 
        logsecondarydelta5=0, 
        logsecondarydelta6=0,
    )
    logr0b = logbasicreproductionratios(
        ukdataleadlagdf, 
        4, 
        [ 1.0; collect(56.0:28:224); 257 ], 
        facialcoveringsrequired, 
        lagleadinterventionsmatrix(facialcoveringsrequired_IM, -21:7:21);
    )
    global ukdatacasesdiff_167 = quantilepredictcumulativedifferenceincases(
        COVIDSERIALINTERVAL, 
        ukdataleadlagdf, 
        logr0a, 
        logr0b,
        maskcovidcases[1:146, :], 
        POPULATION2020
    )
    global ukdatacasesdiff_192 = quantilepredictcumulativedifferenceincases(
        COVIDSERIALINTERVAL, 
        ukdataleadlagdf, 
        logr0a, 
        logr0b,
        maskcovidcases[1:171, :], 
        POPULATION2020
    )
    global ukdatacasesdiff_174 = quantilepredictcumulativedifferenceincases(
        COVIDSERIALINTERVAL, 
        ukdataleadlagdf, 
        logr0a, 
        logr0b,
        maskcovidcases[1:153, :], 
        POPULATION2020
    )
end
ukdataconfounders1leadlagdf = loadanalysisdictsasdf("maskingdatamodel4", 8, maxrounds, 1135)
insertcumulativeeffects!(ukdataconfounders1leadlagdf, -21:7:21; deltaindex=3:8)
ukdataconfounders2leadlagdf = loadanalysisdictsasdf("maskingdatamodel6", 8, maxrounds, 1155)
insertcumulativeeffects!(ukdataconfounders2leadlagdf, -21:7:21; deltaindex=4:9)
ukdataconfounderslogeffectivereproductionratios = quantilelogeffectivereproductionratios(
    ukdataconfounders2leadlagdf,
    4,
    maskcovidcases,
    POPULATION2020,
    [ 1.0; collect(56.0:28:224); 257 ],
    facialcoveringsrequired,
    [
        [ endstayathometimes, somebusinessreopen, facialcoveringsrecommended ];
        lagleadinterventionsmatrix(facialcoveringsrequired_IM, -21:7:21) 
    ]
)
let 
    logr0a = logbasicreproductionratios(
        ukdataconfounders2leadlagdf, 
        4, 
        [ 1.0; collect(56.0:28:224); 257 ], 
        facialcoveringsrequired, 
        [
            [ endstayathometimes, somebusinessreopen, facialcoveringsrecommended ];
            lagleadinterventionsmatrix(facialcoveringsrequired_IM, -21:7:21) 
        ];
        logdelta=0, 
        logsecondarydelta4=0, 
        logsecondarydelta5=0, 
        logsecondarydelta6=0,
        logsecondarydelta7=0,
        logsecondarydelta8=0, 
        logsecondarydelta9=0, 
    )
    logr0b = logbasicreproductionratios(
        ukdataconfounders2leadlagdf, 
        4, 
        [ 1.0; collect(56.0:28:224); 257 ], 
        facialcoveringsrequired, 
        [
            [ endstayathometimes, somebusinessreopen, facialcoveringsrecommended ];
            lagleadinterventionsmatrix(facialcoveringsrequired_IM, -21:7:21) 
        ];
    )
    global ukdataconfounderscasesdiff_167 = quantilepredictcumulativedifferenceincases(
        COVIDSERIALINTERVAL, 
        ukdataconfounders2leadlagdf, 
        logr0a, 
        logr0b,
        maskcovidcases[1:146, :], 
        POPULATION2020
    )
    global ukdataconfounderscasesdiff_192 = quantilepredictcumulativedifferenceincases(
        COVIDSERIALINTERVAL, 
        ukdataconfounders2leadlagdf, 
        logr0a, 
        logr0b,
        maskcovidcases[1:171, :], 
        POPULATION2020
    )
    global ukdataconfounderscasesdiff_174 = quantilepredictcumulativedifferenceincases(
        COVIDSERIALINTERVAL, 
        ukdataconfounders2leadlagdf, 
        logr0a, 
        logr0b,
        maskcovidcases[1:153, :], 
        POPULATION2020
    )
end
usdataleadlagdf = loadanalysisdictsasdf("datamodelus1leadlag", 12, maxrounds, 105)
insertcumulativeeffects!(usdataleadlagdf, -21:7:21)
usstatedataleadlag = let 
    plotstates = String31[ ]
    plotmedians = Float64[ ] 
    plotlcis = Float64[ ]  
    plotucis = Float64[ ] 

    logr0a = logbasicreproductionratios(
        usdataleadlagdf, 
        51, 
        [ collect(1.0:28:113); [ 123 ] ], 
        maskday, 
        lagleadinterventionsmatrix(maskday, -21:7:21);
        logdelta=0, 
        logsecondarydelta1=0, 
        logsecondarydelta2=0, 
        logsecondarydelta3=0,
        logsecondarydelta4=0,
        logsecondarydelta5=0, 
        logsecondarydelta6=0, 
    )
    logr0b = logbasicreproductionratios(
        usdataleadlagdf, 
        51, 
        [ collect(1.0:28:113); [ 123 ] ], 
        maskday, 
        lagleadinterventionsmatrix(maskday, -21:7:21);
    )

    for (i, state) ∈ enumerate(states.state)
        if maskday.starttimes[i] < 123 
            push!(plotstates, state)
            starttime = maskday.starttimes[i] 
            stminus21 = max(0, starttime-21)

            qs = quantilepredictcumulativedifferenceincases(
                COVIDSERIALINTERVAL, 
                usdataleadlagdf, 
                logr0a, 
                logr0b,
                incidence[1:stminus21, :], 
                populations
            )

            med = qs.med[123, i] * 10_000 / populations[i]
            lci = qs.lci[123, i] * 10_000 / populations[i]
            uci = qs.uci[123, i] * 10_000 / populations[i]

            push!(plotmedians, med)
            push!(plotlcis, lci)
            push!(plotucis, uci)
        end
    end
    @ntuple plotstates plotmedians plotlcis plotucis
end
usdataconfoundersleadlagdf = loadanalysisdictsasdf(
    "datamodelus2leadlag", 12, maxrounds, 115; 
)
insertcumulativeeffects!(usdataconfoundersleadlagdf, -21:7:21; deltaindex=6:11)
usstatedataleadlag_confounder = let 
    plotstates = String31[ ]
    plotmedians = Float64[ ] 
    plotlcis = Float64[ ]  
    plotucis = Float64[ ] 

    logr0a = logbasicreproductionratios(
        usdataconfoundersleadlagdf, 
        51, 
        [ collect(1.0:28:113); [ 123 ] ], 
        maskday, 
        [
            [ relaxshelterinplace, reopenbusiness, reopenrestaurants, reopengyms, reopencinemas ];
            lagleadinterventionsmatrix(maskday, -21:7:21)
        ];
        logdelta=0, 
        logsecondarydelta6=0, 
        logsecondarydelta7=0, 
        logsecondarydelta8=0,
        logsecondarydelta9=0,
        logsecondarydelta10=0, 
        logsecondarydelta11=0, 
    )
    logr0b = logbasicreproductionratios(
        usdataconfoundersleadlagdf, 
        51, 
        [ collect(1.0:28:113); [ 123 ] ], 
        maskday, 
        [
            [ relaxshelterinplace, reopenbusiness, reopenrestaurants, reopengyms, reopencinemas ];
            lagleadinterventionsmatrix(maskday, -21:7:21)
        ];
    )

    for (i, state) ∈ enumerate(states.state)
        if maskday.starttimes[i] < 123 
            push!(plotstates, state)
            starttime = maskday.starttimes[i] 
            stminus21 = max(0, starttime-21)

            qs = quantilepredictcumulativedifferenceincases(
                COVIDSERIALINTERVAL, 
                usdataconfoundersleadlagdf, 
                logr0a, 
                logr0b,
                incidence[1:stminus21, :], 
                populations
            )

            med = qs.med[123, i] * 10_000 / populations[i]
            lci = qs.lci[123, i] * 10_000 / populations[i]
            uci = qs.uci[123, i] * 10_000 / populations[i]

            push!(plotmedians, med)
            push!(plotlcis, lci)
            push!(plotucis, uci)
        end
    end
    @ntuple plotstates plotmedians plotlcis plotucis
end
