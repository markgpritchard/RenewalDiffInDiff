

using DrWatson
@quickactivate :RenewalDiffInDiff

import CSV 
using DataFrames
using Dates
using RenewalDiD

# UK masking data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

maskcoviddf = loadukmaskdata(
    datadir("exp_raw", "OxCGRT_GBR_differentiated_withnotes_2020.csv")
)
observedcases = ukobservedcasesmatrix(maskcoviddf)

# do not keep fitting parameters after all groups are treated: 
#   day 160 for masksrecommended
#   day 227 for masksrequired

# H6E_FacialCoverings == 1 means recommended 
masksrecommended = ukinterventions(maskcoviddf, 160, :H6E_FacialCoverings, 1)  
# H6E_FacialCoverings >= 2 means some level of requirement
masksrequired = ukinterventions(maskcoviddf, 226, :H6E_FacialCoverings, 2) 
# C1E_Schoolclosing >= 2 means require some school closing, 1 means recommend, 0 means no 
# measures
# no groups have school reopening before day 160
schoolreopening_226 = ukinterventions(maskcoviddf, 226, :C1E_Schoolclosing, x -> x <= 1, 2) 
# C6E_Stayathome == 0 means no stay at home measures
endstayathome_160 = ukinterventions(maskcoviddf, 160, :C6E_Stayathome, x -> x == 0, 2)    
endstayathome_226 = ukinterventions(maskcoviddf, 226, :C6E_Stayathome, x -> x == 0, 2)    
# H3E_Contacttracing ==2 means "comprehensive contact tracing; done for all identified cases"
comprehensivecontacttracing_160 = ukinterventions(maskcoviddf, 160, :H3E_Contacttracing, 2)    
comprehensivecontacttracing_226 = ukinterventions(maskcoviddf, 226, :H3E_Contacttracing, 2)    

# do not keep fitting parameters after all groups are treated: 
#   day 160 for masksrecommended
#   day 227 for masksrequired

covidmaskdata1 = RenewalDiDData( ;
    observedcases=observedcases[1:161, :],
    interventions=masksrecommended,
    Ns=UKPOPULATION2020,
    sampletime=30,
    id="data1: masks recommended, no other interventions"
)
safesave(datadir("exp_pro", "covidmaskdata1.jld2"), Dict("data" => covidmaskdata1))

covidmaskdata2 = RenewalDiDData( ;
    observedcases=observedcases[1:227, :],
    interventions=masksrequired,
    Ns=UKPOPULATION2020,
    sampletime=30,
    id="data2: masks required, no other interventions"
)
safesave(datadir("exp_pro", "covidmaskdata2.jld2"), Dict("data" => covidmaskdata2))

covidmaskdata3 = RenewalDiDData( ;
    observedcases=observedcases[1:161, :],
    interventions=InterventionArray(masksrecommended; offset=-35:7:35),
    Ns=UKPOPULATION2020,
    sampletime=30,
    id="data3: masks recommended, lead and lag placebos"
)
safesave(datadir("exp_pro", "covidmaskdata3.jld2"), Dict("data" => covidmaskdata3))

covidmaskdata4 = RenewalDiDData( ;
    observedcases=observedcases[1:227, :],
    interventions=InterventionArray(masksrequired; offset=-35:7:35),
    Ns=UKPOPULATION2020,
    sampletime=30,
    id="data4: masks required, lead and lag placebos"
)
safesave(datadir("exp_pro", "covidmaskdata4.jld2"), Dict("data" => covidmaskdata4))

covidmaskdata5 = RenewalDiDData( ;
    observedcases=observedcases[1:161, :],
    interventions=cat(
        masksrecommended, endstayathome_160, comprehensivecontacttracing_160; 
        dims=3
    ),
    Ns=UKPOPULATION2020,
    sampletime=30,
    id="data5: masks recommended, competing interventions"
)
safesave(datadir("exp_pro", "covidmaskdata5.jld2"), Dict("data" => covidmaskdata5))

covidmaskdata6 = RenewalDiDData( ;
    observedcases=observedcases[1:227, :],
    interventions=cat(
        masksrequired,
        schoolreopening_226, 
        endstayathome_226, 
        comprehensivecontacttracing_226; 
        dims=3
    ),
    Ns=UKPOPULATION2020,
    sampletime=30,
    id="data6: masks required, competing interventions"
)
safesave(datadir("exp_pro", "covidmaskdata6.jld2"), Dict("data" => covidmaskdata6))

covidmaskdata7 = RenewalDiDData( ;
    observedcases=observedcases[1:161, :],
    interventions=(
        iv = InterventionArray(masksrecommended; offset=-35:7:35);
        cat(iv, endstayathome_160, comprehensivecontacttracing_160; dims=3)
    ),
    Ns=UKPOPULATION2020,
    sampletime=30,
    id="data7: masks recommended, competing interventions and lead and lag placebos"
)
safesave(datadir("exp_pro", "covidmaskdata7.jld2"), Dict("data" => covidmaskdata7))

covidmaskdata8 = RenewalDiDData( ;
    observedcases=observedcases[1:227, :],
    interventions=(
        iv = InterventionArray(masksrequired; offset=-35:7:35);
        cat(
            iv, schoolreopening_226, endstayathome_226, comprehensivecontacttracing_226; 
            dims=3
        )
    ),
    Ns=UKPOPULATION2020,
    sampletime=30,
    id="data8: masks required, competing interventions and lead and lag placebos"
)
safesave(datadir("exp_pro", "covidmaskdata8.jld2"), Dict("data" => covidmaskdata8))


# Liverpool testing data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Limit to those local authorities near Liverpool that were in the same tiers for control in
# 2020
const LOCATIONINDEXES = [
    15,  # Halton
    17,  # Knowsley
    19,  # Liverpool
    28,  # Sefton
    31,  # St Helens
    35,  # Warrington
    36,  # West Lancashire
    37,  # Wigan
    38,  # Wirral
]

testcoviddf, liverpoolpopulations = let 
    casesdf1 = CSV.read(
        datadir(
            "exp_raw", "North_West_epidemiological_charts__data_set_2021-01-04_cases.csv"
        ),
        DataFrame;
        header=4,
    )
    casesdf2 = CSV.read(
        datadir(
            "exp_raw", "North_West_epidemiological_charts__data_set_2021-06-07_cases.csv"
        ),
        DataFrame;
        header=4,
    )
    testsdf1 = CSV.read(
        datadir(
            "exp_raw", "North_West_epidemiological_charts__data_set_2021-01-04_tests.csv"
        ),
        DataFrame;
        header=4,
    )
    testsdf2 = CSV.read(
        datadir(
            "exp_raw", "North_West_epidemiological_charts__data_set_2021-06-07_tests.csv"
        ),
        DataFrame;
        header=4,
    )
    popdf = CSV.read(
        datadir(
            "exp_raw", 
            "North_West_epidemiological_charts__data_set_2021-06-07_population.csv"
        ),
        DataFrame;
        header=4,
    )
    populations = getproperty(popdf, "Total population")

    for df ∈ [casesdf1, testsdf1]
        insertcols!(df, :date => Date.(df.Date, "dd u yyyy"))
    end

    for df ∈ [casesdf2, testsdf2]
        insertcols!(df, :date => Date.(getproperty(df, "Specimen date"), "dd u yyyy"))
    end

    for df ∈ [casesdf1, casesdf2, testsdf1, testsdf2, popdf]
        rename!(df, Dict(Symbol("Local Authority") => "LocalAuthority"))
    end

    for df ∈ [casesdf1, casesdf2]
        rename!(df, Dict(Symbol("Total Cases") => "cases"))
        rename!(df, Dict(Symbol("Pillar 1 Cases") => "pillar1cases"))
        select!(df, :date, :LocalAuthority, :cases, :pillar1cases)
    end

    for df ∈ [testsdf1, testsdf2]
        rename!(df, Dict(Symbol("Persons tested, 7-day moving average") => "tests"))
        select!(df, :date, :LocalAuthority, :tests)
    end

    rename!(popdf, Dict(Symbol("Total population") => "population"))
    select!(popdf, :LocalAuthority, :population)
    insertcols!(popdf, :location => 1:39)

    # remove cases after 1 December from `casesdf1` and `testsdf1`
    for df ∈ [casesdf1, testsdf1] 
        filter!(:date => x -> x < Date("2020-12-01"), df)
    end

    casesdf = vcat(casesdf1, casesdf2)
    testsdf = vcat(testsdf1, testsdf2)

    joineddf = outerjoin(casesdf, testsdf; on=[:date, :LocalAuthority])
    insertcols!(joineddf, :day => Dates.value.(joineddf.date .- Date("2020-05-31")))
    df = leftjoin(joineddf, popdf; on=:LocalAuthority)

    (df, populations)
end

# universal testing was introduced in January 2021, so remove all rows after this 
filter!(:date => x -> x < Date("2021-01-03"), testcoviddf)

testingstarttimes = [ 
    i == 3 ?  # is Liverpool: start date is 7 November 2020 
        Dates.value(Date("2020-11-07") - Date("2020-05-31")) : 
        i ∈ [ 1, 2, 4, 5, 9 ] ?  # areas where testing introduced on 3 December 
            Dates.value(Date("2020-12-03") - Date("2020-05-31")) :
            nothing  # places where the testing programme was not introduced in 2020
    for i ∈ 1:9 
]

## Convert DataFrame to appropriate matrices 
allcovidcases, pil1covidcases = let 
    # check that each location has the same number of rows 
    for i ∈ 2:39 
        @assert sum(testcoviddf.location .== 1) == sum(testcoviddf.location .== i) 
    end
    
    # how many rows is it?
    covidlength = sum(testcoviddf.location .== 1)
    #304

    allcovidcases = Matrix{Int}(undef, covidlength, 9)
    pil1covidcases = Matrix{Int}(undef, covidlength, 9)

    for i ∈ 1:9
        k = LOCATIONINDEXES[i]
        _tdf = filter(:location => x -> x == k, testcoviddf)
        for j ∈ 1:covidlength
            allcovidcases[j, i] = _tdf.cases[j]
            pil1covidcases[j, i] = _tdf.pillar1cases[j]
        end
    end 

    (allcovidcases, pil1covidcases)
end
selectpops = [liverpoolpopulations[x] for x ∈ LOCATIONINDEXES]

testingintervention = InterventionMatrix(215, testingstarttimes)

covidtestingdata1 = RenewalDiDData( ;
    observedcases=allcovidcases,
    interventions=testingintervention,
    Ns=selectpops,
    sampletime=14,
    id="data1: all cases, no placebo interventions"
)
safesave(datadir("exp_pro", "covidtestingdata1.jld2"), Dict("data" => covidtestingdata1))

covidtestingdata2 = RenewalDiDData( ;
    observedcases=pil1covidcases,
    interventions=testingintervention,
    Ns=selectpops,
    sampletime=14,
    id="data2: pillar 1 cases, no placebo interventions"
)
safesave(datadir("exp_pro", "covidtestingdata2.jld2"), Dict("data" => covidtestingdata2))

covidtestingdata3 = RenewalDiDData( ;
    observedcases=allcovidcases,
    interventions=InterventionArray(testingintervention; offset=-35:7:35),
    Ns=selectpops,
    sampletime=14,
    id="data3: all cases, lead and lag placebos"
)
safesave(datadir("exp_pro", "covidtestingdata3.jld2"), Dict("data" => covidtestingdata3))

covidtestingdata4 = RenewalDiDData( ;
    observedcases=pil1covidcases,
    interventions=InterventionArray(testingintervention; offset=-35:7:35),
    Ns=selectpops,
    sampletime=14,
    id="data4: pillar 1 cases, lead and lag placebos"
)
safesave(datadir("exp_pro", "covidtestingdata4.jld2"), Dict("data" => covidtestingdata4))
