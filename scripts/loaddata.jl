#=
coviddf = let 
    df = CSV.read(
        datadir("exp_raw", "OxCGRT_GBR_differentiated_withnotes_2020.csv"), DataFrame
    )
#=
    regioncodedf = DataFrame(
        :RegionCode => [ "UK_ENG", "UK_NIR", "UK_SCO", "UK_WAL" ],
        :RegionId => 1:4
    )=#

    # Government Response Index includes C1, C2, C3, C4, C5, C6, C7, C8, E1, E2, H1, H2, H3,
    # H6, H7, H8 (NB H6 is face coverings)

    # format date as `Date`
    rename!(
        df, 
        Dict(
            :Date => "_date",  # rename old version of date 
            Symbol("H6E_Facial Coverings") => "H6E_FacialCoverings",  # remove space in name 
            Symbol("C1E_School closing") => "C1E_Schoolclosing",
            Symbol("C2E_Workplace closing") => "C2E_Workplaceclosing",
            Symbol("C3E_Cancel public events") => "C3E_Cancelpublicevents",
            Symbol("C4E_Restrictions on gatherings") => "C4E_Restrictionsongatherings",
            Symbol("C5E_Close public transport") => "C5E_Closepublictransport",
            Symbol("C6E_Stay at home requirements") => "C6E_Stayathome",
            Symbol("C7E_Restrictions on internal movement") => 
                "C7E_Restrictionsoninternalmovement",
            Symbol("C8E_International travel controls") => "C8E_Internationaltravelcontrols",
            Symbol("E1E_Income support") => "E1E_Incomesupport",
            Symbol("E2E_Debt/contract relief") => "E2E_Debtcontractrelief",
            Symbol("H1E_Public information campaigns") => "H1E_Publicinformationcampaigns",   
            Symbol("H2E_Testing policy") => "H2E_Testingpolicy",
            Symbol("H3E_Contact tracing") => "H3E_Contacttracing",
            Symbol("H7E_Vaccination policy") => "H7E_Vaccinationpolicy",
            Symbol("H8E_Protection of elderly people") => "H8E_Protectionofelderlypeople",
        )
    ) 
    insertcols!(df, :Date => [ Date("$d", "yyyymmdd") for d ∈ df._date ])

    # choose necessary columns 
    select!(
        df,
        :RegionName,
        :RegionCode,
        :Jurisdiction,
        :Date,
        :H6E_FacialCoverings,
        :H6E_Flag,
        :C1E_Schoolclosing,
        :C2E_Workplaceclosing,
        :C3E_Cancelpublicevents,
        :C4E_Restrictionsongatherings,
        :C5E_Closepublictransport,
        :C6E_Stayathome,
        :C7E_Restrictionsoninternalmovement,
        :C8E_Internationaltravelcontrols,
        :E1E_Incomesupport,
        :E2E_Debtcontractrelief,
        :H1E_Publicinformationcampaigns,
        :H2E_Testingpolicy,
        :H3E_Contacttracing,
        :H7E_Vaccinationpolicy,
        :H8E_Protectionofelderlypeople,
        :ConfirmedCases,
        :ConfirmedDeaths,
        :GovernmentResponseIndex_SimpleAverage_ForDisplay,
    )

    # H6E_FacialCoverings is missing for each nation on 2020-09-24, but == 3 for UK as a whole 
    for i ∈ axes(df, 1)
        if df.Date[i] == Date("2020-09-24") df.H6E_FacialCoverings[i] = 3 end
    end

    # convert cumulative cases to daily new cases 
    _calcdiff(a::Real, b::Real) = b - a 
    _calcdiff(::Missing, b::Real) = b 
    _calcdiff(::Any, ::Missing) = 0  # includes first case for each country

    insertcols!(
        df,
        :NewConfirmedCases => [ 
            i == 1 ? 0 : _calcdiff(df.ConfirmedCases[i-1], df.ConfirmedCases[i]) 
            for i ∈ axes(df, 1) 
        ],
        :NewConfirmedDeaths => [ 
            i == 1 ? 0 : _calcdiff(df.ConfirmedDeaths[i-1], df.ConfirmedDeaths[i]) 
            for i ∈ axes(df, 1) 
        ],
    )

    # choose only "state-level" data 
    filter!(:Jurisdiction => x -> x == "STATE_TOTAL", df)

    # dates that facial covering guidance changed 
    _tdf = deepcopy(df)
    insertcols!(
        _tdf,
        :RegulationChange => [ 
            i == 1 ? true : _tdf.H6E_FacialCoverings[i] != _tdf.H6E_FacialCoverings[i-1] 
            for i ∈ axes(_tdf, 1) 
        ]  
    )
    filter!(:RegulationChange => x -> x, _tdf)

    #julia> _tdf
    #16×10 DataFrame
    # Row │ RegionName        RegionCode  Jurisdiction  Date        H6E_FacialCoverings  ConfirmedCases  ConfirmedDeaths  NewConfirmedCases  NewConfirmedDeaths  RegulationChange 
    #     │ String31?         String7?    String15      Date        Float64?             Int64?          Int64?           Int64              Int64               Bool
    #─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    #   1 │ England           UK_ENG      STATE_TOTAL   2020-01-01                  0.0         missing          missing                  0                   0              true 
    #   2 │ England           UK_ENG      STATE_TOTAL   2020-05-12                  1.0          185781            28702               2940                 539              true 
    #   3 │ England           UK_ENG      STATE_TOTAL   2020-06-15                  2.0          234576            34958                944                  25              true 
    #   4 │ England           UK_ENG      STATE_TOTAL   2020-08-26                  3.0          289136            36818               1181                  13              true
    #   5 │ Northern Ireland  UK_NIR      STATE_TOTAL   2020-01-01                  0.0         missing          missing                  0                   0              true
    #   6 │ Northern Ireland  UK_NIR      STATE_TOTAL   2020-07-10                  2.0            4774              554                  3                   0              true
    #   7 │ Northern Ireland  UK_NIR      STATE_TOTAL   2020-09-24                  3.0            9250              577                290                   0              true
    #   8 │ Northern Ireland  UK_NIR      STATE_TOTAL   2020-09-25                  2.0            9499              577                249                   0              true
    #   9 │ Northern Ireland  UK_NIR      STATE_TOTAL   2020-10-14                  3.0           23211              602               1071                   4              true
    #  10 │ Scotland          UK_SCO      STATE_TOTAL   2020-01-01                  0.0         missing          missing                  0                   0              true
    #  11 │ Scotland          UK_SCO      STATE_TOTAL   2020-04-28                  1.0           11543             1332                383                  70              true
    #  12 │ Scotland          UK_SCO      STATE_TOTAL   2020-06-22                  2.0           18207             2472                 10                   0              true
    #  13 │ Scotland          UK_SCO      STATE_TOTAL   2020-08-17                  3.0           19489             2491                 72                   0              true
    #  14 │ Wales             UK_WAL      STATE_TOTAL   2020-01-01                  0.0         missing          missing                  0                   0              true
    #  15 │ Wales             UK_WAL      STATE_TOTAL   2020-06-09                  1.0           15173             1410                 62                   9              true
    #  16 │ Wales             UK_WAL      STATE_TOTAL   2020-09-14                  3.0           20139             1597                192                   0              true

    # H6E_FacialCoverings == 1 means recommended 
    # 2 <= H6E_FacialCoverings <= 4 means some level of requirement

    # Northern Ireland never had a separaete "recommended" phase so classify recommended from 
    # 10 July. This is the latest recommendation, so stop the analysis of recommended on 9 July. 
    # The mandate in Wales on 14 September was the latest start of "required", so stop analysis 
    # of required on 13 September.

    insertcols!(
        df,
        :FacialCoveringRecommended => [ 
            df.H6E_FacialCoverings[i] >= 1 ? 1 : 0 
            for i ∈ axes(df, 1) 
        ],
        :FacialCoveringRequired => [ 
            df.H6E_FacialCoverings[i] >= 2 ? 1 : 0 
            for i ∈ axes(df, 1) 
        ],
    )

    # the contribution of facial coverings to the government response index
    fcgri = [ 
        h6 == 0 ? 
            0.0 : 
            100 * (h6 - 0.5 * (1 - flag)) / 4  
        for (h6, flag) ∈ zip(df.H6E_FacialCoverings, df.H6E_Flag)
    ] 

    # government response index without the contribution of facial coverings
    gri_nofc = df.GovernmentResponseIndex_SimpleAverage_ForDisplay .- fcgri ./ 13
    insertcols!(df, :Gri_nofc => gri_nofc)

    # remove all dates after 13 September from DataFrame 
    filter!(:Date => x -> x <= Date("2020-09-13"), df)

    #unique(filter(:ConfirmedCases => x -> ismissing(x), df).Date)
    # last day when there are no cases is 27 February 
    #filter!(:Date => x -> x >= Date("2020-02-28"), df)

    # add numeric code for region 
    leftjoin!(df, REGIONCODEDF; on=:RegionCode)

    df
end

=#





coviddf = let 
    casesdf1 = CSV.read(
        datadir("exp_raw", "North_West_epidemiological_charts__data_set_2021-01-04_cases.csv"),
        DataFrame;
        header=4,
    )
    casesdf2 = CSV.read(
        datadir("exp_raw", "North_West_epidemiological_charts__data_set_2021-06-07_cases.csv"),
        DataFrame;
        header=4,
    )
    testsdf1 = CSV.read(
        datadir("exp_raw", "North_West_epidemiological_charts__data_set_2021-01-04_tests.csv"),
        DataFrame;
        header=4,
    )
    testsdf2 = CSV.read(
        datadir("exp_raw", "North_West_epidemiological_charts__data_set_2021-06-07_tests.csv"),
        DataFrame;
        header=4,
    )
    popdf = CSV.read(
        datadir("exp_raw", "North_West_epidemiological_charts__data_set_2021-06-07_population.csv"),
        DataFrame;
        header=4,
    )
    global populations = getproperty(popdf, "Total population")

    for df ∈ [ casesdf1, testsdf1 ]
        insertcols!(df, :date => Date.(df.Date, "dd u yyyy"))
    end

    for df ∈ [ casesdf2, testsdf2 ]
        insertcols!(df, :date => Date.(getproperty(df, "Specimen date"), "dd u yyyy"))
    end

    for df ∈ [ casesdf1, casesdf2, testsdf1, testsdf2, popdf ]
        rename!(df, Dict(Symbol("Local Authority") => "LocalAuthority"))
    end

    for df ∈ [ casesdf1, casesdf2 ]
        rename!(df, Dict(Symbol("Total Cases") => "cases"))
        rename!(df, Dict(Symbol("Pillar 1 Cases") => "pillar1cases"))
        select!(df, :date, :LocalAuthority, :cases, :pillar1cases)
    end

    for df ∈ [ testsdf1, testsdf2 ]
        rename!(df, Dict(Symbol("Persons tested, 7-day moving average") => "tests"))
        select!(df, :date, :LocalAuthority, :tests)
    end

    rename!(popdf, Dict(Symbol("Total population") => "population"))
    select!(popdf, :LocalAuthority, :population)
    insertcols!(popdf, :location => 1:39)

    # remove cases after 1 December from `casesdf1` and `testsdf1`
    for df ∈ [ casesdf1, testsdf1 ] 
        filter!(:date => x -> x < Date("2020-12-01"), df)
    end

    casesdf = vcat(casesdf1, casesdf2)
    testsdf = vcat(testsdf1, testsdf2)

    joineddf = outerjoin(casesdf, testsdf; on=[ :date, :LocalAuthority ])
    insertcols!(joineddf, :day => Dates.value.(joineddf.date .- Date("2020-05-31")))
    leftjoin(joineddf, popdf; on=:LocalAuthority)
end
#=
# universal testing was introduced in April 2021, so remove all rows after this 
filter!(:date => x -> x < Date("2021-04-01"), coviddf)
=#

# universal testing was introduced in January 2021, so remove all rows after this 
filter!(:date => x -> x < Date("2021-01-03"), coviddf)
