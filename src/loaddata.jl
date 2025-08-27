

# constants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

const _MASKINGCOLUMNNAMES = [
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
]

# Government Response Index includes C1, C2, C3, C4, C5, C6, C7, C8, E1, E2, H1, H2, H3, H6, 
# H7, H8 (NB H6 is face coverings)
const _UNSPACEDMASKINGCOLUMNNAMES = Dict(
    Symbol("H6E_Facial Coverings") => "H6E_FacialCoverings",   
    Symbol("C1E_School closing") => "C1E_Schoolclosing",
    Symbol("C2E_Workplace closing") => "C2E_Workplaceclosing",
    Symbol("C3E_Cancel public events") => "C3E_Cancelpublicevents",
    Symbol("C4E_Restrictions on gatherings") => "C4E_Restrictionsongatherings",
    Symbol("C5E_Close public transport") => "C5E_Closepublictransport",
    Symbol("C6E_Stay at home requirements") => "C6E_Stayathome",
    Symbol("C7E_Restrictions on internal movement") => "C7E_Restrictionsoninternalmovement",
    Symbol("C8E_International travel controls") => "C8E_Internationaltravelcontrols",
    Symbol("E1E_Income support") => "E1E_Incomesupport",
    Symbol("E2E_Debt/contract relief") => "E2E_Debtcontractrelief",
    Symbol("H1E_Public information campaigns") => "H1E_Publicinformationcampaigns",   
    Symbol("H2E_Testing policy") => "H2E_Testingpolicy",
    Symbol("H3E_Contact tracing") => "H3E_Contacttracing",
    Symbol("H7E_Vaccination policy") => "H7E_Vaccinationpolicy",
    Symbol("H8E_Protection of elderly people") => "H8E_Protectionofelderlypeople",
)

# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## UK masking data

loadukmaskdata(filename; startdate="2020-01-31") = _loadukmaskdata(filename, startdate)

_loadukmaskdata(filename, startdate::String) = _loadukmaskdata(filename, Date(startdate))

function _loadukmaskdata(filename, startdate::Date)
    df = CSV.read(filename, DataFrame)
    rename!(df, _UNSPACEDMASKINGCOLUMNNAMES)  # remove spaces in names
    rename!(df, Dict(:Date => "_date"))  # rename old version of date 
    insertcols!(df, :Date => [ Date("$d", "yyyymmdd") for d ∈ df._date ])  # new date 
    select!(df, _MASKINGCOLUMNNAMES...)

    # H6E_FacialCoverings is missing for each nation on 2020-09-24, but == 3 for UK as a whole 
    for i ∈ axes(df, 1)
        if df.Date[i] == Date("2020-09-24") df.H6E_FacialCoverings[i] = 3 end
    end

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
    # remove first 30 days when there are no infections 
    filter!(:Date => x -> x >= startdate, df)

    return df
end

function ukobservedcasesmatrix(df)
    cases = zeros(Int, 336, 4)
    for (g, c) in enumerate(UKNATIONS)
        _tdf = filter(:RegionName => x -> x == c, df)
        for t in axes(_tdf, 1)
            cases[t, g] = _tdf.NewConfirmedCases[t]
        end
    end
    return cases
end

function ukinterventions(df, interventioncolumn::Symbol, interventionvalue, eventindex=Val(1))
    return ukinterventions(df, 335, interventioncolumn, interventionvalue, eventindex)
end

function ukinterventions(
    df, duration::Integer, interventioncolumn::Symbol, interventionvalue, eventindex=Val(1)
)
    return _ukinterventions(df, duration, interventioncolumn, interventionvalue, eventindex)
end

function _ukinterventions(df, duration, interventioncolumn, interventionvalue, eventindex)
    return InterventionMatrix(
        duration,
        [
            (
                x = _ukinterventionsfunction(
                    interventionvalue, getproperty(tdf, interventioncolumn), eventindex
                );
                _sub1(x)
            )
            for tdf in [filter(:RegionName => x -> x == c, df) for c in UKNATIONS]
        ]
    )
end

function _ukinterventionsfunction(interventionvalue::Number, vec, eventindex)
    return _ukinterventionsfunction(x -> x >= interventionvalue, vec, eventindex)
end

_ukinterventionsfunction(f::Function, vec, ::Val{1}) = findfirst(f, vec)

function _ukinterventionsfunction(f::Function, vec, eventindex::Integer)
    eventindex == 1 && return _ukinterventionsfunction(f, vec, Val(1))
    return _findnth(f, vec, eventindex)
end

_calcdiff(a::Real, b::Real) = b - a 
_calcdiff(::Missing, b::Real) = b 
_calcdiff(::Any, ::Missing) = 0  # includes first case for each country

_sub1(::Nothing) = nothing 
_sub1(x::Real) = x - 1

function _findnth(f, v, n::Integer)
    n > 1 || throw(ArgumentError("$n: n must be greater than 1; for n == 1 use `findfirst`"))
    allindexes = findall(f, v)
    counter = 1
    for (j, index) in enumerate(allindexes)
        j == 1 && continue 
        allindexes[j] == allindexes[j-1] + 1 && continue  # consequtive indexes counted as one
        counter += 1
        counter == n && return index
    end
    return nothing
end
