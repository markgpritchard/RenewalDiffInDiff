
function insertcumulativeeffects!(df; leadlagtimes, kwargs...)
    insertcumulativeeffects!(df, leadlagtimes; kwargs...)
end

function insertcumulativeeffects!(df, leadlagtimes; deltaindex=automatic, kwargs...)
    insertcumulativeeffects!(df, leadlagtimes, deltaindex; kwargs...)
end

function insertcumulativeeffects!(df, leadlagtimes, ::Automatic; kwargs...)
    deltaindex = eachindex(findall(_notzero, leadlagtimes))
    insertcumulativeeffects!(df, leadlagtimes, deltaindex; kwargs...)
end

function insertcumulativeeffects!(df, leadlagtimes, deltaindex; kwargs...)
    previoust = -1 
    previousphi = nothing
    for (i, t) ∈ enumerate(leadlagtimes)
        if t == 0
            # insert the value for the time of the intervention 
            newphi = :phi_0
            _insertcumulativeeffects!(df, newphi, previousphi, :logdelta)
        else 
            if t < 0 
                newphi = Symbol("phi_minus$(abs(t))")
                ind = deltaindex[i]
            else 
                newphi = Symbol("phi_plus$(abs(t))") 
                ind = deltaindex[i-1]
            end
            _insertcumulativeeffects!(df, newphi, previousphi, ind)
        end
        previoust = t
        previousphi = newphi
    end
end

function _insertcumulativeeffects!(df, newphi, previousphi, i::Int)
    colname = Symbol("logsecondarydelta$i.logsecondarydelta")
    _insertcumulativeeffects!(df, newphi, previousphi, colname)
end

function _insertcumulativeeffects!(df, newphi, ::Nothing, colname::Symbol)
    insertcols!(df, newphi => getproperty(df, colname))
end

function _insertcumulativeeffects!(df, newphi, previousphi, colname::Symbol)
    insertcols!(df, newphi => getproperty(df, previousphi) .+ getproperty(df, colname))
end

function phiquantiles(df, leadlagtimes; deltaindex=automatic, kwargs...)
    return _phiquantiles(df, leadlagtimes, deltaindex; kwargs...)
end 

function _phiquantiles(df, leadlagtimes, ::Automatic; kwargs...)
    deltaindex = eachindex(leadlagtimes)
    return _phiquantiles(df, leadlagtimes, deltaindex; kwargs...)
end 

function _phiquantiles(df, leadlagtimes, deltaindex; CrI=( 0.05, 0.95 ))
    lci = zeros(length(leadlagtimes))
    med = zeros(length(leadlagtimes))
    uci = zeros(length(leadlagtimes))
    for (i, t) ∈ zip(deltaindex, leadlagtimes) 
        if t < 0 
            phisymbol = Symbol("phi_minus$(abs(t))")
        elseif t == 0 
            phisymbol = :phi_0
        else
            phisymbol = Symbol("phi_plus$(abs(t))") 
        end
        lcii, medi, ucii = quantile(getproperty(df, phisymbol), [ CrI[1], 0.5, CrI[2] ])
        lci[i] = lcii 
        med[i] = medi 
        uci[i] = ucii 
    end 
    return @ntuple lci med uci 
end 

function quantilelogeffectivereproductionratios(
    df, nlocations, cases, ns, timeknots, interventions, secondaryinterventions=nothing;
    times=automatic, kwargs...
)
    return _quantilelogeffectivereproductionratios(
        df, nlocations, cases, ns, timeknots, interventions, secondaryinterventions, times;
        kwargs...
    )
end

function _quantilelogeffectivereproductionratios(
    df, nlocations, cases, ns, timeknots, interventions, secondaryinterventions, ::Automatic;
    kwargs...
)
    times = round(Int, timeknots[1]):1:round(Int, last(timeknots))
    return _quantilelogeffectivereproductionratios(
        df, nlocations, cases, ns, timeknots, interventions, secondaryinterventions, times;
        kwargs...
    )
end

function _quantilelogeffectivereproductionratios(
    df, nlocations, cases, ns, timeknots, interventions, secondaryinterventions, times;
    quantiles=[ 0.05, 0.5, 0.95 ], kwargs...
)
    logre = logeffectivereproductionratios(
        df, nlocations, cases, ns, timeknots, interventions, secondaryinterventions;
        times, kwargs...
    )

    lci = zeros(length(times), nlocations)
    med = zeros(length(times), nlocations)
    uci = zeros(length(times), nlocations)

    Threads.@threads for j ∈ 1:nlocations
        for t ∈ eachindex(times) 
            lcii, medi, ucii = quantile(logre[:, t, j], quantiles)
            lci[t, j] = lcii 
            med[t, j] = medi 
            uci[t, j] = ucii 
        end
    end 

    return @ntuple lci med uci
end

function quantilepredictcases(
    g, df, logR0, initialcases, ns; quantiles=[ 0.05, 0.5, 0.95 ], 
    kwargs...
)
    cases = predictcases(g, df, logR0, initialcases, ns; kwargs...)

    lci = zeros(size(cases, 2), size(cases, 3))
    med = zeros(size(cases, 2), size(cases, 3))
    uci = zeros(size(cases, 2), size(cases, 3))

    Threads.@threads for j ∈ axes(cases, 3)
        for t ∈ axes(cases, 2) 
            lcii, medi, ucii = quantile(cases[:, t, j], quantiles)
            lci[t, j] = lcii 
            med[t, j] = medi 
            uci[t, j] = ucii 
        end
    end 

    return @ntuple lci med uci
end

function quantilepredictdifferenceincases(
    g, df, logR0a, logR0b, initialcases, ns; quantiles=[ 0.05, 0.5, 0.95 ], 
    kwargs...
)
    casesa = predictcases(g, df, logR0a, initialcases, ns; kwargs...)
    casesb = predictcases(g, df, logR0b, initialcases, ns; kwargs...)
    casesdiff = casesb .- casesa

    lci = zeros(size(casesdiff, 2), size(casesdiff, 3))
    med = zeros(size(casesdiff, 2), size(casesdiff, 3))
    uci = zeros(size(casesdiff, 2), size(casesdiff, 3))

    Threads.@threads for j ∈ axes(casesdiff, 3)
        for t ∈ axes(casesdiff, 2) 
            lcii, medi, ucii = quantile(casesdiff[:, t, j], quantiles)
            lci[t, j] = lcii 
            med[t, j] = medi 
            uci[t, j] = ucii 
        end
    end 

    return @ntuple lci med uci
end

@memoize function quantilepredictcumulativedifferenceincases(
    g, df, logR0a, logR0b, initialcases, ns; 
    quantiles=[ 0.05, 0.5, 0.95 ], kwargs...
)
    casesa = predictcases(g, df, logR0a, initialcases, ns; kwargs...)
    casesb = predictcases(g, df, logR0b, initialcases, ns; kwargs...)
    casesdiff = similar(casesa)
    lci = zeros(size(casesdiff, 2), size(casesdiff, 3))
    med = zeros(size(casesdiff, 2), size(casesdiff, 3))
    uci = zeros(size(casesdiff, 2), size(casesdiff, 3))

    Threads.@threads for j ∈ axes(casesdiff, 3)
        for t ∈ axes(casesdiff, 2)
            for k ∈ axes(casesdiff, 1)
                if t == 1 
                    casesdiff[k, t, j] = casesb[k, t, j] - casesa[k, t, j]
                else
                    casesdiff[k, t, j] = (
                        sum(@view casesb[k, 1:t, j]) - sum(@view casesa[k, 1:t, j])
                    )
                end
            end
            lcii, medi, ucii = quantile(casesdiff[:, t, j], quantiles)
            lci[t, j] = lcii 
            med[t, j] = medi 
            uci[t, j] = ucii 
        end
    end 

    return @ntuple lci med uci
end

function logeffectivereproductionratio(
    logzeta_g_vec::AbstractVector,
    timespline,
    logdelta,
    interventions, 
    s,
    t, 
    g
)
    return logeffectivereproductionratio(
        logzeta_g_vec[g],
        timespline,
        logdelta,
        interventions, 
        s,
        t,
        g
    )
end

function logeffectivereproductionratio(
    zeta,
    timespline,
    logdelta,
    interventions, 
    s,
    t,
    g
)
    return +(
        zeta, 
        timespline[t], 
        logdelta * interventions[t, g], 
        s
    ) 
end

function logeffectivereproductionratio(
    logzeta_g_vec,
    timespline,
    logdelta,
    interventions, 
    s,
    logsecondarydelta,
    secondaryinterventions,
    t, 
    g
)
    return logeffectivereproductionratio(
        logzeta_g_vec,
        timespline,
        logdelta,
        interventions, 
        s,
        t, 
        g
    ) + _logeffectivereproductionratiosecondaryinterventions(
        logsecondarydelta, secondaryinterventions, t, g
    )
end

_logeffectivereproductionratiosecondaryinterventions(::Any, ::Nothing, ::Any, ::Any) = 0

function _logeffectivereproductionratiosecondaryinterventions(
    logsecondarydelta, secondaryinterventions::AbstractMatrix, t, g
)
    return logsecondarydelta * secondaryinterventions[t, g] 
end

function _logeffectivereproductionratiosecondaryinterventions(
    logsecondarydelta::AbstractVector, secondaryinterventions::AbstractVector, t, g
)
    return sum(
        [ 
            logsecondarydelta[i] * secondaryinterventions[i][t, g] 
            for i ∈ eachindex(secondaryinterventions) 
        ]
    )
end

function logbasicreproductionratios(
    df, nlocations, timeknots, interventions, secondaryinterventions;
    kwargs...
)
    # get the basic reproduction ratio by inputting `nothing` for cases 
    return logeffectivereproductionratios(
        df, nlocations, nothing, nothing, timeknots, interventions, secondaryinterventions;
        kwargs...
    )
end

function logeffectivereproductionratios(
    df, nlocations, cases, ns, timeknots, interventions, secondaryinterventions=nothing;
    times=automatic, kwargs...
)
    return _logeffectivereproductionratios(
        df, nlocations, cases, ns, timeknots, interventions, secondaryinterventions, times;
        kwargs...
    )
end

function _logeffectivereproductionratios(
    df, nlocations, cases, ns, timeknots, interventions, secondaryinterventions, ::Automatic;
    kwargs...
)
    times = round(Int, timeknots[1]):1:round(Int, last(timeknots))
    return _logeffectivereproductionratios(
        df, nlocations, cases, ns, timeknots, interventions, secondaryinterventions, times;
        kwargs...
    )
end

function _logeffectivereproductionratios(
    df, nlocations, cases, ns, timeknots, interventions, secondaryinterventions, times;
    kwargs...
)
    nsamples = size(df, 1)
    ntimes = length(times)
    logre = Array{Float64, 3}(undef, nsamples, ntimes, nlocations)

    Threads.@threads for j ∈ 1:nlocations
        logetavector = zeros(length(timeknots))
        for i ∈ 1:nsamples 
            _logeffectivereproductionratio!(
                logre,
                logetavector, 
                df, 
                cases, 
                ns, 
                timeknots, 
                interventions, 
                secondaryinterventions, 
                times, 
                i, 
                j; 
                kwargs...
            )
        end
    end
    return logre
end
#=
_initiallogsecondarydeltavector(::Nothing) = nothing
_initiallogsecondarydeltavector(::AbstractMatrix) = zeros(1)

function _initiallogsecondarydeltavector(secondaryinterventions::Vector{<:AbstractMatrix})
    return zeros(length(secondaryinterventions))
end
=#
function _logeffectivereproductionratio!(
    logre,
    logetavector, 
    df, 
    cases, 
    Ns, 
    timeknots, 
    interventions, 
    secondaryinterventions, 
    times, 
    i, 
    j;
    extrapl=zeros(1), extrapr=zeros(1),
    zeroperiod=2,
    kwargs...
)
    intercept = getproperty(df, :intercept)[i]
    zeta_sigma2 = getproperty(df, :zeta_sigma2)[i]
    zeta = intercept + getproperty(df, Symbol("logzeta_g$j.logzeta"))[i] * zeta_sigma2
    for k ∈ eachindex(logetavector)
        if isa(zeroperiod, Number) && k == zeroperiod 
            if k == 1 
                _otherk = 2 
            else 
                _otherk = k - 1
            end
            logetavector[k] = zero(getproperty(df, Symbol("eta_t$_otherk.logeta"))[i])
        else
            logetavector[k] = getproperty(df, Symbol("eta_t$k.logeta"))[i]
        end
    end
    eta_sigma2 = getproperty(df, :eta_sigma2)[i]
    timespline = CubicSpline(timeknots, logetavector * eta_sigma2; extrapl, extrapr)
    logdelta = getproperty(df, :logdelta)[i]
    logsecondarydelta = _logsecondarydeltav(df, secondaryinterventions, i)
    psi = getproperty(df, :psi)[i]

    for t ∈ times 
        s = _calcs(cases, Ns, psi, t, j)
        logre[i, t, j] = logeffectivereproductionratio(
            zeta,
            timespline,
            logdelta,
            interventions, 
            s,
            logsecondarydelta,
            secondaryinterventions, 
            t,
            j
        )  
        #println("logre[$t] = $(logre[t])")      
    end
end

function _calcs(cases, Ns, psi, t, j)
    if t == 1 
        return zero(log(1 / (Ns[j] * psi))) 
    else 
        return log(1 - sum(cases[1:(t - 1), j]) / (Ns[j] * psi))
    end
end

_calcs(::Nothing, ::Nothing, ::Any, ::Any, ::Any) = 0

_logsecondarydeltav(::Any, ::Nothing, ::Any) = nothing
_logsecondarydeltav(df, ::AbstractMatrix, i) = getproperty(df, :logsecondarydelta)[i]

function _logsecondarydeltav(df, secondaryinterventions::Vector{<:AbstractMatrix}, i)
    ℓ = length(secondaryinterventions)
    return [ 
        getproperty(df, Symbol("logsecondarydelta$k.logsecondarydelta"))[i] 
        for k ∈ 1:ℓ 
    ]
end

#=
function _logeffectivereproductionratiocumulativecasesvalue(
    ::Any, ::Nothing, ::Any, ::Any, ::Any
)
    # get the basic reproduction ratio by inputting `nothing` for cases 
    return 0
end

function _logeffectivereproductionratiocumulativecasesvalue(
    identifiedproportion, cases, ns, j, m
)
    return log(1 - sum(@view cases[1:m, j]) / (identifiedproportion * ns[j]))
end

function _logeffectivereproductionratiosecondaryinterventions(
    ::Nothing, ::Any, ::Any, ::Any;
)
    return 0
end

function _logeffectivereproductionratiosecondaryinterventions(
    secondaryinterventions::AbstractMatrix, logsecondarydeltavector, j, m;
)
    return _logeffectivereproductionratiosecondaryinterventions(
        [ secondaryinterventions ], logsecondarydeltavector, j, m;
    )
end

function _logeffectivereproductionratiosecondaryinterventions(
    secondaryinterventions::Vector{<:AbstractMatrix}, logsecondarydeltavector, j, m;
)
    return sum(
        [ 
            logsecondarydeltavector[k] * secondaryinterventions[k][m, j] 
            for k ∈ eachindex(secondaryinterventions) 
        ]
    )     
end

function _logeffectivereproductionratio_timespline!(
    logetavector, df, timeknots, i;
    extrapl=zeros(1), extrapr=zeros(1), kwargs...
)
    _logeffectivereproductionratio_etavector!(logetavector, df, i; kwargs...)
    return CubicSpline(timeknots, logetavector; extrapl, extrapr)
end

function _logeffectivereproductionratio_etavector!(
    logetavector, df, i;
    kwargs...
)
    for k ∈ eachindex(logetavector)
        _logeffectivereproductionratio_etavalue!(logetavector, df, i, k; kwargs...)
    end
end

function _logeffectivereproductionratio_etavalue!(logetavector, df, i, k; kwargs...)
    logetavector[k] = get(kwargs, Symbol("eta_t$k")) do 
        getproperty(df, Symbol("eta_t$k.logeta"))[i]
    end
end

function _logeffectivereproductionratio_zetavalue(df, i, j; kwargs...)
    get(kwargs, Symbol("logzeta_g$j")) do 
        getproperty(df, Symbol("logzeta_g$j.logzeta"))[i]
    end
end

function _logeffectivereproductionratio_secondarydeltavector!(
    ::Nothing, ::Any, ::Any; 
    kwargs...
)
    nothing
end

function _logeffectivereproductionratio_secondarydeltavector!(
    logsecondarydeltavector, df, i;
    kwargs...
)
    for k ∈ eachindex(logsecondarydeltavector)
        _logeffectivereproductionratio_secondarydeltavalue!(
            logsecondarydeltavector, df, i, k; 
            kwargs...
        )
    end
end

function _logeffectivereproductionratio_secondarydeltavalue!(
    logsecondarydeltavector, df, i, k; 
    kwargs...
)
    logsecondarydeltavector[k] = get(kwargs, Symbol("logsecondarydelta$k")) do 
        getproperty(df, Symbol("logsecondarydelta$k.logsecondarydelta"))[i]
    end
end
=#
function _logeffectivereproductionratio_value(df, symbol, i; kwargs...)
    get(kwargs, symbol) do 
        getproperty(df, symbol)[i]
    end
end

function predictcases(g::Function, df, logR0, initialcases, ns; kwargs...)
    ninitialcases = size(initialcases, 1)
    predictedcases = similar(logR0)
    Threads.@threads for k ∈ axes(predictedcases, 1)
        identifiedproportion = _logeffectivereproductionratio_value(df, :psi, k; kwargs...)
        for t ∈ axes(predictedcases, 2)
            if t <= ninitialcases
                predictedcases[k, t, :] = initialcases[t, :]
            else
                for j ∈ axes(predictedcases, 3)
                    predictedcases[k, t, j] = (
                        exp(logR0[k, t, j]) * 
                        (
                            1 - 
                            sum(@view predictedcases[k, 1:(t-1), j]) / 
                            (identifiedproportion * ns[j])
                        ) * 
                        sum([ predictedcases[k, x, j] * g(t - x) for x ∈ 1:(t-1) ]) 
                    )
                end
            end
        end
    end
    return predictedcases
end

function predictcases(g::Vector, df, logR0, initialcases, ns; kwargs...)
    ℓ = length(g) 
    return predictcases(x -> x > ℓ ? 0 : g[x], df, logR0, initialcases, ns; kwargs...)
end
