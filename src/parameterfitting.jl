#=
function generatew_gt(f, cases, psi, Ns; blankn=0)
    testf(f)
    w_gt = Matrix{Float64}(undef, size(cases))
    for g ∈ axes(w_gt, 2)
        for t ∈ axes(w_gt, 1)
            if t <= blankn 
                w_gt[t, g] = NaN
            elseif t == 1 
                w_gt[t, g] = _generatew_gtrow_cases(cases, t, g)
            else
                w_gt[t, g] = _generatew_gtrow_cases(cases, t, g) +
                    #_generatew_gtrow_susceptibles(cases, psi, Ns, t, g) +
                    _generatew_gtrow_gentime(f, cases, t, g)
            end
        end
    end
    return w_gt
end
=#

function generatew_gt(f, cases, Ns; blankn=0)
    testf(f)
    w_gt = Matrix{Float64}(undef, size(cases))
    for g ∈ axes(w_gt, 2)
        for t ∈ axes(w_gt, 1)
            if t <= blankn 
                w_gt[t, g] = NaN
            elseif t == 1 
                w_gt[t, g] = _generatew_gtrow_cases(cases, t, g)
            else
                w_gt[t, g] = _generatew_gtrow_cases(cases, t, g) +
                    _generatew_gtrow_gentime(f, cases, t, g)
            end
        end
    end
    return w_gt
end

_generatew_gtrow_cases(cases, t, g) = log(cases[t, g])

#function _generatew_gtrow_susceptibles(cases, psi, Ns, t, g)
#    -log(1 - (sum(cases[1:(t - 1), g])) / (psi * Ns[g]))
#end

function _generatew_gtrow_gentime(f::Function, cases, t, g)
    return -log(sum([ f(τ) * cases[t-τ, g] for τ ∈ 1:(t - 1) ]))
end

function _generatew_gtrow_gentime(f::AbstractVector, cases, t, g)
    taumax = min(t - 1, length(f))
    return -log(sum([ f[τ] * cases[t-τ, g] for τ ∈ 1:taumax ]))
end

@model _estimatelogzeta_g(μ, σ2) = logzeta ~ Normal(μ, sqrt(σ2))
@model _estimatelogeta_t(μ, σ2) = logeta ~ Normal(μ, sqrt(σ2))
@model _estimatelogsecondarydelta(deltaprior) = logsecondarydelta ~ deltaprior

function diffindiffparameters_discretetimes(; w, y, interventions, timeperiods, Ns, kwargs...)
    return diffindiffparameters_discretetimes(w, y, interventions, timeperiods, Ns; kwargs...)
end

@model function diffindiffparameters_discretetimes(
    w, y, interventions, timeperiods, Ns;
    etavarianceprior=Exponential(1),
    logdeltaprior=Normal(0, 1),
    logmeanetaprior=Normal(0, 1),
    logmeanzetaprior=Normal(0, 1),
    peakperiod=1,
    psiprior=0.5,
    secondaryinterventions=nothing,
    sigmasquareprior=Exponential(1),
    zetavarianceprior=Exponential(1),
)
    ntimeperiods = _nunique(timeperiods)
    nlocations = size(w, 2)

    logzeta_mean ~ logmeanzetaprior
    zeta_sigma2 ~ zetavarianceprior
    logzeta_g_vec = Vector{typeof(logzeta_mean)}(undef, nlocations) 
    for i ∈ 1:nlocations
        @submodel prefix="logzeta_g$i" logzeta = _estimatelogzeta_g(
            logzeta_mean, zeta_sigma2
        )
        logzeta_g_vec[i] = logzeta
    end

    logeta_mean ~ logmeanetaprior
    eta_sigma2 ~ etavarianceprior
    logeta_t_vec = Vector{typeof(logeta_mean)}(undef, ntimeperiods) 
    for i ∈ 1:ntimeperiods
        if i == peakperiod
            logeta = zero(typeof(logeta_mean))
        else
            @submodel prefix="eta_t$i" logeta = _estimatelogeta_t(logeta_mean, eta_sigma2)
        end
        logeta_t_vec[i] = logeta
    end

    logdelta ~ logdeltaprior
    σ2 ~ sigmasquareprior 

    psiminimum = maximum([ sum(y[:, g]) / Ns[g] for g ∈ axes(y, 2) ])

    if isa(psiprior, Distribution)
        psi ~ truncated(psiprior, psiminimum, 1) 
    else
        @assert 0 < psiprior <= 1
        a = 10 * psiprior 
        b = 10 - a 
        psi ~ truncated(Beta(a, b), psiminimum, 1) 
    end

    if isnothing(secondaryinterventions)
        modeloutput1 = logzeta_g_vec[1] + logeta_t_vec[1] + logdelta * interventions[1, 1] 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = logzeta_g_vec[g] + logeta_t_vec[timeperiods[t]] + logdelta * interventions[t, g] + s
        end
    elseif isa(secondaryinterventions, AbstractMatrix)
        logsecondarydelta ~ logdeltaprior
        modeloutput1 = logzeta_g_vec[1] + logeta_t_vec[1] + logdelta * interventions[1, 1] + logsecondarydelta * secondaryinterventions[1, 1]
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = logzeta_g_vec[g] + logeta_t_vec[timeperiods[t]] + logdelta * interventions[t, g] + s + logsecondarydelta * secondaryinterventions[t, g] 
        end
    else
        logsecondarydelta_vec = Vector{typeof(logdelta)}(
            undef, length(secondaryinterventions)
        ) 
        for i ∈ eachindex(secondaryinterventions)
            @submodel prefix="logsecondarydelta$i" logsecondarydelta = 
                _estimatelogsecondarydelta(logdeltaprior)
                logsecondarydelta_vec[i] = logsecondarydelta
        end

        modeloutput1 = logzeta_g_vec[1] + logeta_t_vec[1] + logdelta * interventions[1, 1] + + sum([ logsecondarydelta_vec[i] * secondaryinterventions[i][1, 1] for i ∈ eachindex(secondaryinterventions) ])
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = logzeta_g_vec[g] + logeta_t_vec[timeperiods[t]] + logdelta * interventions[t, g] + s + sum([ logsecondarydelta_vec[i] * secondaryinterventions[i][t, g] for i ∈ eachindex(secondaryinterventions) ]) 
        end
    end

    for i ∈ eachindex(w)
        _skip(w[i]) && continue 
        w[i] ~ Normal(modeloutput[i], σ2)
    end 
end

function diffindiffparameters_splinetimes(; w, y, interventions, timeknots, Ns, kwargs...)
    return diffindiffparameters_splinetimes(w, y, interventions, timeknots, Ns; kwargs...)
end

@model function diffindiffparameters_splinetimes(
    w, y, interventions, timeknots, Ns;
    etavarianceprior=Exponential(1),
    logdeltaprior=Normal(0, 1),
    logmeanetaprior=Normal(0, 1),
    logmeanzetaprior=Normal(0, 1),
    peakperiod=2,
    psiprior=0.5,
    secondaryinterventions=nothing,
    sigmasquareprior=Exponential(1),
    zetavarianceprior=Exponential(1),
    extrapl=zeros(1), extrapr=zeros(1), 
)
    ntimeknots = length(timeknots)
    nlocations = size(w, 2)

    logzeta_mean ~ logmeanzetaprior
    zeta_sigma2 ~ zetavarianceprior
    logzeta_g_vec = Vector{typeof(logzeta_mean)}(undef, nlocations) 
    for i ∈ 1:nlocations
        @submodel prefix="logzeta_g$i" logzeta = _estimatelogzeta_g(
            logzeta_mean, zeta_sigma2
        )
        logzeta_g_vec[i] = logzeta
    end

    logeta_mean ~ logmeanetaprior
    eta_sigma2 ~ etavarianceprior
    logeta_t_vec = Vector{typeof(logeta_mean)}(undef, ntimeknots) 
    for i ∈ 1:ntimeknots
        if i == peakperiod
            logeta = zero(typeof(logeta_mean))
        else
            @submodel prefix="eta_t$i" logeta = _estimatelogeta_t(logeta_mean, eta_sigma2)
        end
        logeta_t_vec[i] = logeta
    end
    timespline = CubicSpline(timeknots, ForwardDiff.value.(logeta_t_vec); extrapl, extrapr)

    logdelta ~ logdeltaprior
    σ2 ~ sigmasquareprior 

    psiminimum = maximum([sum(y[:, g]) / Ns[g] for g ∈ axes(w, 2)])

    if isa(psiprior, Distribution)
        psi ~ truncated(psiprior, psiminimum, 1) 
    else
        @assert 0 < psiprior <= 1
        a = 10 * psiprior 
        b = 10 - a 
        psi ~ truncated(Beta(a, b), psiminimum, 1) 
    end

    if isnothing(secondaryinterventions)
        modeloutput1 = logzeta_g_vec[1] + timespline[1] + logdelta * interventions[1, 1] 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = logzeta_g_vec[g] + timespline[t] + logdelta * interventions[t, g] + s
        end
    elseif isa(secondaryinterventions, AbstractMatrix)
        logsecondarydelta ~ logdeltaprior
        modeloutput1 = logzeta_g_vec[1] + timespline[1] + logdelta * interventions[1, 1] + logsecondarydelta * secondaryinterventions[1, 1]
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = logzeta_g_vec[g] + timespline[t] + logdelta * interventions[t, g] + s + logsecondarydelta * secondaryinterventions[t, g] 
        end
    else
        logsecondarydelta_vec = Vector{typeof(logdelta)}(
            undef, length(secondaryinterventions)
        ) 
        for i ∈ eachindex(secondaryinterventions)
            @submodel prefix="logsecondarydelta$i" logsecondarydelta = 
                _estimatelogsecondarydelta(logdeltaprior)
                logsecondarydelta_vec[i] = logsecondarydelta
        end

        modeloutput1 = logzeta_g_vec[1] + timespline[1] + logdelta * interventions[1, 1] + + sum([ logsecondarydelta_vec[i] * secondaryinterventions[i][1, 1] for i ∈ eachindex(secondaryinterventions) ])
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = logzeta_g_vec[g] + timespline[t] + logdelta * interventions[t, g] + s + sum([ logsecondarydelta_vec[i] * secondaryinterventions[i][t, g] for i ∈ eachindex(secondaryinterventions) ]) 
        end
    end

    for i ∈ eachindex(w)
        _skip(w[i]) && continue 
        w[i] ~ Normal(modeloutput[i], σ2)
    end 
end

function diffindiffparameters_polytimes(; w, y, interventions, timeknots, Ns, kwargs...)
    return diffindiffparameters_polytimes(w, y, interventions, timeknots, Ns; kwargs...)
end

@model function diffindiffparameters_polytimes(
    w, y, interventions, timeknots, Ns;
    etapowerprior=Uniform(-3, 3),
    etasprior=Exponential(1),
    etavarianceprior=Exponential(1),
    logdeltaprior=Normal(0, 1),
    logmeanetaprior=Normal(0, 1),
    logmeanzetaprior=Normal(0, 1),
    peakperiod=2,
    psiprior=0.5,
    secondaryinterventions=nothing,
    sigmasquareprior=Exponential(1),
    zetavarianceprior=Exponential(1),
    extrapl=zeros(1), extrapr=zeros(1), 
)
    ntimeknots = length(timeknots)
    nlocations = size(w, 2)

    logzeta_mean ~ logmeanzetaprior
    zeta_sigma2 ~ zetavarianceprior
    logzeta_g_vec = Vector{typeof(logzeta_mean)}(undef, nlocations) 
    for i ∈ 1:nlocations
        @submodel prefix="logzeta_g$i" logzeta = _estimatelogzeta_g(
            logzeta_mean, zeta_sigma2
        )
        logzeta_g_vec[i] = logzeta
    end

    logeta_mean ~ logmeanetaprior
    eta_sigma2 ~ etavarianceprior
    logeta_t_vec = Vector{typeof(logeta_mean)}(undef, 4) 
    for i ∈ 1:4
        if i == peakperiod
            logeta = zero(typeof(logeta_mean))
        else
            @submodel prefix="eta_t$i" logeta = _estimatelogeta_t(logeta_mean, eta_sigma2)
        end
        logeta_t_vec[i] = logeta
    end

    etapower1 ~ Uniform(-3, -0.5)
    etapower2 ~ Uniform(etapower1, 0.5)
    etapower3 ~ Uniform(etapower2, 3)
    etapowers = [ etapower1, etapower2, etapower3 ]

    eta(t) = sum([ ltv * (t)^etp for (ltv, etp) ∈ zip(logeta_t_vec[1:3], etapowers) ]) + logeta_t_vec[4] * log(t)
    
    logdelta ~ logdeltaprior
    σ2 ~ sigmasquareprior 

    psiminimum = maximum([sum(y[:, g]) / Ns[g] for g ∈ axes(w, 2)])

    if isa(psiprior, Distribution)
        psi ~ truncated(psiprior, psiminimum, 1) 
    else
        @assert 0 < psiprior <= 1
        a = 10 * psiprior 
        b = 10 - a 
        psi ~ truncated(Beta(a, b), psiminimum, 1) 
    end

    if isnothing(secondaryinterventions)
        modeloutput1 = logzeta_g_vec[1] + eta(1) + logdelta * interventions[1, 1] 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = logzeta_g_vec[g] + eta(t) + logdelta * interventions[t, g] + s
        end
    elseif isa(secondaryinterventions, AbstractMatrix)
        logsecondarydelta ~ logdeltaprior
        modeloutput1 = logzeta_g_vec[1] + eta(1) + logdelta * interventions[1, 1] + logsecondarydelta * secondaryinterventions[1, 1]
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = logzeta_g_vec[g] + eta(t) + logdelta * interventions[t, g] + s + logsecondarydelta * secondaryinterventions[t, g] 
        end
    else
        logsecondarydelta_vec = Vector{typeof(logdelta)}(
            undef, length(secondaryinterventions)
        ) 
        for i ∈ eachindex(secondaryinterventions)
            @submodel prefix="logsecondarydelta$i" logsecondarydelta = 
                _estimatelogsecondarydelta(logdeltaprior)
                logsecondarydelta_vec[i] = logsecondarydelta
        end

        modeloutput1 = logzeta_g_vec[1] + eta(1) + logdelta * interventions[1, 1] + + sum([ logsecondarydelta_vec[i] * secondaryinterventions[i][1, 1] for i ∈ eachindex(secondaryinterventions) ])
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = logzeta_g_vec[g] + eta(t) + logdelta * interventions[t, g] + s + sum([ logsecondarydelta_vec[i] * secondaryinterventions[i][t, g] for i ∈ eachindex(secondaryinterventions) ]) 
        end
    end

    for i ∈ eachindex(w)
        _skip(w[i]) && continue 
        w[i] ~ Normal(modeloutput[i], σ2)
    end 
end

#=
function diffindiffparameters(
    w, interventions; 
    timeknots=nothing, timeperiods=nothing, kwargs...
)
    return diffindiffparameters(w, interventions, timeknots, timeperiods; kwargs...)
end

function diffindiffparameters(w, interventions, ::Nothing, timeperiods; kwargs...)
    return diffindiffparameters_discretetimes(w, interventions, timeperiods; kwargs...)
end

function diffindiffparameters(w, interventions, timeknots, ::Nothing; kwargs...)
    return diffindiffparameters_splinetimes(w, interventions, timeknots; kwargs...)
end
=#

_skip(x) = isnan(x) || x == -Inf || x == Inf 
