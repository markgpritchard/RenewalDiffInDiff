
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
                w_gt[t, g] = +(
                    _generatew_gtrow_cases(cases, t, g),
                    _generatew_gtrow_gentime(f, cases, t, g)
                ) 
            end
        end
    end
    return w_gt
end

_generatew_gtrow_cases(cases, t, g) = log(cases[t, g])

function _generatew_gtrow_gentime(f::Function, cases, t, g)
    return -log(sum([ f(τ) * cases[t-τ, g] for τ ∈ 1:(t - 1) ]))
end

function _generatew_gtrow_gentime(f::AbstractVector, cases, t, g)
    taumax = min(t - 1, length(f))
    return -log(sum([ f[τ] * cases[t-τ, g] for τ ∈ 1:taumax ]))
end

@model _estimatelogzeta_g(μ, σ2) = logzeta ~ Normal(μ, sqrt(σ2))
@model _estimatelogeta_t(μ, σ2) = logeta ~ Normal(μ, sqrt(σ2))
@model _estimatelogeta_t(etaprior) = logeta ~ etaprior
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
            modeloutput[t, g] = +(
                logzeta_g_vec[g], 
                logeta_t_vec[timeperiods[t]], 
                logdelta * interventions[t, g], 
                s
            ) 
        end
    elseif isa(secondaryinterventions, AbstractMatrix)
        logsecondarydelta ~ logdeltaprior
        modeloutput1 = +(
            logzeta_g_vec[1],
            logeta_t_vec[1], 
            logdelta * interventions[1, 1], 
            logsecondarydelta * secondaryinterventions[1, 1]
        ) 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = +(
                logzeta_g_vec[g],
                logeta_t_vec[timeperiods[t]],
                logdelta * interventions[t, g],
                s,
                logsecondarydelta * secondaryinterventions[t, g] 
            ) 
        end
    else
        logsecondarydelta_vec = Vector{typeof(logdelta)}(
            undef, length(secondaryinterventions)
        ) 
        for i ∈ eachindex(secondaryinterventions)
            @submodel prefix="logsecondarydelta$i" logsecondarydelta = _estimatelogsecondarydelta(
                logdeltaprior
            )
            logsecondarydelta_vec[i] = logsecondarydelta
        end

        modeloutput1 = +(
            logzeta_g_vec[1],
            logeta_t_vec[1], 
            logdelta * interventions[1, 1],
            sum(
                [ 
                    logsecondarydelta_vec[i] * secondaryinterventions[i][1, 1] 
                    for i ∈ eachindex(secondaryinterventions) 
                ]
            )
        ) 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = +(
                logzeta_g_vec[g], 
                logeta_t_vec[timeperiods[t]],
                logdelta * interventions[t, g],
                s,
                sum(
                    [ 
                        logsecondarydelta_vec[i] * secondaryinterventions[i][t, g] 
                        for i ∈ eachindex(secondaryinterventions) 
                    ]
                )
            )  
        end
    end

    for i ∈ eachindex(w)
        _skip(w[i]) && continue 
        w[i] ~ Normal(modeloutput[i], σ2)
    end 
end

function diffindiffparameters_twodiscretetimes(; w, y, interventions, timeperiods, Ns, kwargs...)
    return diffindiffparameters_twodiscretetimes(w, y, interventions, timeperiods, Ns; kwargs...)
end

@model function diffindiffparameters_twodiscretetimes(
    w, y, interventions, timeperiods, Ns;
    logdeltaprior=Normal(0, 1),
    logeta2eprior=Normal(0, 1),
    logmeanzetaprior=Normal(0, 1),
    peakperiod=1,
    psiprior=0.5,
    secondaryinterventions=nothing,
    sigmasquareprior=Exponential(1),
    zetavarianceprior=Exponential(1),
)
    ntimeperiods = _nunique(timeperiods)
    @assert ntimeperiods == 2
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

    logeta_t_vec = Vector{typeof(logzeta_mean)}(undef, ntimeperiods) 
    for i ∈ 1:ntimeperiods
        if i == peakperiod
            logeta = zero(typeof(logzeta_mean))
        else
            @submodel prefix="eta_t$i" logeta = _estimatelogeta_t(logeta2eprior)
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
            modeloutput[t, g] = +(
                logzeta_g_vec[g], 
                logeta_t_vec[timeperiods[t]], 
                logdelta * interventions[t, g], 
                s
            ) 
        end
    elseif isa(secondaryinterventions, AbstractMatrix)
        logsecondarydelta ~ logdeltaprior
        modeloutput1 = +(
            logzeta_g_vec[1],
            logeta_t_vec[1], 
            logdelta * interventions[1, 1], 
            logsecondarydelta * secondaryinterventions[1, 1]
        ) 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = +(
                logzeta_g_vec[g],
                logeta_t_vec[timeperiods[t]],
                logdelta * interventions[t, g],
                s,
                logsecondarydelta * secondaryinterventions[t, g] 
            ) 
        end
    else
        logsecondarydelta_vec = Vector{typeof(logdelta)}(
            undef, length(secondaryinterventions)
        ) 
        for i ∈ eachindex(secondaryinterventions)
            @submodel prefix="logsecondarydelta$i" logsecondarydelta = _estimatelogsecondarydelta(
                logdeltaprior
            )
            logsecondarydelta_vec[i] = logsecondarydelta
        end

        modeloutput1 = +(
            logzeta_g_vec[1],
            logeta_t_vec[1], 
            logdelta * interventions[1, 1],
            sum(
                [ 
                    logsecondarydelta_vec[i] * secondaryinterventions[i][1, 1] 
                    for i ∈ eachindex(secondaryinterventions) 
                ]
            )
        ) 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = +(
                logzeta_g_vec[g], 
                logeta_t_vec[timeperiods[t]],
                logdelta * interventions[t, g],
                s,
                sum(
                    [ 
                        logsecondarydelta_vec[i] * secondaryinterventions[i][t, g] 
                        for i ∈ eachindex(secondaryinterventions) 
                    ]
                )
            )  
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

    psiminimum = maximum([ sum(y[:, g]) / Ns[g] for g ∈ axes(w, 2) ])

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
            modeloutput[t, g] = +(
                logzeta_g_vec[g], 
                timespline[t],
                logdelta * interventions[t, g],
                s
            ) 
        end
    elseif isa(secondaryinterventions, AbstractMatrix)
        logsecondarydelta ~ logdeltaprior
        modeloutput1 = +(
            logzeta_g_vec[1], 
            timespline[1],
            logdelta * interventions[1, 1],
            logsecondarydelta * secondaryinterventions[1, 1]
        ) 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = +(
                logzeta_g_vec[g],
                timespline[t],
                logdelta * interventions[t, g],
                s,
                logsecondarydelta * secondaryinterventions[t, g] 
            )
        end
    else
        logsecondarydelta_vec = Vector{typeof(logdelta)}(
            undef, length(secondaryinterventions)
        ) 
        for i ∈ eachindex(secondaryinterventions)
            @submodel prefix="logsecondarydelta$i" logsecondarydelta = _estimatelogsecondarydelta(
                logdeltaprior
            )
            logsecondarydelta_vec[i] = logsecondarydelta
        end

        modeloutput1 = +(
            logzeta_g_vec[1], 
            timespline[1],
            logdelta * interventions[1, 1],
            sum(
                [ 
                    logsecondarydelta_vec[i] * secondaryinterventions[i][1, 1] 
                    for i ∈ eachindex(secondaryinterventions) 
                ]
            )
        ) 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = +(
                logzeta_g_vec[g],
                timespline[t],
                logdelta * interventions[t, g],
                s,
                sum(
                    [ 
                        logsecondarydelta_vec[i] * secondaryinterventions[i][t, g] 
                        for i ∈ eachindex(secondaryinterventions) 
                    ]
                )
            )  
        end
    end

    for i ∈ eachindex(w)
        _skip(w[i]) && continue 
        w[i] ~ Normal(modeloutput[i], σ2)
    end 
end

#=
function diffindiffparameters_splinetimes_twopsis( ; 
    w, y, interventions, timeknots, Ns, psiratio, kwargs...
)
    return diffindiffparameters_splinetimes_twopsis(
        w, y, interventions, timeknots, Ns, psiratio; 
        kwargs...
    )
end

@model function diffindiffparameters_splinetimes_twopsis(
    w, y, interventions, timeknots, Ns, psiratio;
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

    @assert psiratio >= 1
    psiminimum = maximum([ sum(y[:, g]) / Ns[g] for g ∈ axes(w, 2) ])
    psimaximum = max(1 / psiratio, psiminimum)

    if isa(psiprior, Distribution)
        psi ~ truncated(psiprior, psiminimum, psimaximum) 
    else
        @assert 0 < psiprior <= 1
        a = 10 * psiprior 
        b = 10 - a 
        psi ~ truncated(Beta(a, b), psiminimum, psimaximum) 
    end
    inflatedpsi = min(1, psi * psiratio)
    newmultiplier = inflatedpsi / psi

    if isnothing(secondaryinterventions)
        modeloutput1 = logzeta_g_vec[1] + timespline[1] + logdelta * interventions[1, 1] 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            elseif interventions[t, g] == 0
                s = log(1 - sum(@view y[1:(t - 1), g]) / (Ns[g] * psi))
            else
                A = interventions.starttimes[g]
                s = log(
                    +(
                        1,
                        -sum(@view y[1:(A - 1), g]) / (Ns[g] * psi),
                        -sum(@view y[A:(t - 1), g]) / (Ns[g] * inflatedpsi)
                    ) 
                )
            end
            modeloutput[t, g] = +(
                logzeta_g_vec[g], 
                timespline[t],
                logdelta * interventions[t, g],
                s
            ) 
        end
    elseif isa(secondaryinterventions, AbstractMatrix)
        logsecondarydelta ~ logdeltaprior
        modeloutput1 = +(
            logzeta_g_vec[1], 
            timespline[1],
            logdelta * interventions[1, 1],
            logsecondarydelta * secondaryinterventions[1, 1]
        ) 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            elseif interventions[t, g] == 0
                s = log(1 - sum(@view y[1:(t - 1), g]) / (Ns[g] * psi))
            else
                A = interventions.starttimes[g]
                s = log(
                    +(
                        1,
                        -sum(@view y[1:(A - 1), g]) / (Ns[g] * psi),
                        -sum(@view y[A:(t - 1), g]) / (Ns[g] * inflatedpsi)
                    ) 
                )
            end
            modeloutput[t, g] = +(
                logzeta_g_vec[g],
                timespline[t],
                logdelta * interventions[t, g],
                s,
                logsecondarydelta * secondaryinterventions[t, g] 
            )
        end
    else
        logsecondarydelta_vec = Vector{typeof(logdelta)}(
            undef, length(secondaryinterventions)
        ) 
        for i ∈ eachindex(secondaryinterventions)
            @submodel prefix="logsecondarydelta$i" logsecondarydelta = _estimatelogsecondarydelta(
                logdeltaprior
            )
            logsecondarydelta_vec[i] = logsecondarydelta
        end

        modeloutput1 = +(
            logzeta_g_vec[1], 
            timespline[1],
            logdelta * interventions[1, 1],
            sum(
                [ 
                    logsecondarydelta_vec[i] * secondaryinterventions[i][1, 1] 
                    for i ∈ eachindex(secondaryinterventions) 
                ]
            )
        ) 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            elseif interventions[t, g] == 0
                s = log(1 - sum(@view y[1:(t - 1), g]) / (Ns[g] * psi))
            else
                A = interventions.starttimes[g]
                s = log(
                    +(
                        1,
                        -sum(@view y[1:(A - 1), g]) / (Ns[g] * psi),
                        -sum(@view y[A:(t - 1), g]) / (Ns[g] * inflatedpsi)
                    ) 
                )
            end
            modeloutput[t, g] = +(
                logzeta_g_vec[g],
                timespline[t],
                logdelta * interventions[t, g],
                s,
                sum(
                    [ 
                        logsecondarydelta_vec[i] * secondaryinterventions[i][t, g] 
                        for i ∈ eachindex(secondaryinterventions) 
                    ]
                )
            )  
        end
    end

    for i ∈ eachindex(w)
        _skip(w[i]) && continue 
        w[i] ~ Normal(modeloutput[i], σ2)
    end 
end
=#

@model function diffindiffparameters_fittocurve_splinetimes(
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
    nlocations = size(y, 2)

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
        modeloutput1 = logzeta_g_vec[1] + timespline[1] + logdelta * interventions[1, 1] 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(y))
        for t ∈ axes(y, 1), g ∈ axes(y, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = +(
                logzeta_g_vec[g], 
                timespline[t],
                logdelta * interventions[t, g],
                s
            ) 
        end
    elseif isa(secondaryinterventions, AbstractMatrix)
        logsecondarydelta ~ logdeltaprior
        modeloutput1 = +(
            logzeta_g_vec[1], 
            timespline[1],
            logdelta * interventions[1, 1],
            logsecondarydelta * secondaryinterventions[1, 1]
        ) 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(y))
        for t ∈ axes(y, 1), g ∈ axes(y, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = +(
                logzeta_g_vec[g],
                timespline[t],
                logdelta * interventions[t, g],
                s,
                logsecondarydelta * secondaryinterventions[t, g] 
            )
        end
    else
        logsecondarydelta_vec = Vector{typeof(logdelta)}(
            undef, length(secondaryinterventions)
        ) 
        for i ∈ eachindex(secondaryinterventions)
            @submodel prefix="logsecondarydelta$i" logsecondarydelta = _estimatelogsecondarydelta(
                logdeltaprior
            )
            logsecondarydelta_vec[i] = logsecondarydelta
        end

        modeloutput1 = +(
            logzeta_g_vec[1], 
            timespline[1],
            logdelta * interventions[1, 1],
            sum(
                [ 
                    logsecondarydelta_vec[i] * secondaryinterventions[i][1, 1] 
                    for i ∈ eachindex(secondaryinterventions) 
                ]
            )
        ) 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(y))
        for t ∈ axes(y, 1), g ∈ axes(y, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
                modeloutput[t, g] = zy[t, g]
            elseif t <= 20 
                s = log(1 - sum(modeloutput[1:(t - 1), g]) / (Ns[g] * psi))
                modeloutput[t, g] = zy[t, g]
            else 
                s = log(1 - sum(modeloutput[1:(t - 1), g]) / (Ns[g] * psi))
            
                modeloutput[t, g] = *(
                    exp(logzeta_g_vec[g]),
                    exp(timespline[t]),
                    exp(logdelta)^interventions[t, g],
                    s,
                    sum(
                        [ 
                            logsecondarydelta_vec[i] * secondaryinterventions[i][t, g] 
                            for i ∈ eachindex(secondaryinterventions) 
                        ]
                    )
                ) 
            end
        end
    end

    for i ∈ eachindex(w)
        #_skip(w[i]) && continue 
        y[i] ~ Normal(modeloutput[i], σ2)
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

    eta(t) = +(
        sum(
            [ 
                ltv * (t)^etp 
                for (ltv, etp) ∈ zip(logeta_t_vec[1:3], etapowers) 
            ]
        ),
        logeta_t_vec[4] * log(t)
    ) 
    
    logdelta ~ logdeltaprior
    σ2 ~ sigmasquareprior 

    psiminimum = maximum([ sum(y[:, g]) / Ns[g] for g ∈ axes(w, 2) ])

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
            modeloutput[t, g] = +(
                logzeta_g_vec[g], eta(t), logdelta * interventions[t, g], s
            ) 
        end
    elseif isa(secondaryinterventions, AbstractMatrix)
        logsecondarydelta ~ logdeltaprior
        modeloutput1 = +(
            logzeta_g_vec[1],
            eta(1),
            logdelta * interventions[1, 1],
            logsecondarydelta * secondaryinterventions[1, 1]
        ) 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = +(
                logzeta_g_vec[g],
                eta(t),
                logdelta * interventions[t, g],
                s,
                logsecondarydelta * secondaryinterventions[t, g]
            )  
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

        modeloutput1 = +(
            logzeta_g_vec[1],
            eta(1),
            logdelta * interventions[1, 1],
            sum(
                [ 
                    logsecondarydelta_vec[i] * secondaryinterventions[i][1, 1] 
                    for i ∈ eachindex(secondaryinterventions) 
                ]
            )
        ) 
        modeloutput = Matrix{typeof(modeloutput1)}(undef, size(w))
        for t ∈ axes(w, 1), g ∈ axes(w, 2) 
            if t == 1 
                s = zero(log(1 / (Ns[g] * psi))) 
            else 
                s = log(1 - sum(y[1:(t - 1), g]) / (Ns[g] * psi))
            end
            modeloutput[t, g] = +(
                logzeta_g_vec[g],
                eta(t),
                logdelta * interventions[t, g],
                s,
                sum(
                    [ 
                        logsecondarydelta_vec[i] * secondaryinterventions[i][t, g] 
                        for i ∈ eachindex(secondaryinterventions) 
                    ]
                ) 
            ) 
        end
    end

    for i ∈ eachindex(w)
        _skip(w[i]) && continue 
        w[i] ~ Normal(modeloutput[i], σ2)
    end 
end

_skip(x) = isnan(x) || x == -Inf || x == Inf 

function loadanalysisdictsasdf(modelname, n_chains, maxrounds, seedstart; nsims=4)
    fn1 = _findanalysisfilename(modelname, n_chains, maxrounds, seedstart + 1)
    chain1 = load(fn1)["chain"]
    df = DataFrame(chain1)
    for i ∈ 2:nsims 
        fn = _findanalysisfilename(modelname, n_chains, maxrounds, seedstart + i)
        isnothing(fn) && continue
        chain = load(fn)["chain"]
        _tdf = DataFrame(chain)
        _tdf.chain = [ i for _ ∈ axes(_tdf, 1) ]
        df = vcat(df, _tdf) 
    end
    return df
end

function _findanalysisfilename(modelname, n_chains, maxrounds, seed)
    n_rounds = maxrounds
    while true 
        filename = "modelname=$(modelname)_n_chains=$(n_chains)_n_rounds=$(n_rounds)_seed=$(seed).jld2"
        if isfile(datadir("sims", filename))
            return datadir("sims", filename)
        else 
            n_rounds += -1 
            if n_rounds <= 0 return nothing end
        end
    end
end

function keyvalues(fitteddf, fittedvaluesset; cri=[ 0.05, 0.95 ])
    deltamean = exp(mean(fitteddf.logdelta))
    deltap05_95 = exp.(quantile(fitteddf.logdelta, cri))
    totalcases = Matrix{Int}(
        undef, 
        length(fittedvaluesset.y_matrix_poisson_vec), 
        size(fittedvaluesset.y_matrix_poisson_vec[1], 2)
    )
    peakcases = Matrix{Int}(
        undef, 
        length(fittedvaluesset.y_matrix_poisson_vec),
        size(fittedvaluesset.y_matrix_poisson_vec[1], 2)
    )
    peakcasesdate = Matrix{Int}(
        undef, 
        length(fittedvaluesset.y_matrix_poisson_vec), 
        size(fittedvaluesset.y_matrix_poisson_vec[1], 2)
    )
    for i ∈ eachindex(fittedvaluesset.y_matrix_poisson_vec)
        totalcases[i, :] = [
            sum(fittedvaluesset.y_matrix_poisson_vec[i][:, j]) - 
                sum(fittedvaluesset.y_matrix_poisson_vec_counterfactual[i][:, j])
            for j ∈ axes(fittedvaluesset.y_matrix_poisson_vec[i], 2)
        ]
        for j ∈ axes(fittedvaluesset.y_matrix_poisson_vec[i], 2)
            x, ind = findmax(fittedvaluesset.y_matrix_poisson_vec[i][:, j])
            x_cf, ind_cf = findmax(
                fittedvaluesset.y_matrix_poisson_vec_counterfactual[i][:, j]
            )
            peakcases[i, j] = x - x_cf 
            peakcasesdate[i, j] = ind - ind_cf
        end
    end
    totalcasesmeaneffect = [ mean(totalcases[:, i]) for i ∈ axes(totalcases, 2) ]
    totalcasesp05_90 = [ quantile(totalcases[:, i], cri) for i ∈ axes(totalcases, 2) ]
    peakcasesmeaneffect = [ mean(peakcases[:, i]) for i ∈ axes(peakcases, 2) ]
    peakcasesp05_90 = [ quantile(peakcases[:, i], cri) for i ∈ axes(peakcases, 2) ]
    peakcasesdatemeaneffect = [ mean(peakcasesdate[:, i]) for i ∈ axes(peakcasesdate, 2) ]
    peakcasesdatep05_90 = [ quantile(peakcasesdate[:, i], cri) for i ∈ axes(peakcasesdate, 2) ]
    return (
        deltamean=deltamean,
        deltap05_95=deltap05_95,
        #totalcases=totalcases,
        #peakcases=peakcases,
        #peakcasesdate=peakcasesdate,
        totalcasesmeaneffect=totalcasesmeaneffect,
        totalcasesp05_90=totalcasesp05_90,
        peakcasesmeaneffect=peakcasesmeaneffect,
        peakcasesp05_90=peakcasesp05_90,
        peakcasesdatemeaneffect=peakcasesdatemeaneffect,
        peakcasesdatep05_90=peakcasesdatep05_90,
    )
end
