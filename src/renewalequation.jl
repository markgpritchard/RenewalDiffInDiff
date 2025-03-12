



function renewalequation_expectedcases(ρ::Number, s::Number, f_output::Number; kwargs...)
    return *(ρ, s, f_output)
end
 
function renewalequation_expectedcases(
    f, y_vector::AbstractVector, ρ::Number, N::Number; 
    s0=1, kwargs...
)
    @assert 0 <= s0 <= 1
    s = s0 - sum(y_vector) / N
    f_output = _renewalfoutput(f, y_vector)
    return renewalequation_expectedcases(ρ, s, f_output; kwargs...)
end

function renewalequation_expectedcases!(
    f, y_matrix::AbstractMatrix, rho_matrix::AbstractMatrix, 
    initialvalues::AbstractMatrix, Ns::AbstractVector, psi; 
    kwargs...
)
    for g ∈ axes(y_matrix, 2)
        for t ∈ axes(y_matrix, 1)
            if t <= size(initialvalues, 1)
                y_matrix[t, g] = initialvalues[t, g]
            else
                y_matrix[t, g] = renewalequation_expectedcases(
                    f, y_matrix[1:(t-1), g], rho_matrix[t, g], Ns[g] * psi; 
                    kwargs...
                )
            end
        end
    end
end

renewalequation_poissoncases(expctcases::Real) = rand(Poisson(expctcases))

function renewalequation_poissoncases(ρ::Number, s::Number, f_output::Number)
    expctcases = renewalequation_expectedcases(ρ, s, f_output) 
    return renewalequation_poissoncases(expctcases)
end

function renewalequation_poissoncases(
    f, y_vector::AbstractVector, ρ::Number, N::Number; 
    psi=1, s0=1, kwargs...
)
    @assert 0 < psi <= 1
    @assert 0 <= s0 <= 1
    s = s0 - sum(y_vector) / (psi * N) 
    f_output = _renewalfoutput(f, y_vector)
    pcases = renewalequation_poissoncases(ρ, s, f_output)
    return min(round(Int, s * (N * psi), RoundDown), pcases)
end

function renewalequation_poissoncases!(
    f, y_matrix::AbstractMatrix, rho_matrix::AbstractMatrix, 
    initialvalues::AbstractMatrix, Ns::AbstractVector, psi; 
    columnseeds=nothing, kwargs...
)
    for g ∈ axes(y_matrix, 2)
        _setcolumnseed!(columnseeds, g)
        for t ∈ axes(y_matrix, 1)
            if t <= size(initialvalues, 1)
                y_matrix[t, g] = initialvalues[t, g]
            else
                pcases = renewalequation_poissoncases(
                    f, y_matrix[1:(t-1), g], rho_matrix[t, g], Ns[g] * psi; 
                    kwargs...
                )
                y_matrix[t, g] = pcases
            end
        end
    end
end

function _renewalfoutput(f::Function, i_vector)
    t = length(i_vector) + 1
    sumfi = zero(f(1) * i_vector[1])
    for τ in 1:(t - 1)
        sumfi += f(τ) * i_vector[t - τ]
    end
    return sumfi
end

function _renewalfoutput(f::AbstractVector, i_vector) 
    t = length(i_vector) + 1
    maxrun = min(t - 1, length(f))
    sumfi = zero(f[1] * i_vector[1])
    for τ in 1:maxrun
        sumfi += f[τ] * i_vector[t - τ]
    end
    return sumfi
end

function samplerenewalequation(f, chaindf, interventions; inds=axes(chaindf, 1), kwargs...)
    # matrices that will be updated / overwritten for each set of parameter values on the chain
    rho_matrix = zeros(size(interventions))
    y_matrix_det = zeros(size(interventions))
    y_matrix_poisson = zeros(Int, size(interventions))

    # a vector of these matrices, stored to allow calculation of mean and credible intervals
    rho_matrix_vec = Vector{Matrix{Float64}}(undef, size(chaindf, 1))
    y_matrix_det_vec = Vector{Matrix{Float64}}(undef, size(chaindf, 1))
    y_matrix_poisson_vec = Vector{Matrix{Int}}(undef, size(chaindf, 1))

    psi_vec = [ _samplepsi(chaindf, i) for i ∈ inds ]

    samplerenewalequation!(
        f, 
        rho_matrix, y_matrix_det, y_matrix_poisson, 
        rho_matrix_vec, y_matrix_det_vec, y_matrix_poisson_vec, 
        chaindf, interventions, psi_vec; 
        inds, kwargs...
    )
    return (
        psi_vec=psi_vec,
        rho_matrix=rho_matrix,
        rho_matrix_vec=rho_matrix_vec,
        y_matrix_det=y_matrix_det,
        y_matrix_det_vec=y_matrix_det_vec,
        y_matrix_poisson=y_matrix_poisson,
        y_matrix_poisson_vec=y_matrix_poisson_vec
    )
end

function samplerenewalequation!(
    f, 
    rho_matrix, y_matrix_det, y_matrix_poisson, 
    rho_matrix_vec, y_matrix_det_vec, y_matrix_poisson_vec, 
    chaindf, interventions, psi_vec; 
    initialvalues, Ns, 
    inds=axes(chaindf, 1), 
    secondaryinterventions=nothing, timeknots=nothing, timeperiods=nothing, 
    kwargs...
)
    testf(f)
    for ind ∈ inds
        samplerhomatrix!(
            rho_matrix, 
            chaindf, ind, interventions, secondaryinterventions, timeknots, timeperiods; 
            kwargs...
        )
        renewalequation_expectedcases!(
            f, y_matrix_det, rho_matrix, initialvalues, Ns, psi_vec[ind]; 
            kwargs...
        )
        renewalequation_poissoncases!(
            f, y_matrix_poisson, rho_matrix, initialvalues, Ns, psi_vec[ind]; 
            kwargs...
        )
        rho_matrix_vec[ind] = deepcopy(rho_matrix)
        y_matrix_det_vec[ind] = deepcopy(y_matrix_det)
        y_matrix_poisson_vec[ind] = deepcopy(y_matrix_poisson)
    end
end

function samplerhomatrix(chaindf, ind, interventions, args...; kwargs...)
    rho_matrix = zeros(size(interventions))
    samplerhomatrix!(rho_matrix, chaindf, ind, interventions, args...; kwargs...)
    return rho_matrix
end

function samplerhomatrix!(
    rho_matrix, chaindf, ind,
    interventions, secondaryinterventions, timeknots, timeperiods;
    logdelta=automatic, kwargs...
)
    nlocations = size(rho_matrix, 2)
    logzetas = _samplezetas(chaindf, ind, nlocations)
    logetas = _sampleetas(chaindf, timeknots, timeperiods, ind; kwargs...)
    logdeltavalue = _sampledelta(chaindf, logdelta, ind)
    for g ∈ axes(rho_matrix, 2), t ∈ axes(rho_matrix, 1)
        log_ρ = logzetas[g] + 
            _samplerhomatrixetavalue(logetas, timeperiods, t) +
            logdeltavalue * interventions[t, g] +
            _rhomatrixsecondaryinterventions(
                chaindf, secondaryinterventions, ind, t, g
            )
        rho_matrix[t, g] = exp(log_ρ)
    end
end

function samplerhomatrix!(
    rho_matrix, predictedw, chaindf, ind,
    interventions, secondaryinterventions, timeknots, timeperiods;
    logdelta=automatic, kwargs...
)
    nlocations = size(rho_matrix, 2)
    logzetas = _samplezetas(chaindf, ind, nlocations)
    logetas = _sampleetas(chaindf, timeknots, timeperiods, ind; kwargs...)
    logdeltavalue = _sampledelta(chaindf, logdelta, ind)
    for g ∈ axes(rho_matrix, 2), t ∈ axes(rho_matrix, 1)
        log_ρ = logzetas[g] + 
            _samplerhomatrixetavalue(logetas, timeperiods, t) +
            logdeltavalue * interventions[t, g] +
            _rhomatrixsecondaryinterventions(
                chaindf, secondaryinterventions, ind, t, g
            )
        rho_matrix[t, g] = exp(log_ρ)
    end
end

function _samplerhomatrixetavalue(logetas::AbstractVector, timeperiods::AbstractVector, t) 
    return logetas[timeperiods[t]]
end

_samplerhomatrixetavalue(logetas::CubicSpline, ::Any, t) = logetas[t]

function _samplezetas(chaindf, ind, nlocations)
    logzetas = [ 
        getproperty(chaindf, Symbol("logzeta_g$i.logzeta"))[ind] 
        for i ∈ 1:nlocations 
    ]
    return logzetas
end

function _sampleetas(
    chaindf, ::Nothing, timeperiods::AbstractVector, ind; 
    peakperiod=1, kwargs...
)
    ntimeperiods = _nunique(timeperiods)
    logetas = __sampleetas(chaindf, timeperiods, ntimeperiods, ind; peakperiod, kwargs...)
    return logetas
end

function _sampleetas(
    chaindf, timeknots::AbstractVector, ::Nothing, ind; 
    peakperiod=2, extrapl=zeros(1), extrapr=zeros(1), kwargs...
)
    ntimeperiods = length(timeknots)
    logetas = __sampleetas(chaindf, timeknots, ntimeperiods, ind; peakperiod, kwargs...)
    logtimespline = CubicSpline(timeknots, logetas; extrapl, extrapr)
    return logtimespline
end

function _sampleetas(
    chaindf, ::Nothing, ::Nothing, ind; 
    peakperiod=2, kwargs...
)
    ntimeperiods = length(timeknots)
    logetas = __sampleetas(chaindf, timeknots, ntimeperiods, ind; peakperiod, kwargs...)
    logtimespline = CubicSpline(timeknots, logetas; extrapl, extrapr)
    return logtimespline
end

function _sampleetas(::Any, ::Any, ::Nothing, ::Nothing, ::Any; kwargs...) 
    @error "Either timeknots or timeperiods must be supplied"
end

function __sampleetas(chaindf, timeperiods, ntimeperiods, ind; peakperiod, kwargs...)
    logetas = [
        i == peakperiod ? 
            0.0 : 
            getproperty(chaindf, Symbol("eta_t$i.logeta"))[ind] 
        for i ∈ 1:ntimeperiods
    ]
    return logetas
end

function _sampledelta(chaindf, ::Automatic, ind)
    logd = chaindf.logdelta[ind]
    return logd
end

_sampledelta(::Any, logdelta::Number, ::Any) = logdelta

_samplepsi(chaindf, ind) = chaindf.psi[ind]

_samplesecondarydeltas(::Any, ::Nothing, ::Any) = nothing

function _samplesecondarydeltas(chaindf, ::AbstractMatrix, ind)
    logd = chaindf.logsecondarydelta[ind]
    return exp(logd)
end

function _samplesecondarydeltas(
    chaindf, secondaryinterventions::Vector{<:AbstractMatrix}, ind
)
    ℓ = length(secondaryinterventions)
    logd = [ 
        getproperty(chaindf, Symbol("logsecondarydelta$i.logsecondarydelta"))[ind]
        for i ∈ 1:ℓ
    ]
    return exp.(logd)
end

_rhomatrixsecondaryinterventions(::Any, ::Nothing, ::Any, ::Any, ::Any) = 0.0

function _rhomatrixsecondaryinterventions(
    chaindf, secondaryinterventions::AbstractMatrix, ind, t, g
)
    logd = chaindf.logsecondarydelta[ind]
    return logd * secondaryinterventions[t, g]
end

function _rhomatrixsecondaryinterventions(
    chaindf, secondaryinterventions::Vector{<:AbstractMatrix}, ind, t, g
)
    ℓ = length(secondaryinterventions)
    logds = [ 
        getproperty(chaindf, Symbol("logsecondarydelta$i.logsecondarydelta"))[ind]
        for i ∈ 1:ℓ
    ]
    logsecondartinterventionsrho = sum([
        logd * secondint[t, g]
        for (logd, secondint) ∈ zip(logds, secondaryinterventions)
    ])
    return logsecondartinterventionsrho
end

_setcolumnseed!(::Nothing, ::Any) = nothing
_setcolumnseed!(columnseed::Number, ::Any) = Random.seed!(columnseed)
_setcolumnseed!(columnseeds::AbstractVector, g) = Random.seed!(columnseeds[g])

function testf(f::Function)
    m = minimum([ f(t) for t ∈ 1:1000 ])
    @assert m >= 0 "Function $f gives negative values in the first 1000 days"
    s = sum([ f(t) for t ∈ 1:1000 ])
    @info "Function $f sums to $s in the first 1000 days"
    return
end

function testf(f::AbstractVector)
    m = minimum(f)
    @assert m >= 0 "Vector f contains negative values"
    s = sum(f)
    @info "Vector f sums to $s"
    return
end

function samplerenewalequation_counterfactual(f, chaindf, interventions; kwargs...)
    samples = samplerenewalequation(f, chaindf, interventions; logdelta=0, kwargs...)
    return (
        rho_matrix_counterfactual=samples.rho_matrix,
        rho_matrix_vec_counterfactual=samples.rho_matrix_vec,
        y_matrix_det_counterfactual=samples.y_matrix_det,
        y_matrix_det_vec_counterfactual=samples.y_matrix_det_vec,
        y_matrix_poisson_counterfactual=samples.y_matrix_poisson,
        y_matrix_poisson_vec_counterfactual=samples.y_matrix_poisson_vec
    )
end

function samplerenewalequation_2sets(f, chaindf, interventions; columnseeds=1, kwargs...)
    results1 = samplerenewalequation(f, chaindf, interventions; columnseeds, kwargs...)
    results_counterfactual = samplerenewalequation_counterfactual(
        f, chaindf, interventions; 
        columnseeds, kwargs...
    )
    names1 = keys(results1)
    names2 = keys(results_counterfactual)
    allnames = ( names1..., names2... )
    values1 = values(results1)
    values2 = values(results_counterfactual)
    return NamedTuple{allnames}(( values1..., values2... ))
end
