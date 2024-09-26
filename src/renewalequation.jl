
#=
function expectedinfections(f, i_vector, ρ_gt, s_gt) 
    sumfi = calculatesumfi(f, i_vector)
    return expectedinfections(sumfi, ρ_gt, s_gt)
end

#expectedinfections(sumfi, ρ_gt, s_gt) = *(sumfi, ρ_gt, s_gt)

function expectedinfections(sumfi, ρ_gt, s_gt) 
    results = *(sumfi, ρ_gt, s_gt)
    return results
end


function expectedinfections(f, i_vector, ρ_gt, N_g; s_g0=1) 
    sumfi = calculatesumfi(f, i_vector)
    s_gt = s_g0 - sum(i_vector) / N_g 
    return expectedinfections(sumfi, ρ_gt, s_gt)
end

function expectedinfections(f, y_vector, ρ_gt, N_g, ψ; kwargs...) 
    N_prop = N_g * ψ
    # note, this will return expected number of *reported* infections
    return expectedinfections(f, y_vector, ρ_gt, N_prop; kwargs...) 
end

expectedinfections(sumfi, ρ_gt, s_gt) = *(sumfi, ρ_gt, s_gt)

function calculatesumfi(f::Function, i_vector)
    t = length(i_vector) + 1
    sumfi = zero(f(1) * i_vector[1])
    for τ in 1:(t - 1)
        sumfi += f(τ) * i_vector[t - τ]
    end
    return sumfi
end

function calculatesumfi(f::AbstractVector, i_vector) 
    t = length(i_vector) + 1
    maxrun = min(t - 1, length(f))
    sumfi = zero(f[1] * i_vector[1])
    for τ in 1:maxrun
        sumfi += f[τ] * i_vector[t - τ]
    end
    return sumfi
end

poissoninfections(expctdinfections::Number) = rand(Poisson(expctdinfections))

function poissoninfections(args...; kwargs...) 
    expctdinfections = expectedinfections(args...; kwargs...)
    return poissoninfections(expctdinfections)
end

const InfectionFunction = Union{typeof(expectedinfections), typeof(poissoninfections)}

function runrenewalequation(inffunc::InfectionFunction, f, interventions; kwargs...)
    rho_matrix = zeros(size(interventions))
    expectedinfections_matrix = zeros(size(interventions))
    i_matrix = zeros(size(interventions))
    runrenewalequation!(
        inffunc, f, rho_matrix, expectedinfections_matrix, i_matrix, interventions; 
        kwargs...
    )
    return (
        expectedinfections_matrix=expectedinfections_matrix, 
        i_matrix=i_matrix, 
        rho_matrix=rho_matrix
    )
end

function runrenewalequation(args...; kwargs...)
    return runrenewalequation(poissoninfections, args...; kwargs...)
end

function runrenewalequation!(
    inffunc::InfectionFunction, f, rho_matrix, expectedinfections_matrix, i_matrix, interventions; 
    etas, initialvalues, Ns, zetas,
    columnseeds=nothing, delta=0, psi=1, S0s=Ns, timeknots=nothing, timeperiods=nothing, 
    kwargs...
)
    testf(f)
    n_initialvalues = size(initialvalues, 1)
    _runrenewalequation!(
        inffunc, f, rho_matrix, expectedinfections_matrix, i_matrix, interventions,
        etas, initialvalues, Ns, zetas,
        delta, n_initialvalues, psi, S0s, timeknots, timeperiods, columnseeds;
        kwargs...
    )
end

function runrenewalequation!(args...; kwargs...)
    return runrenewalequation!(poissoninfections, args...; kwargs...)
end

function _runrenewalequation!(
    inffunc::InfectionFunction, f, rho_matrix, expectedinfections_matrix, i_matrix, interventions,
    etas::AbstractVector{T}, initialvalues, Ns, zetas::AbstractVector{T},
    delta, n_initialvalues, psi, S0s, ::Nothing, timeperiods, columnseeds;
    secondaryinterventions=nothing, secondarydeltavalues=nothing,
) where T
    _calculaterhomatrix!(
        rho_matrix, i_matrix, interventions, 
        secondaryinterventions, timeperiods, zetas, etas, delta, secondarydeltavalues
    )
    __runrenewalequation!(
        inffunc, f, rho_matrix, expectedinfections_matrix, i_matrix, 
        S0s, initialvalues, n_initialvalues, Ns, psi, columnseeds
    )
end

function _runrenewalequation!(
    inffunc::InfectionFunction, f, rho_matrix, expectedinfections_matrix, i_matrix, interventions,
    etas::AbstractVector{T}, initialvalues, Ns, zetas::AbstractVector{T},
    delta, n_initialvalues, psi, S0s, timeknots::AbstractVector, ::Nothing, columnseeds;
    secondaryinterventions=nothing, secondarydeltavalues=nothing,
    extrapl=zeros(1), extrapr=zeros(1), 
) where T
    spline = CubicSpline(timeknots, etas; extrapl, extrapr)
    
    _calculaterhomatrix!(
        rho_matrix, i_matrix, interventions, 
        secondaryinterventions, nothing, zetas, spline, delta, secondarydeltavalues
    )
    __runrenewalequation!(
        inffunc, f, rho_matrix, expectedinfections_matrix, i_matrix, 
        S0s, initialvalues, n_initialvalues, Ns, psi, columnseeds
    )
end

function __runrenewalequation!(
    inffunc::InfectionFunction, f, rho_matrix, expectedinfections_matrix, i_matrix, 
    S0s, initialvalues, n_initialvalues, Ns, psi, columnseeds
)
    for g ∈ axes(i_matrix, 2)
        _setcolumnseed!(columnseeds, g)  # re-seed for each column of the matrix
        for t ∈ axes(i_matrix, 1)
            if t <= n_initialvalues
                expectedinfections_matrix[t, g] = NaN
                i_matrix[t, g] = initialvalues[t, g]
            else
                expectedinfections_matrix[t, g] = expectedinfections(
                    f, 
                    view(i_matrix, 1:(t - 1), g), 
                    rho_matrix[t, g], 
                    S0s[g] / Ns[g] - sum(i_matrix[1:(t - 1), g]) / (psi * Ns[g])
                )
                i_matrix[t, g] = _renewalequationobservation(
                    inffunc, expectedinfections_matrix[t, g]
                )
            end
        end
    end
end

function _renewalequationobservation(::typeof(expectedinfections), expctdinfections)
    return expctdinfections
end

function _renewalequationobservation(inffunc::typeof(poissoninfections), expctdinfections)
    return inffunc(expctdinfections)
end

_setcolumnseed!(::Nothing, ::Any) = nothing
_setcolumnseed!(columnseed::Number, ::Any) = Random.seed!(columnseed)
_setcolumnseed!(columnseeds::AbstractVector, g) = Random.seed!(columnseeds[g])

function _calculaterhomatrix!(
    rho_matrix, i_matrix, interventions, 
    secondaryinterventions, timeperiods::AbstractVector, 
    zetas, etas::Vector, deltavalue, secondarydeltavalues
)
    for g ∈ axes(i_matrix, 2), t ∈ axes(i_matrix, 1)
        rho_matrix[t, g] = zetas[g] * 
            etas[timeperiods[t]] * 
            deltavalue^interventions[t, g] * 
            _rhomatrixsecondaryinterventions(
                secondaryinterventions, secondarydeltavalues, t, g
            )
    end
end

function _calculaterhomatrix!(
    rho_matrix, i_matrix, interventions, 
    secondaryinterventions, ::Nothing,
    zetas, spline::CubicSpline, deltavalue, secondarydeltavalues
)
    for g ∈ axes(i_matrix, 2), t ∈ axes(i_matrix, 1)
        rho_matrix[t, g] = zetas[g] * 
            spline[t] * 
            deltavalue^interventions[t, g] *
            _rhomatrixsecondaryinterventions(
                secondaryinterventions, secondarydeltavalues, t, g
            )
    end
end

function _calculaterhomatrix!(
    rho_matrix, i_matrix, interventions, 
    secondaryinterventions, ::Nothing, 
    zetas, etas::Vector, deltavalue, secondarydeltavalues
)
    for g ∈ axes(i_matrix, 2), t ∈ axes(i_matrix, 1)
        rho_matrix[t, g] = zetas[g] * 
            etas[t] * 
            deltavalue^interventions[t, g] * 
            _rhomatrixsecondaryinterventions(
                secondaryinterventions, secondarydeltavalues, t, g
            )
    end
end

_rhomatrixsecondaryinterventions(::Nothing, ::Nothing, ::Any, ::Any) = 1.0

function _rhomatrixsecondaryinterventions(
    secondaryinterventions::AbstractMatrix, secondarydeltavalues::Number, t, g
)
    return secondarydeltavalues^secondaryinterventions[t, g]
end

function _rhomatrixsecondaryinterventions(
    secondaryinterventions::Vector{<:AbstractMatrix}, secondarydeltavalues::Vector, t, g
)
    multiplier = prod([
        _rhomatrixsecondaryinterventions(secint, secdelta, t, g)
        for (secint, secdelta) ∈ zip(secondaryinterventions, secondarydeltavalues)
    ])
    return multiplier
end

function runrenewalequationsamples(args...; kwargs...)
    return runrenewalequationsamples(poissoninfections, args...; kwargs...)
end

function runrenewalequationsamples(inffunc::InfectionFunction, f, interventions, chain::Chains, args...; kwargs...)
    chaindf = DataFrame(chain)
    return runrenewalequationsamples(inffunc, f, interventions, chaindf, args...; kwargs...)
end

function runrenewalequationsamples(
    inffunc::InfectionFunction, f, interventions, chaindf::DataFrame, psi, nsamples=500; 
    initialseed=10, kwargs...
)
    expectinf_matrix = zeros(Float64, size(interventions))
    expectinf_matrixvector = Vector{typeof(expectinf_matrix)}(undef, nsamples)
    expectinf_matrix_counterfactual = Vector{typeof(expectinf_matrix)}(undef, nsamples)
    i_matrix = zeros(Float64, size(interventions))
    i_matrixvector = Vector{typeof(i_matrix)}(undef, nsamples)
    i_matrixvector_counterfactual = Vector{typeof(i_matrix)}(undef, nsamples)
    ρ_matrix = zeros(Float64, size(interventions))
    ρ_matrixvector = Vector{typeof(ρ_matrix)}(undef, nsamples)
    ρ_matrixvector_counterfactual = Vector{typeof(ρ_matrix)}(undef, nsamples)
    predictedw_matrix = zeros(Float64, size(interventions))
    predictedw_matrixvector = Vector{typeof(predictedw_matrix)}(undef, nsamples)
    ind_vector = Vector{Int}(undef, nsamples)
    seed = initialseed

    for i ∈ 1:nsamples
        columnseeds = [ seed + j for j ∈ axes(interventions, 2) ]
        ind = sample(axes(chaindf, 1))
        Random.seed!(seed)
        runrenewalequationsample!(
            inffunc, f, expectinf_matrix, i_matrix, ρ_matrix, predictedw_matrix, 
            interventions, chaindf, psi, ind; 
            columnseeds, kwargs...
        )
        expectinf_matrixvector[i] = deepcopy(expectinf_matrix)
        i_matrixvector[i] = deepcopy(i_matrix)
        ρ_matrixvector[i] = deepcopy(ρ_matrix)
        predictedw_matrixvector[i] = deepcopy(predictedw_matrix)
        # repeat with the same seed and delta=0 for the counterfactual
        Random.seed!(seed)
        runrenewalequationsample!(
            f, expectinf_matrix, i_matrix, ρ_matrix, predictedw_matrix, 
            interventions, chaindf, psi, ind; 
            columnseeds, delta=0.0, kwargs...
        )
        expectinf_matrix_counterfactual[i] = deepcopy(expectinf_matrix)
        i_matrixvector_counterfactual[i] = deepcopy(i_matrix)
        ρ_matrixvector_counterfactual[i] = deepcopy(ρ_matrix)
        ind_vector[i] = ind
        seed += 10
    end
    return (
        expectinf_matrix=expectinf_matrix,
        expectinf_matrixvector=expectinf_matrixvector,
        expectinf_matrix_counterfactual=expectinf_matrix_counterfactual,
        ind_vector=ind_vector,
        i_matrixvector=i_matrixvector, 
        i_matrixvector_counterfactual=i_matrixvector_counterfactual, 
        predictedw_matrixvector=predictedw_matrixvector, 
        ρ_matrixvector=ρ_matrixvector, 
        ρ_matrixvector_counterfactual=ρ_matrixvector_counterfactual
    )
end

function runrenewalequationsample!(args...; kwargs...)
    return runrenewalequationsample!(poissoninfections, args...; kwargs...)
end

function runrenewalequationsample!(
    inffunc::InfectionFunction, f, expectinf_matrix, i_matrix, ρ_matrix, predictedw_matrix, 
    interventions, chaindf::DataFrame, psi, ind; 
    initialvalues, Ns, columnseeds=nothing, delta=automatic, 
    n_initialvalues=size(initialvalues, 1), S0s=Ns, 
    secondaryinterventions=nothing,
    timeknots=nothing, timeperiods=nothing, kwargs...
)
    zetas = _samplezetas(interventions, chaindf, ind)
    etas = _sampleetas(i_matrix, chaindf, timeknots, timeperiods, ind; kwargs...)
    deltavalue = _sampledelta(chaindf, delta, ind)
    secondarydeltavalues = _samplesecondarydeltas(chaindf, secondaryinterventions, ind)

    _calculaterhomatrix!(
        ρ_matrix, i_matrix, interventions, secondaryinterventions, timeperiods, zetas, etas, deltavalue, secondarydeltavalues
    )
    _runrenewalequation!(
        inffunc, f, ρ_matrix, expectinf_matrix, i_matrix, interventions,
        etas, initialvalues, Ns, zetas,
        deltavalue, n_initialvalues, psi, S0s, nothing, timeperiods, columnseeds;
        secondaryinterventions, secondarydeltavalues,
    )
    _predictw!(
        predictedw_matrix, i_matrix, 
        interventions, timeperiods, Ns, zetas, etas, deltavalue, psi
    )
end

function _samplezetas(interventions, chaindf, ind)
    nlocations = size(interventions, 2)
    logzetas = [ 
        getproperty(chaindf, Symbol("logzeta_g$i.logzeta"))[ind] 
        for i ∈ 1:nlocations 
    ]
    return exp.(logzetas)
end

function _sampleetas(
    i_matrix, chaindf, ::Nothing, timeperiods::AbstractVector, ind; 
    peakperiod=1, kwargs...
)
    ntimeperiods = _nunique(timeperiods)
    logetas = __sampleetas(chaindf, timeperiods, ntimeperiods, ind; peakperiod, kwargs...)
    return exp.(logetas)
end

function _sampleetas(
    i_matrix, chaindf, timeknots::AbstractVector, ::Nothing, ind; 
    peakperiod=2, extrapl=zeros(1), extrapr=zeros(1), kwargs...
)
    ntimeperiods = length(timeknots)
    logetas = __sampleetas(chaindf, timeknots, ntimeperiods, ind; peakperiod, kwargs...)
    timespline = CubicSpline(timeknots, logetas; extrapl, extrapr)
    returntimeperiods = collect(axes(i_matrix, 1))
    etas = [ exp(timespline(t)) for t ∈ returntimeperiods ]
    return etas
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
    return exp(logd)
end

_sampledelta(::Any, delta::Number, ::Any) = exp(delta)

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

function _predictw!(
    predictedw_matrix, i_matrix, 
    interventions, timeperiods::AbstractVector, Ns, zetas, etas::Vector, delta, psi
)
    for g ∈ axes(i_matrix, 2), t ∈ axes(i_matrix, 1)
        predictedw_matrix[t, g] = log(zetas[g]) + 
            log(etas[timeperiods[t]]) + 
            interventions[t, g] * log(delta) #
    end
end

function _predictw!(
    predictedw_matrix, i_matrix, 
    interventions, ::Nothing, Ns, zetas, spline::CubicSpline, delta, psi
)
    for g ∈ axes(i_matrix, 2), t ∈ axes(i_matrix, 1)
        predictedw_matrix[t, g] = log(zetas[g]) + 
            log(spline[t]) + 
            interventions[t, g] * log(delta) 
    end
end

function _predictw!(
    predictedw_matrix, i_matrix, 
    interventions, ::Nothing, Ns, zetas, etas::Vector, delta, psi
)
    for g ∈ axes(i_matrix, 2), t ∈ axes(i_matrix, 1)
        predictedw_matrix[t, g] = log(zetas[g]) + 
            log(etas[t]) + 
            interventions[t, g] * log(delta) 
    end
end

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
=#

#=

renewalequation_expectedcases(ρ::Number, s::Number, f_output::Number) = *(ρ, s, f_output)

function renewalequation_expectedcases(args...; kwargs...) 
    return _genericrenewalequation(renewalequation_expectedcases, args...; kwargs...)
end

function renewalequation_expectedcases!(
    f, y_matrix::AbstractMatrix, rho_matrix::AbstractMatrix, 
    initialvalues::AbstractMatrix, Ns::AbstractVector; 
    kwargs...
)
    for g ∈ axes(y_matrix, 2)
        for t ∈ axes(y_matrix, 1)
            if t <= size(initialvalues, 1)
                y_matrix[t, g] = initialvalues[t, g]
            else
                y_matrix[t, g] = renewalequation_expectedcases(
                    f, y_matrix[1:(t-1), g], rho_matrix[t, g], Ns[g]; 
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

function renewalequation_poissoncases(args...; kwargs...) 
    return _genericrenewalequation(renewalequation_poissoncases, args...; kwargs...)
end

function renewalequation_poissoncases!(
    f, y_matrix::AbstractMatrix, rho_matrix::AbstractMatrix, 
    initialvalues::AbstractMatrix, Ns::AbstractVector; 
    columnseeds=nothing, kwargs...
)
    for g ∈ axes(y_matrix, 2)
        _setcolumnseed!(columnseeds, g)
        for t ∈ axes(y_matrix, 1)
            if t <= size(initialvalues, 1)
                y_matrix[t, g] = initialvalues[t, g]
            else
                expctcases = renewalequation_expectedcases(
                    f, y_matrix[1:(t-1), g], rho_matrix[t, g], Ns[g]; 
                    kwargs...
                )
                y_matrix[t, g] = renewalequation_poissoncases(expctcases)
            end
        end
    end
end

function _genericrenewalequation(
    renfunc, f, y_vector::AbstractVector, ρ::Number, N::Number; 
    psi=1, s0=1, kwargs...
)
    @assert 0 < psi <= 1
    @assert 0 <= s0 <= 1
    s = s0 - sum(y_vector) / (psi * N) 
    f_output = _renewalfoutput(f, y_vector)
    return renfunc(ρ, s, f_output; kwargs...)
end

function _genericrenewalequation(
    frenfunc, f, 
    rho_matrix::AbstractMatrix{T}, initialvalues::AbstractMatrix, Ns::AbstractVector; 
    kwargs...
) where T
    y_matrix = Matrix{T}(undef, size(rho_matrix))
    _genericrenewalequation!(renfunc, f, y_matrix, rho_matrix, initialvalues, Ns; kwargs...)
    return y_matrix
end

function _genericrenewalequation(renfunc, args...; kwargs...)
    return renfunc(args...; kwargs...)
end

function _genericrenewalequation!(::typeof(renewalequation_expectedcases), args...; kwargs...)
    renewalequation_expectedcases!(args...; kwargs...)
end

function _genericrenewalequation!(::typeof(renewalequation_poissoncases), args...; kwargs...)
    renewalequation_poissoncases!(args...; kwargs...)
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

function samplerenewalequation(f, chaindf, interventions; nsamples=:all, kwargs...)
    rho_matrix = zeros(size(interventions))
    y_matrix_det = zeros(size(interventions))
    y_matrix_poisson = zeros(Int, size(interventions))
    return samplerenewalequation(
        f, chaindf, interventions, rho_matrix, y_matrix_det, y_matrix_poisson, nsamples; 
        kwargs...
    )
end

function samplerenewalequation(
    f, chaindf, interventions, rho_matrix, y_matrix_det, y_matrix_poisson, nsamples::Symbol; 
    kwargs...
)
    @assert nsamples == :all
    inds = axes(chaindf, 1)
    rho_matrix_vec = Vector{typeof(rho_matrix)}(undef, size(chaindf, 1))
    y_matrix_det_vec = Vector{typeof(y_matrix_det)}(undef, size(chaindf, 1))
    y_matrix_poisson_vec = Vector{typeof(y_matrix_poisson)}(undef, size(chaindf, 1))
    samplerenewalequation!(
        f, 
        rho_matrix, y_matrix_det, y_matrix_poisson, 
        rho_matrix_vec, y_matrix_det_vec, y_matrix_poisson_vec, 
        chaindf, interventions, inds; 
        kwargs...
    )
    return (
        rho_matrix=rho_matrix,
        rho_matrix_vec=rho_matrix_vec,
        y_matrix_det=y_matrix_det,
        y_matrix_det_vec=y_matrix_det_vec,
        y_matrix_poisson=y_matrix_poisson,
        y_matrix_poisson_vec=y_matrix_poisson_vec
    )
end

function samplerenewalequation(
    f, chaindf, interventions, rho_matrix, y_matrix_det, y_matrix_poisson, nsamples::Integer; 
    kwargs...
)
    inds = [ sample(axes(chaindf, 1)) for _ ∈ nsamples ]
    rho_matrix_vec = Vector{typeof(rho_matrix)}(undef, nsamples)
    y_matrix_det_vec = Vector{typeof(y_matrix_det)}(undef, nsamples)
    y_matrix_poisson_vec = Vector{typeof(y_matrix_poisson)}(undef, nsamples)   
    rho_matrix_vec = Vector{typeof(rho_matrix)}(undef, nsamples)
    y_matrix_det_vec = Vector{typeof(y_matrix_det)}(undef, nsamples)
    y_matrix_poisson_vec = Vector{typeof(y_matrix_poisson)}(undef, nsamples)
    samplerenewalequation!(
        f, 
        rho_matrix, y_matrix_det, y_matrix_poisson, 
        rho_matrix_vec, y_matrix_det_vec, y_matrix_poisson_vec, 
        chaindf, interventions, inds; 
        kwargs...
    )
    return (
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
    chaindf, interventions, inds; 
    initialvalues, Ns, 
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
            f, y_matrix_det, rho_matrix, initialvalues, Ns; 
            kwargs...
        )
        renewalequation_poissoncases!(
            f, y_matrix_poisson, rho_matrix, initialvalues, Ns; 
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
    return exp(logd)
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
    @unpack rho_matrix,rho_matrix_vec, y_matrix_det, y_matrix_det_vec, y_matrix_poisson, y_matrix_poisson_vec =
        samplerenewalequation(f, chaindf, interventions; logdelta=0, kwargs...)
    return (
        rho_matrix_counterfactual=rho_matrix,
        rho_matrix_vec_counterfactual=rho_matrix_vec,
        y_matrix_det_counterfactual=y_matrix_det,
        y_matrix_det_vec_counterfactual=y_matrix_det_vec,
        y_matrix_poisson_counterfactual=y_matrix_poisson,
        y_matrix_poisson_vec_counterfactual=y_matrix_poisson_vec
    )
end

function samplerenewalequation_2sets(f, chaindf, interventions; kwargs...)
    results1 = samplerenewalequation(f, chaindf, interventions; kwargs...)
    results_counterfactual = samplerenewalequation_counterfactual(
        f, chaindf, interventions; 
        kwargs...
    )
    names1 = keys(results1)
    names2 = keys(results_counterfactual)
    values1 = values(results1)
    values2 = values(results_counterfactual)
    return NamedTuple{( names1..., names2... )}( values1..., values2... )
end

=#

renewalequation_expectedcases(ρ::Number, s::Number, f_output::Number; kwargs...) = *(ρ, s, f_output)

function renewalequation_expectedcases(
    f, y_vector::AbstractVector, ρ::Number, N::Number; 
    psi=1, s0=1, kwargs...
)
    @assert 0 < psi <= 1
    @assert 0 <= s0 <= 1
    s = s0 - sum(y_vector) / (psi * N) 
    f_output = _renewalfoutput(f, y_vector)
    return renewalequation_expectedcases(ρ, s, f_output; kwargs...)
end

function renewalequation_expectedcases(
    f, rho_matrix::AbstractMatrix{T}, initialvalues::AbstractMatrix, Ns::AbstractVector; 
    kwargs...
) where T
    y_matrix = Matrix{T}(undef, size(rho_matrix))
    renewalequation_expectedcases!(f, y_matrix, rho_matrix, initialvalues, Ns; kwargs...)
    return y_matrix
end

function renewalequation_expectedcases!(
    f, y_matrix::AbstractMatrix, rho_matrix::AbstractMatrix, 
    initialvalues::AbstractMatrix, Ns::AbstractVector; 
    kwargs...
)
    for g ∈ axes(y_matrix, 2)
        for t ∈ axes(y_matrix, 1)
            if t <= size(initialvalues, 1)
                y_matrix[t, g] = initialvalues[t, g]
            else
                y_matrix[t, g] = renewalequation_expectedcases(
                    f, y_matrix[1:(t-1), g], rho_matrix[t, g], Ns[g]; 
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
    #println("rho=$ρ, s=$s, f_output=$f_output, N=$N, sn=$(s*N), min1=$(round(Int, s * N, RoundDown)), psf = $(*(ρ, s, f_output))")
    pcases = renewalequation_poissoncases(ρ, s, f_output)
    return min(round(Int, s * (N * psi), RoundDown), pcases)
end

function renewalequation_poissoncases(
    f, rho_matrix::AbstractMatrix{T}, initialvalues::AbstractMatrix, Ns::AbstractVector; 
    kwargs...
) where T
    y_matrix = Matrix{T}(undef, size(rho_matrix))
    renewalequation_poissoncases!(f, y_matrix, rho_matrix, initialvalues, Ns; kwargs...)
    return y_matrix
end

function renewalequation_poissoncases!(
    f, y_matrix::AbstractMatrix, rho_matrix::AbstractMatrix, 
    initialvalues::AbstractMatrix, Ns::AbstractVector; 
    columnseeds=nothing, kwargs...
)
    for g ∈ axes(y_matrix, 2)
        _setcolumnseed!(columnseeds, g)
        for t ∈ axes(y_matrix, 1)
            if t <= size(initialvalues, 1)
                y_matrix[t, g] = initialvalues[t, g]
            else
                #expctcases = renewalequation_expectedcases(
                #    f, y_matrix[1:(t-1), g], rho_matrix[t, g], Ns[g]; 
                #    kwargs...
                #)
                pcases = renewalequation_poissoncases(
                    f, y_matrix[1:(t-1), g], rho_matrix[t, g], Ns[g]; 
                    kwargs...
                )
                #println(pcases)
                y_matrix[t, g] = pcases# renewalequation_poissoncases(expctcases)
            end
        end
    end
end
#=
function _genericrenewalequation(
    renfunc, f, y_vector::AbstractVector, ρ::Number, N::Number; 
    psi=1, s0=1, kwargs...
)
    @assert 0 < psi <= 1
    @assert 0 <= s0 <= 1
    s = s0 - sum(y_vector) / (psi * N) 
    f_output = _renewalfoutput(f, y_vector)
    return renfunc(ρ, s, f_output; kwargs...)
end

function _genericrenewalequation(
    frenfunc, f, 
    rho_matrix::AbstractMatrix{T}, initialvalues::AbstractMatrix, Ns::AbstractVector; 
    kwargs...
) where T
    y_matrix = Matrix{T}(undef, size(rho_matrix))
    _genericrenewalequation!(renfunc, f, y_matrix, rho_matrix, initialvalues, Ns; kwargs...)
    return y_matrix
end

function _genericrenewalequation(renfunc, args...; kwargs...)
    return renfunc(args...; kwargs...)
end

function _genericrenewalequation!(::typeof(renewalequation_expectedcases), args...; kwargs...)
    renewalequation_expectedcases!(args...; kwargs...)
end

function _genericrenewalequation!(::typeof(renewalequation_poissoncases), args...; kwargs...)
    renewalequation_poissoncases!(args...; kwargs...)
end
=#
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

function samplerenewalequation(f, chaindf, interventions; nsamples=:all, kwargs...)
    rho_matrix = zeros(size(interventions))
    y_matrix_det = zeros(size(interventions))
    y_matrix_poisson = zeros(Int, size(interventions))
    return samplerenewalequation(
        f, chaindf, interventions, rho_matrix, y_matrix_det, y_matrix_poisson, nsamples; 
        kwargs...
    )
end

function samplerenewalequation(
    f, chaindf, interventions, rho_matrix, y_matrix_det, y_matrix_poisson, nsamples::Symbol; 
    kwargs...
)
    @assert nsamples == :all
    inds = axes(chaindf, 1)
    rho_matrix_vec = Vector{typeof(rho_matrix)}(undef, size(chaindf, 1))
    y_matrix_det_vec = Vector{typeof(y_matrix_det)}(undef, size(chaindf, 1))
    y_matrix_poisson_vec = Vector{typeof(y_matrix_poisson)}(undef, size(chaindf, 1))
    samplerenewalequation!(
        f, 
        rho_matrix, y_matrix_det, y_matrix_poisson, 
        rho_matrix_vec, y_matrix_det_vec, y_matrix_poisson_vec, 
        chaindf, interventions, inds; 
        kwargs...
    )
    psi_vec = [ _samplepsi(chaindf, i) for i ∈ inds ]
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

function samplerenewalequation(
    f, chaindf, interventions, rho_matrix, y_matrix_det, y_matrix_poisson, nsamples::Integer; 
    kwargs...
)
    inds = [ sample(axes(chaindf, 1)) for _ ∈ nsamples ]
    rho_matrix_vec = Vector{typeof(rho_matrix)}(undef, nsamples)
    y_matrix_det_vec = Vector{typeof(y_matrix_det)}(undef, nsamples)
    y_matrix_poisson_vec = Vector{typeof(y_matrix_poisson)}(undef, nsamples)   
    rho_matrix_vec = Vector{typeof(rho_matrix)}(undef, nsamples)
    y_matrix_det_vec = Vector{typeof(y_matrix_det)}(undef, nsamples)
    y_matrix_poisson_vec = Vector{typeof(y_matrix_poisson)}(undef, nsamples)
    samplerenewalequation!(
        f, 
        rho_matrix, y_matrix_det, y_matrix_poisson, 
        rho_matrix_vec, y_matrix_det_vec, y_matrix_poisson_vec, 
        chaindf, interventions, inds; 
        kwargs...
    )
    psi_vec = [ _samplepsi(chaindf, i) for i ∈ inds ]
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
    chaindf, interventions, inds; 
    initialvalues, Ns, 
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
            f, y_matrix_det, rho_matrix, initialvalues, Ns; 
            kwargs...
        )
        renewalequation_poissoncases!(
            f, y_matrix_poisson, rho_matrix, initialvalues, Ns; 
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
    @unpack rho_matrix,rho_matrix_vec, y_matrix_det, y_matrix_det_vec, y_matrix_poisson, y_matrix_poisson_vec =
        samplerenewalequation(f, chaindf, interventions; logdelta=0, kwargs...)
    return (
        rho_matrix_counterfactual=rho_matrix,
        rho_matrix_vec_counterfactual=rho_matrix_vec,
        y_matrix_det_counterfactual=y_matrix_det,
        y_matrix_det_vec_counterfactual=y_matrix_det_vec,
        y_matrix_poisson_counterfactual=y_matrix_poisson,
        y_matrix_poisson_vec_counterfactual=y_matrix_poisson_vec
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

## Version to estimate R_t

renewalequation_expectedcases_rt(ρ::Number, f_output::Number; kwargs...) = *(ρ, f_output)

function renewalequation_expectedcases_rt(
    f, y_vector::AbstractVector, ρ::Number; 
    kwargs...
)
    f_output = _renewalfoutput(f, y_vector)
    return renewalequation_expectedcases(ρ, f_output; kwargs...)
end

function renewalequation_expectedcases_rt(
    f, rho_matrix::AbstractMatrix{T}, initialvalues::AbstractMatrix; 
    kwargs...
) where T
    y_matrix = Matrix{T}(undef, size(rho_matrix))
    renewalequation_expectedcases!(f, y_matrix, rho_matrix, initialvalues; kwargs...)
    return y_matrix
end

function renewalequation_expectedcases_rt!(
    f, y_matrix::AbstractMatrix, rho_matrix::AbstractMatrix, 
    initialvalues::AbstractMatrix; 
    kwargs...
)
    for g ∈ axes(y_matrix, 2)
        for t ∈ axes(y_matrix, 1)
            if t <= size(initialvalues, 1)
                y_matrix[t, g] = initialvalues[t, g]
            else
                y_matrix[t, g] = renewalequation_expectedcases_rt(
                    f, y_matrix[1:(t-1), g], rho_matrix[t, g]; 
                    kwargs...
                )
            end
        end
    end
end

renewalequation_poissoncases_rt(expctcases::Real) = renewalequation_poissoncases(expctcases)

function renewalequation_poissoncases_rt(ρ::Number, f_output::Number)
    expctcases = renewalequation_expectedcases_rt(ρ, f_output) 
    return renewalequation_poissoncases_rt(expctcases)
end

function renewalequation_poissoncases_rt(
    f, y_vector::AbstractVector, ρ::Number; 
    kwargs...
)
    f_output = _renewalfoutput(f, y_vector)
    #println("rho=$ρ, s=$s, f_output=$f_output, N=$N, sn=$(s*N), min1=$(round(Int, s * N, RoundDown)), psf = $(*(ρ, s, f_output))")
    return renewalequation_poissoncases_rt(ρ, f_output)
end

function renewalequation_poissoncases_rt(
    f, rho_matrix::AbstractMatrix{T}, initialvalues::AbstractMatrix; 
    kwargs...
) where T
    y_matrix = Matrix{T}(undef, size(rho_matrix))
    renewalequation_poissoncases_rt!(f, y_matrix, rho_matrix, initialvalues; kwargs...)
    return y_matrix
end

function renewalequation_poissoncases_rt!(
    f, y_matrix::AbstractMatrix, rho_matrix::AbstractMatrix, 
    initialvalues::AbstractMatrix; 
    columnseeds=nothing, kwargs...
)
    for g ∈ axes(y_matrix, 2)
        _setcolumnseed!(columnseeds, g)
        for t ∈ axes(y_matrix, 1)
            if t <= size(initialvalues, 1)
                y_matrix[t, g] = initialvalues[t, g]
            else
                #expctcases = renewalequation_expectedcases(
                #    f, y_matrix[1:(t-1), g], rho_matrix[t, g], Ns[g]; 
                #    kwargs...
                #)
                pcases = renewalequation_poissoncases_rt(
                    f, y_matrix[1:(t-1), g], rho_matrix[t, g]; 
                    kwargs...
                )
                #println(pcases)
                y_matrix[t, g] = pcases# renewalequation_poissoncases(expctcases)
            end
        end
    end
end
