
function stochasticiteration(
    transitionfunction::Function, u::AbstractVector, args...; 
    kwargs...
)
    return stochasticiteration(transitionfunction, default_rng(), u, args...; kwargs...)
end

function stochasticiteration(
    transitionfunction::Function, rng::AbstractRNG, 
    u::AbstractVector, t::Number, p, transitionmatrix::AbstractMatrix;
    broadcast_t=true, kwargs...
)
    # calculate expected rate of each transition
    if broadcast_t
        rates = transitionfunction(u, t, p)
    else 
        rates = transitionfunction(u, p)
    end

    return stochasticiteration(rates, rng, u, t, transitionmatrix; kwargs...)
end

function stochasticiteration(
    transitionfunction::Function, rng::AbstractRNG, 
    u::AbstractVector, t::Number, transitionmatrix::AbstractMatrix;
    broadcast_t=true, kwargs...
)
    # calculate expected rate of each transition
    if broadcast_t
        rates = transitionfunction(u, t)
    else 
        rates = transitionfunction(u)
    end

    return stochasticiteration(rates, rng, u, t, transitionmatrix; kwargs...)
end

function stochasticiteration(rates::AbstractVector, u::AbstractVector, args...; kwargs...)
    return stochasticiteration(rates, default_rng(), u, args...; kwargs...)
end

function stochasticiteration(
    rates::AbstractVector, rng::AbstractRNG, 
    u::AbstractVector, t::Number, transitionmatrix::AbstractMatrix;
    maxtstep=Inf,
)
    # check that there are the same number of compartments as columns in transitionmatrix
    @assert length(u) == size(transitionmatrix, 2)

    # check there are the same number of rates as rows in transitionmatrix
    @assert length(rates) == size(transitionmatrix, 1)

    # when is the change 
    rnum = rand(rng) 
    tstep = -log(rnum) / sum(rates)

    if tstep > maxtstep  # then nothing happens
        t += maxtstep
        return ( u, t )
    end

    # else, what happens 
    eventid = sample(rng, eachindex(rates), Weights(rates))

    # make the change 
    u += transitionmatrix[eventid, :]
    t += tstep 

    return ( u, t )
end

function stochasticmodel(transitionfunction::Function, u0::AbstractVector, args...; kwargs...) 
    return stochasticmodel(transitionfunction, default_rng(), u0, args...; kwargs...)
end

function stochasticmodel(
    transitionfunction::Function, rng::AbstractRNG, u0::AbstractVector, duration::T, args...;
    kwargs...
) where T <: Number
    tspan = ( one(T), duration )
    return stochasticmodel(transitionfunction, rng, u0, tspan, args...; kwargs...)
end

function stochasticmodel(
    transitionfunction::Function, rng::AbstractRNG, 
    u0::AbstractVector, tspan::Tuple{S, T}, args...;
    saveat=one(S), kwargs...
) where {S <: Number, T <: Number}
    tspanvector = tspan[1]:saveat:tspan[2] 
    return stochasticmodel(transitionfunction, rng, u0, tspanvector, args...; kwargs...)
end

function stochasticmodel(
    transitionfunction::Function, rng::AbstractRNG, 
    u0::AbstractVector{T}, tspanvector::AbstractVector, args...;
    saveall=false, kwargs...
) where T
    if saveall 
        return _saveall_stochasticmodel(
            transitionfunction, rng, u0, tspanvector, args...; 
            kwargs...
        )
    else 
        return _nosaveall_stochasticmodel(
            transitionfunction, rng, u0, tspanvector, args...; 
            kwargs...
        )
    end
end

function _saveall_stochasticmodel(
    transitionfunction::Function, rng::AbstractRNG, 
    u0::AbstractVector{T}, tspanvector::AbstractVector, args...;
    broadcast_t=true, maxtstep=Inf, kwargs...
) where T
    outputs = zeros(T, 1, length(u0))
    outputtimes = [ tspanvector[1] + 0.0 ]
    outputs[1, :] = u0 
    t = tspanvector[1]
    u = u0

    while t < last(tspanvector)
        timetolastt = last(tspanvector) - t
        u, t = stochasticiteration(
            transitionfunction, rng, u, t, args...; 
            broadcast_t, maxtstep=min(maxtstep, timetolastt)
        )
        outputs = vcat(outputs, u')
        push!(outputtimes, t)
    end

    return ( outputs=outputs, outputtimes=outputtimes )
end

function _nosaveall_stochasticmodel(
    transitionfunction::Function, rng::AbstractRNG, 
    u0::AbstractVector{T}, tspanvector::AbstractVector, args...;
    steptosaveat=true, kwargs...
) where T
    outputs = zeros(T, length(tspanvector), length(u0))
    outputs[1, :] = u0
    u = u0 
    if steptosaveat 
        _steptodaveat_stochasticmodel!(
            transitionfunction, rng, outputs, u, tspanvector, args...; 
            kwargs...
        )
    else
        _nosteptodaveat_stochasticmodel!(
            transitionfunction, rng, outputs, u, tspanvector, args...; 
            kwargs...
        )
    end
    return outputs
end

function _steptodaveat_stochasticmodel!(
    transitionfunction, rng, outputs, u, tspanvector, args...; 
    broadcast_t=true, maxtstep=Inf,
) 
    i = 1 
    t = tspanvector[i]
    while t < last(tspanvector)
        timetonextt = tspanvector[i+1] - t
        u, t = stochasticiteration(
            transitionfunction, rng, u, t, args...; 
            broadcast_t, maxtstep=min(maxtstep, timetonextt)
        )
        if t == tspanvector[i+1]
            i += 1 
            outputs[i, :] = u 
        end
    end
end

function _nosteptodaveat_stochasticmodel!(
    transitionfunction, rng, outputs, u, tspanvector, args...; 
    broadcast_t=true, maxtstep=Inf,
) 
    i = 1 
    t = tspanvector[i]
    while t < last(tspanvector)
        timetolastt = last(tspanvector) - t
        u, t = stochasticiteration(
            transitionfunction, rng, u, t, args...; 
            broadcast_t, maxtstep=min(maxtstep, timetonextt)
        )
        while t >= tspanvector[i+1]
            i += 1 
            outputs[i, :] = u 
        end
    end
end
