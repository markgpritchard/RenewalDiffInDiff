#

const largepop = Uniform(2_000_000, 60_000_000) 
const smallpop = Uniform(20_000, 60_000)

function simu0(population, proportionexposed; kwargs...)
    return simu0(default_rng(), population, proportionexposed; kwargs...)
end

function simu0(rng::AbstractRNG, population::Distribution, proportionexposed; kwargs...)
    n = round(Int, rand(rng, population))
    return simu0(rng, n, proportionexposed; kwargs...)
end

function simu0(rng::AbstractRNG, population::Integer, proportionexposed::Float64; minexposed=1)
    0 <= proportionexposed <= proportionexposed || throw(_propexposederror(proportionexposed))
    numberexposed = max(rand(rng, Binomial(population, proportionexposed)), minexposed)
    return simu0(rng, population, numberexposed)
end

function simu0(::AbstractRNG, population::Integer, numberexposed::Integer)
    return simulationu0(; s=(population - numberexposed), e=numberexposed)
end

function _propexposederror(proportionexposed)
    m = "$proportionexposed: proportion exposed must be between 0 and 1. If you intended to \
        pass `numberexposed`, this must be an `Integer` value"
    return ArgumentError(m)
end
