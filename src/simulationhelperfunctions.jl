#

const largepop = Uniform(2_000_000, 60_000_000) 
const smallpop = Uniform(20_000, 60_000)

simu0(population, proportionexposed) = simu0(default_rng(), population, proportionexposed)

function simu0(rng::AbstractRNG, population::Distribution, proportionexposed)
    n = round(Int, rand(rng, population))
    return simu0(rng, n, proportionexposed)
end

function simu0(rng::AbstractRNG, population::Integer, proportionexposed::Float64)
    0 <= proportionexposed <= proportionexposed || throw(_propexposederror(proportionexposed))
    numberexposed = rand(rng, Binomial(population, proportionexposed))
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
