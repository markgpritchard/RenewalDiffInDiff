
module AnalysisFunctions 

using DataFrames, Distributions, DrWatson, Random
import Base: getindex, length, size, sum

export SEIRCompartments, SEIRParameters
export COVIDSERIALINTERVAL, POPULATION2020, REGIONCODEDF
export runseir_deterministic, runseir_noisy, seir_deterministic, seir_noisy
export anyminimum, getindex, length, size, sum

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Types 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

abstract type AbstractModelParameters{S, T} end 

struct SEIRParameters{S, T} <: AbstractModelParameters{S, T}
    β                       :: S
    γ                       :: T
    ν                       :: T 
    ψ                       :: T
end

SEIRParameters(; beta, gamma, nu, psi) = SEIRParameters(beta, gamma, nu, psi)

abstract type AbstractCompartments{T} <: AbstractVector{T} end

struct SEIRCompartments{T} <: AbstractCompartments{T}
    S                       :: T 
    E                       :: T 
    I                       :: T 
    R                       :: T 
    IdentifiedCases         :: T 
    Cumulative_I            :: T 
    Cumulative_Identified   :: T
end

function SEIRCompartments(args::Vararg{T}) where T
    u = zeros(T, 7)
    for (i, a) ∈ enumerate(args) u[i] = a end 
    return SEIRCompartments(u...)
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Constants 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

const REGIONCODEDF = DataFrame(
    :RegionCode => [ "UK_ENG", "UK_NIR", "UK_SCO", "UK_WAL" ],
    :RegionId => 1:4
)

# result from https://github.com/mrc-ide/EpiEstim/blob/master/data/covid_deaths_2020_uk.rda
COVIDSERIALINTERVAL = [  
    0.0000000000,
    0.0440204506,
    0.1298284450,
    0.1397552873,
    0.1277718301,
    0.1100166556,
    0.0917470443,
    0.0749977679,
    0.0604725660,
    0.0482765015,
    0.0382484935,
    0.0301228893,
    0.0236092441,
    0.0184305583,
    0.0143398489,
    0.0111254255,
    0.0086104507,
    0.0066498251,
    0.0051260438,
    0.0039448946,
    0.0030314300,
    0.0023264019,
    0.0017832132,
    0.0013653739,
    0.0010444113,
    0.0007981781,
    0.0006094926,
    0.0004650564,
    0.0003545982,
    0.0002701988,
    0.0002057625
]

const POPULATION2020 = [  
    # values from 
    # https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland/mid2020/ukpopestimatesmid2020on2021geography.xls
    56_550_138,  # ENGLAND
    1_895_510,  # NORTHERN IRELAND
    5_466_000,  # SCOTLAND
    3_169_586,  # WALES
]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SEIR model  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function seir_noisy(u, p, t)
    @unpack S, E, I, R = u
    N = +(S, E, I, R)
    @unpack ψ = p
    probabilities = _seirprobabilities(I, N, p, t)
    events = [ 
        (
            c = rand(Poisson(probabilities[i] * u[i]));  
            c < u[i] ? c : u[i]  # sampled value may be larger than compartment
        )
        for i ∈ 1:3
    ]
    recordedevents = rand(Binomial(events[2], ψ))
    return _seiroutputs(u, events, recordedevents)
end

function seir_deterministic(u, p, t)
    @unpack S, E, I, R = u
    N = +(S, E, I, R)
    @unpack ψ = p
    probabilities = _seirprobabilities(I, N, p, t)
    events = [ probabilities[i] * u[i] for i ∈ 1:3 ]
    recordedevents = events[2] * ψ
    return _seiroutputs(u, events, recordedevents)
end

function _seirprobabilities(I, N, p, t)
    β = _unpackbeta(p, t)
    @unpack γ, ν = p
    probabilities = [ 
        1 - exp(-β * I / N),  # S -> E 
        1 - exp(-ν),  # E -> I 
        1 - exp(-γ)  # I -> R
    ]
    return probabilities
end

function _seiroutputs(u, events, recordedevents)
    for v ∈ [ events; [ recordedevents ] ]
        if isnan(v) 
            return SEIRCompartments(
                zero(u[1] - events[1]),
                zero(u[2] + events[1] - events[2]),
                zero(u[3] + events[2] - events[3]),
                zero(u[4] + events[3]),
                zero(recordedevents),
                u[6],
                u[7]
            )
        end
    end

    return SEIRCompartments(
        u[1] - events[1],
        u[2] + events[1] - events[2],
        u[3] + events[2] - events[3],
        u[4] + events[3],
        recordedevents,
        u[6] + events[2],
        u[7] + recordedevents
    )
end

_unpackbeta(p::SEIRParameters{S, T}, t) where {S <: Function, T <: Any} = p.β(t)
_unpackbeta(p::SEIRParameters{S, T}, ::Any) where {S <: Number, T <: Any} = p.β
_unpackbeta(p::SEIRParameters{S, T}, t) where {S <: Vector, T <: Any} = p.β[t]

function runseir_noisy(u0::S, p, duration) where S <: AbstractVector{T} where T 
    compartments = Vector{S}(undef, duration)
    reportedcases = Vector{T}(undef, duration)
    u1 = seir_noisy(u0, p, 1)
    compartments[1] = u1 
    reportedcases[1] = u1[5]
    for t ∈ 2:duration 
        compartments[t] = seir_noisy(compartments[t-1], p, t)
        reportedcases[t] = compartments[t][5]
    end
    return @ntuple compartments reportedcases
end

function runseir_deterministic(u0::S, p, duration) where S <: AbstractVector{T} where T 
    # note that `T` must not be an Integer type
    compartments = Vector{S}(undef, duration)
    reportedcases = Vector{T}(undef, duration)
    u1 = seir_deterministic(u0, p, 1)
    compartments[1] = u1 
    reportedcases[1] = u1[5]
    for t ∈ 2:duration 
        compartments[t] = seir_deterministic(compartments[t-1], p, t)
        reportedcases[t] = compartments[t][5]
    end
    return @ntuple compartments reportedcases
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extras 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function getindex(a::SEIRCompartments{T}, i) where T 
    @boundscheck checkbounds(a, i) 
    v = [ a.S, a.E, a.I, a.R, a.IdentifiedCases, a.Cumulative_I, a.Cumulative_Identified ]
    return getindex(v, i)
end

length(::SEIRCompartments) = 7
size(::SEIRCompartments) = ( 7, )

sum(a::SEIRCompartments) = sum(a[1:4])
anyminimum(a::SEIRCompartments) = minimum([ a[i] for i ∈ 1:7 ])
    
end  # module AnalysisFunctions
