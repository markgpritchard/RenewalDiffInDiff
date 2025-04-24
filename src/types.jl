

struct Automatic end  # not exported
const automatic = Automatic()  # not exported

abstract type AbstractInterventionsArray{T, N} <: AbstractArray{T, N} end

struct InterventionsMatrix{T} <: AbstractInterventionsArray{T, 2}
    starttimes  :: Vector{Int}
    duration    :: Int
end

function InterventionsMatrix(T, starttimes::Vector{<:Union{Int, Nothing}}, duration) 
    # set any `missing` start times to be after the end of the study duration
    newstarttimes = Vector{Int}(undef, length(starttimes)) 
    for (i, t) ∈ enumerate(starttimes) 
        newstarttimes[i] = _setinterventionsmatrixstarttimes(t, duration)
    end
    return InterventionsMatrix{T}(newstarttimes, duration) 
end

function InterventionsMatrix(starttimes, duration)
    InterventionsMatrix(Int, starttimes, duration)
end

_setinterventionsmatrixstarttimes(t::Int, ::Any) = t 
_setinterventionsmatrixstarttimes(::Nothing, duration) = duration + 1



struct InterventionsVector{T} <: AbstractInterventionsArray{T, 1}
    starttime   :: Int
    duration    :: Int
end

function InterventionsVector{T}(::Nothing, duration) where T
    return InterventionsVector{T}(duration + 1, duration) 
end

InterventionsVector(args...) = InterventionsVector{Int}(args...)

abstract type AbstractModelParameters{S, T, U} end 

struct SEIRParameters{S, T, U} <: AbstractModelParameters{S, T, U} 
    β           :: S 
    μ           :: T 
    γ           :: T 
    θ           :: U

    function SEIRParameters(β, μ, γ, θ) 
        # The generation interval calculations used in this analysis assume that μ != γ
        @assert μ != γ
        S = typeof(β)
        T = typeof(μ)
        U = typeof(θ)
        return new{S, T, U}(β, μ, γ, θ) 
    end
end

SEIRParameters(; beta, gamma, nu, psi) = SEIRParameters(beta, gamma, nu, psi)
