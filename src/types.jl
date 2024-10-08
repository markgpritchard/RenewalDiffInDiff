

struct Automatic end  # not exported
const automatic = Automatic()  # not exported

abstract type AbstractInterventionsArray{T, N} <: AbstractArray{T, N} end

struct InterventionsMatrix{T} <: AbstractInterventionsArray{T, 2}
    starttimes  :: Vector{Int}
    duration    :: Int
end

function InterventionsMatrix(T::DataType, starttimes::Vector{Int}, duration) 
    return InterventionsMatrix{T}(starttimes, duration)
end

function InterventionsMatrix(
    T::DataType, starttimes::Vector{U}, duration
) where U <: Union{Int, Nothing} 
    newstarttimes = Vector{Int}(undef, length(starttimes)) 
    for (i, t) ∈ enumerate(starttimes) 
        if isnothing(t)
            newstarttimes[i] = duration + 1 
        else 
            newstarttimes[i] = t 
        end
    end
    return InterventionsMatrix(T, newstarttimes, duration) 
end

InterventionsMatrix(args...) = InterventionsMatrix(Int, args...)

struct InterventionsVector{T} <: AbstractInterventionsArray{T, 1}
    starttime   :: Int
    duration    :: Int
end

function InterventionsVector(T::DataType, starttime::Int, duration) 
    return InterventionsVector{T}(starttime, duration)
end

function InterventionsVector(T::DataType, ::Nothing, duration) 
    return InterventionsVector(T, duration + 1, duration) 
end

InterventionsVector(args...) = InterventionsVector(Int, args...)

abstract type AbstractModelParameters{S, T, U} end 

struct SEIRParameters{S, T, U} <: AbstractModelParameters{S, T, U} 
    β           :: S 
    μ           :: T 
    γ           :: T 
    θ           :: U
end

SEIRParameters(; beta, gamma, nu, psi) = SEIRParameters(beta, gamma, nu, psi)
