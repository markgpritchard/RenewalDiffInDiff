
abstract type AbstractInterventionsArray{T, N} <: AbstractArray{T, N} end

struct InterventionsMatrix{T} <: AbstractInterventionsArray{T, 2}
    starttimes::Vector{Int}
    duration::Int

    function InterventionsMatrix{T}(ts, d; warnalltreated=true) where T
        length(ts) >= 2 || throw(_starttimeargumenterror())
        duration = convert(Int, d)
        duration >= 2 || throw(_durationdomainerror(duration))
        starttimes = _convertstarttimes(ts, duration)
        if warnalltreated && maximum(starttimes) <= duration 
            _untreatedwarning(maximum(starttimes), duration)
        end
        return new{T}(starttimes, duration)
    end
end

function InterventionsMatrix(T::DataType, starttimes, duration; kwargs...) 
    return InterventionsMatrix{T}(starttimes, duration; kwargs...)
end

function InterventionsMatrix(starttimes, duration; kwargs...)
    return InterventionsMatrix{Int}(starttimes, duration; kwargs...)
end 

struct InterventionsVector{T} <: AbstractInterventionsArray{T, 1}
    starttime::Int
    duration::Int

    function InterventionsVector{T}(t, d) where T
        duration = convert(Int, d)
        duration >= 2 || throw(_durationdomainerror(duration)) 
        starttime = _convertstarttime(t, duration)
        return new{T}(starttime, duration)
    end
end

function InterventionsVector(T::DataType, starttime, duration) 
    return InterventionsVector{T}(starttime, duration)
end

InterventionsVector(starttime, duration) = InterventionsVector{Int}(starttime, duration)

_convertstarttimes(ts::AbstractVector, d) = [ _convertstarttime(t, d) for t ∈ ts ]
_convertstarttimes(ts::Tuple, d) = _convertstarttimes([ ts... ], d)

function _convertstarttime(t, ::Any)
    t > 1 || throw(_starttimedomainerror(t))
    return convert(Int, t)
end 

_convertstarttime(::Nothing, d) = d + 1

_durationdomainerror(d) = DomainError(d, "duration must be ≥2")

function _starttimeargumenterror()
    return ArgumentError(
        """
        at least two start times must be provided to `InterventionsMatrix` (at least two \
        groups must be represented). To represent a single group use `InterventionsVector`
        """
    )
end

function _starttimedomainerror(t)
    return DomainError(t, "start time must be >1 (all groups must initially be untreated)")
end

function _untreatedwarning(maxst, d)
    if maxst == d 
        @warn "no untreated groups at time $maxst"
    else
        @warn "no untreated groups between times $maxst and $d"
    end
    return nothing
end

duration(A::AbstractInterventionsArray) = A.duration

Base.size(A::InterventionsMatrix) = ( duration(A), length(A.starttimes) )
Base.size(v::InterventionsVector) = ( duration(v), )

function Base.getindex(A::AbstractInterventionsArray, I...)
    @boundscheck checkbounds(A, I...) 
    return _getindex(A, to_indices(A, I)...)
end

function _getindex(A::AbstractInterventionsArray{T, N}, i, j=nothing) where {T, N} 
    if i < _jthstarttime(A, j)
        return zero(T)
    else
        return one(T)
    end
end

_jthstarttime(A::InterventionsMatrix, j) = A.starttimes[j] 
_jthstarttime(v::InterventionsVector, ::Any) = v.starttime

function Base.getindex(A::InterventionsMatrix{T}, ::Colon, i) where T 
    return InterventionsVector{T}(A.starttimes[i], duration(A))
end

function Base.show(io::IO, ::MIME"text/plain", A::AbstractInterventionsArray) 
    data = _showinterventionsarraydata(A)
    return pretty_table(
        io, data; 
        header=_showinterventionsarrayheader(A),
        hlines=_showinterventionsarrayhlines(A), 
        show_row_number=false, 
        title=summary(A), 
        vlines=[ 1 ], 
    )
end

function _showinterventionsarraydata(A) 
    allt = _showinterventionsarraytimes(A)
    data = _showinterventionsarraydatarow(A, allt, 1)
    
    for i ∈ eachindex(allt)
        i == 1 && continue 
        if allt[i] != allt[i-1] + 1 
            data = vcat(data, _showinterventionsarraygaprow(A))
        end
        data = vcat(data, _showinterventionsarraydatarow(A, allt, i)) 
    end

    return data
end

function _showinterventionsarraytimes(A)
    allt = _showinterventionsarrayalltimes(A)
    filter!(x -> x <= duration(A), allt)
    uniquet = unique(allt)
    sort!(uniquet)
    return uniquet
end

_showinterventionsarrayalltimes(A::InterventionsMatrix) = [ A.starttimes; 1; duration(A) ]
_showinterventionsarrayalltimes(v::InterventionsVector) = [ v.starttime, 1, duration(v) ]

function _showinterventionsarraydatarow(A, allt, i)
    datarow = Matrix{String}(undef, 1, size(A, 2) + 1)
    datarow[1, 1] = "$(allt[i])"

    for c ∈ axes(A, 2)
        datarow[1, c+1] = "$(getindex(A, allt[i], c))"
    end

    return datarow 
end

_showinterventionsarraygaprow(A) = [ "⋮" for _ ∈ 1:1, _ ∈ 1:(size(A, 2) + 1) ]

function _showinterventionsarrayheader(A::InterventionsMatrix)
    return [ "time"; [ "$x" for x ∈ axes(A, 2) ] ]
end

_showinterventionsarrayheader(::InterventionsVector) = [ "time", " " ]

_showinterventionsarrayhlines(::InterventionsMatrix) = [ 1 ]
_showinterventionsarrayhlines(::InterventionsVector) = Int[ ]

function offsetinterventionsmatrix(M::InterventionsMatrix{T}, offset) where T 
    starttimes = [ _offsetinterventiontime(t, duration(M), offset) for t ∈ M.starttimes ]
    return InterventionsMatrix{T}(starttimes, duration(M))
end

function _offsetinterventiontime(t, duration, offset)
    if t > duration || t + offset <= 1 || t + offset > duration 
        return duration + 1 
    else
        return t + convert(Int, offset) 
    end
end
