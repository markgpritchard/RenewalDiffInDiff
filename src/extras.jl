
function interventionsoffset(A::InterventionsMatrix{T}, v::Int) where T
    newstime = [ _interventionsoffset(stime, v, A.duration) for stime ∈ A.starttimes ]
    return InterventionsMatrix(T, newstime, A.duration)
end

function interventionsoffset(A::Matrix{T}, v::Int) where T
    tmax = size(A, 1)
    outputmatrix = zeros(T, size(A))
    for g ∈ axes(outputmatrix, 2)
        for t ∈ axes(outputmatrix, 1)
            if t - v < 1 
                outputmatrix[t, g] = zero(T)
            elseif t - v > tmax 
                outputmatrix[t, g] = A[tmax, g]
            else 
                outputmatrix[t, g] = A[(t - v), g]
            end
        end
    end
    return outputmatrix
end

_interventionsoffset(::Nothing, ::Any, dur) = dur + 1

function _interventionsoffset(stime::Int, v, dur)
    if stime > dur 
        return dur + 1 
    elseif stime + v < 1 
        return 1 
    else
        return stime + v 
    end
end

function getindex(A::InterventionsMatrix{T}, ::Colon, i) where T 
    return InterventionsVector{T}(A.starttimes[i], A.duration)
end

function getindex(A::InterventionsMatrix{T}, i, j) where T 
    @boundscheck checkbounds(A, i, j) 
    if i < A.starttimes[j] 
        return zero(T)
    else
        return one(T)
    end
end

function getindex(v::InterventionsVector{T}, i) where T 
    @boundscheck checkbounds(v, i) 
    if i < v.starttime 
        return zero(T)
    else
        return one(T)
    end
end

length(A::InterventionsMatrix) = length(A.starttimes) * A.duration
length(v::InterventionsVector) = v.duration

function size(A::InterventionsMatrix)
    l = A.duration
    w = length(A.starttimes)
    return ( l, w )
end

size(v::InterventionsVector) = ( v.duration, )

function _nunique(x)
    y = unique(x)
    return length(y)
end

function seir_deterministic(u, p, ::Number)
    s, e, i, r = u
    n = sum(u)
    return [
        s * (exp(-p.β * i / n)),  # s 
        s * (1 - exp(-p.β * i / n)) + e * (exp(-p.μ)),  # e 
        e * (1 - exp(-p.μ)) + i * (exp(-p.γ)),  # i 
        r + i * (1 - exp(-p.γ))  # r
    ]
end

function seir_deterministic(u0, p, tvec::AbstractVector{<:Number})
    outputs = zeros(Float64, length(tvec), 4)
    u = deepcopy(u0)
    outputs[1, :] = u 
    for (i, t) ∈ enumerate(tvec)
        u = seir_deterministic(u, p, t)
        outputs[i, :] = u 
    end 
    return outputs 
end
