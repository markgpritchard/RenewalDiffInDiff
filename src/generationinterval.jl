
function vectorg(t::Integer, vec::AbstractVector{<:T}; max_t=length(vec)) where T 
    t <= 0 && return zero(T)
    t > max_t && return zero(T)
    return getindex(vec, t)
end

covidvectorg(t) = vectorg(t, COVIDSERIALINTERVAL)
