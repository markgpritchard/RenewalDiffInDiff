
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
