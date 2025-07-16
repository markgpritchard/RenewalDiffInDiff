
function dataforstanmodel(; g, incidence, interventions, nseedtimes, Ns, glength=automatic)
    return _dataforstanmodel(g, incidence, interventions, nseedtimes, Ns, glength)
end

function _dataforstanmodel(
    ::Function, 
    ::Matrix{<:Integer}, 
    ::Any, 
    ::Integer, 
    ::Vector{<:Integer}, 
    ::RenewalDiffInDiff.Automatic,
)
    throw(ArgumentError("If `g` is provided as a function, `glength` must be an integer"))
end

function _dataforstanmodel(
    g::Function, 
    incidence::Matrix{<:Integer}, 
    interventions, 
    nseedtimes::Integer, 
    Ns::Vector{<:Integer}, 
    glength::Integer,
)
    _g = [ g(t) for t ∈ 1:glength ]
    return _dataforstanmodel(_g, incidence, interventions, nseedtimes, Ns, glength,)
end

function _dataforstanmodel(
    g::AbstractVector, 
    incidence::Matrix{<:Integer}, 
    interventions, 
    nseedtimes::Integer, 
    Ns::Vector{<:Integer}, 
    ::RenewalDiffInDiff.Automatic,
)
    return _dataforstanmodel(g, incidence, interventions, nseedtimes, Ns, length(g))
end

function _dataforstanmodel(
    g::AbstractVector, 
    incidence::Matrix{<:Integer}, 
    interventions::InterventionsMatrix{<:Integer}, 
    nseedtimes::Integer, 
    Ns::Vector{<:Integer}, 
    glength::Integer,
)
    _interventions = collect(interventions)
    return _dataforstanmodel(g, incidence, _interventions, nseedtimes, Ns, glength)
end

function _dataforstanmodel(
    g::AbstractVector, 
    incidence::Matrix{<:Integer}, 
    interventions::Matrix{<:Integer}, 
    nseedtimes::Integer, 
    Ns::Vector{<:Integer}, 
    glength::Integer,
)
    length(g) == glength || throw(_glengtherror(g, glength))
    size(incidence) == size(interventions) || throw(_matrixsizeerr(incidence, interventions))

    return Dict(
        "ntimes" => size(incidence, 1),
        "nlocations" => size(incidence, 2),
        "nseedtimes" => nseedtimes,
        "glength" => glength,
        "incidence" => incidence,
        "interventions" => interventions,
        "Ns" => Ns,
        "g" => g,
    )
end

function _glengtherror(g, glength)
    return DimensionMismatch(
        "$glength != $(length(g)), `glength` must equal the length of vector `g`"
    )
end

function _matrixsizeerr(incidence, interventions)
    sic = size(incidence)
    siv = size(interventions)
    return DimensionMismatch(
        "$sic != $siv, `incidence` and `interventions` must have equal sizes"
    )
end

