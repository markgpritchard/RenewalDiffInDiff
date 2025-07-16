
# result from https://github.com/mrc-ide/EpiEstim/blob/master/data/covid_deaths_2020_uk.rda
const COVIDSERIALINTERVAL = [  # not exported
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

function generationproportion(t::Integer, g; max_t=automatic, kwargs...)
    # `kwargs` are passed to `g` if `g` is a function; if `g` is not a function then 
    # additional keyword arguments will lead to an error
    return _generationproportion(t, g, max_t, automatic, automatic; kwargs...)
    # the two `automatic` arguments are used to indicate tests with respect to `max_t`; they 
    # are removed after tests are performed or if not needed
end

function _generationproportion(t, g::AbstractVector, ::Automatic, ::Any, ::Any)
    return _generationproportion(t, g, length(vec), nothing, automatic)  
end

function _generationproportion(t, g::Function, ::Automatic, ::Any, ::Any; kwargs...)
    return _generationproportion(t, g; kwargs...)
end

function _generationproportion(t, g::AbstractVector, max_t::Integer, ::Automatic, ::Any)
    max_t <= length(g) || throw(_maxlimit_generationproportionerror(g, max_t))
    return _generationproportion(t, g, length(vec), nothing, automatic)
end

function _generationproportion(t, g::Function, max_t::Integer, ::Automatic, ::Any; kwargs...)
    return _generationproportion(t, g, max_t, nothing, automatic; kwargs...)
end

function _generationproportion(t, g, max_t::Integer, ::Nothing, ::Automatic; kwargs...)  
    # first `nothing` in signature when there is no need to check `max_t <= length(g)`
    max_t >= 1 || throw(ArgumentError("$max_t, `max_t` must be at least 1"))
    return _generationproportion(t, g, max_t, nothing, nothing; kwargs...)  
end

function _generationproportion(t, g, max_t::Integer, ::Nothing, ::Nothing; kwargs...)  
    # second `nothing` in signature when there is no need to check `max_t >= 1`
    if t > max_t 
        return _zeroproportion(g; kwargs...)
    else
        return _generationproportion(t, g; kwargs...) 
    end
end

function _generationproportion(t, g; kwargs...) 
    if t <= 0
        return _zeroproportion(g; kwargs...)
    else
        return __generationproportion(t, g; kwargs...) 
    end
end

__generationproportion(t, g::AbstractVector) = getindex(g, t)
__generationproportion(t, g::Function; kwargs...) = g(t; kwargs...)
_zeroproportion(g; kwargs...) = zero(__generationproportion(1, g; kwargs...))

function _maxlimit_generationproportionerror(g, max_t)
    return ArgumentError(
        "$max_t, `max_t` ($max_t) cannot be longer than the vector `g` (length $(length(g))"
    )
end

g_covid(t::Integer) = _generationproportion(t, COVIDSERIALINTERVAL, 31, nothing, nothing)

gseir(t; gamma=automatic, sigma=automatic) = _gseir(t, sigma, gamma)
gseir(t, sigma, gamma) = _gseir(t, sigma, gamma)
_gseir(t, ::Automatic, ::Automatic) = _gseirunequalgammasigma(t, 0.5, 0.4)  # default values

function _gseir(t, sigma, gamma)
    sigma > 0 || throw(_gseirzeroerror("sigma", sigma))
    gamma > 0 || throw(_gseirzeroerror("gamma", gamma))
    if gamma == sigma 
        return _gseirequalgammasigma(t, gamma)
    else
        return _gseirunequalgammasigma(t, sigma, gamma)
    end
end

_gseirequalgammasigma(t, gamma) = gamma^2 * t * exp(-gamma * t)

function _gseirunequalgammasigma(t, sigma, gamma)
    return sigma * gamma * (exp(-gamma * t) - exp(-sigma * t)) / (sigma - gamma)
end

_gseirzeroerror(name, v) = ArgumentError("$v, `gseir` requires $name to be positive")
