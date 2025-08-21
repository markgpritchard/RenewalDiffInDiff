
function analysisworkflow(model; priorsseed=nothing, sampleseed=nothing, kwargs...)
    priorsrng = _workflowrngs(priorsseed)
    samplerng = _workflowrngs(sampleseed)
    return analysisworkflow(priorsrng, samplerng, model; kwargs...)
end

function analysisworkflow(
    priorsrng::AbstractRNG, samplerng::AbstractRNG, model; 
    chain, name, nsamples, kwargs...
)
    return _analysisworkflow(priorsrng, samplerng, model, nsamples, name, chain; kwargs...)
end

function _analysisworkflow(
    priorsrng, samplerng, model, nsamples::Integer, name::AbstractString, chain::Integer; 
    kwargs...
)
    priorsdf, priorschain = priorsworkflow(
        priorsrng, model; 
        chain, name, kwargs...
    )
    map_df, map_estimate = maximumlikelihoodworkflow(
        model, priorsdf; 
        chain, name, kwargs...
    )
    mcmcdf, mcmcchain = mcmcworkflow(
        samplerng, model, map_estimate; 
        chain, name, nsamples, kwargs...
    )
    d = Dict(
        # `map_estimate` contains anonymous functions and leads to warnings and errors when 
        # saved and loaded
        "priorschain" => priorschain, 
        "priorsdf" => priorsdf,
        "mapdf" => map_df,
        "mcmcchain" => mcmcchain, 
        "mcmcdf" => mcmcdf,
    )
    safesave(datadir("sims", "$(name)_results_$(chain)_$(nsamples)samples.jld2"), d)
    return d
end

function priorsworkflow(model; priorsseed=nothing, kwargs...)
    priorsrng = _workflowrngs(priorsseed)
    return priorsworkflow(priorsrng, model; kwargs...)
end

function priorsworkflow(priorsrng::AbstractRNG, model; chain, name, npriors=10_000, kwargs...)
    return _priorsworkflow(priorsrng, model, chain, name, npriors; kwargs...)
end

function _priorsworkflow(
    priorsrng, model, chain::Integer, name::AbstractString, npriors::Integer; 
    savepriors=(chain == 1), kwargs...
)
    priorschain = sample(priorsrng, model, Prior(), npriors)
    priorsdf = DataFrame(priorschain)
  #  if savepriors
  #      safesave(
   #         datadir("sims", "$(name)_prior.jld2"), 
   #         Dict("priorschain" => priorschain, "priorsdf" => priorsdf)
   #     )
   # end
    return (priorsdf, priorschain)
end

function maximumlikelihoodworkflow(
    model, priorsdf=nothing; 
    chain, name, mapmaxtime=600, kwargs...
)
    return _maximumlikelihoodworkflow(model, priorsdf, chain, name, mapmaxtime)
end

function _maximumlikelihoodworkflow(
    model, priorsdf::DataFrame, chain::Integer, name::AbstractString, mapmaxtime::Integer
)
    indexformap = findall(x -> x == chain, ordinalrank(priorsdf.lp; rev=true))[1]
    initparamslastindex = size(priorsdf, 2) - 3
    initparamsformap = [values(priorsdf[indexformap, 3:initparamslastindex])...]
    map_estimate = maximum_likelihood(
        model; 
        adtype=AutoReverseDiff(), initial_params=initparamsformap, maxtime=mapmaxtime,
    )
    return __maximumlikelihoodworkflow(map_estimate, chain, name)
end

function _maximumlikelihoodworkflow(
    model, ::Nothing, chain::Integer, name::AbstractString, mapmaxtime::Integer
)
    map_estimate = maximum_likelihood(
        model; 
        adtype=AutoReverseDiff(), maxtime=mapmaxtime,
    )
    return __maximumlikelihoodworkflow(map_estimate, chain, name)
end

function __maximumlikelihoodworkflow(map_estimate, chain, name)
    map_df = map_DataFrame(map_estimate)
   # safesave(
   #     datadir("sims", "$(name)_map_$chain.jld2"), 
   #     Dict("mapdf" => map_df)
   # )
    return (map_df, map_estimate)
end

function mcmcworkflow(model, map_estimate=nothing; sampleseed=nothing, kwargs...)
    samplerng = _workflowrngs(sampleseed)
    return mcmcworkflow(samplerng, model, map_estimate; kwargs...)
end

function mcmcworkflow(
    samplerng::AbstractRNG, model, map_estimate=nothing; 
    acceptancedelta=0.65, chain, name, nsamples, kwargs...
)
    return _mcmcworkflow(
        samplerng, model, map_estimate, nsamples, acceptancedelta, chain, name
    )
end

function _mcmcworkflow(
    samplerng, 
    model, 
    map_estimate, 
    nsamples::Integer, 
    acceptancedelta::Float64, 
    chain::Integer, 
    name::AbstractString
)
    mcmcchain = sample(
        samplerng, model, NUTS(acceptancedelta; adtype=AutoReverseDiff()), nsamples; 
        initial_params=map_estimate.values.array
    ) 
    return __mcmcworkflow(mcmcchain, nsamples, chain, name)
end

function _mcmcworkflow(
    samplerng,
    model,
    ::Nothing, 
    nsamples::Integer, 
    acceptancedelta::Float64, 
    chain::Integer, 
    name::AbstractString
)
    mcmcchain = sample(
        samplerng, model, NUTS(acceptancedelta; adtype=AutoReverseDiff()), nsamples; 
    ) 
    return __mcmcworkflow(mcmcchain, nsamples, chain, name)
end

function __mcmcworkflow(mcmcchain, nsamples, chain, name)
    mcmcdf = DataFrame(mcmcchain)
   # safesave(
   #     datadir("sims", "$(name)_mcmc_$(chain)_$(nsamples)samples.jld2"), 
   #     Dict("mcmcchain" => mcmcchain, "mcmcdf" => mcmcdf)
   # )
    return (mcmcdf, mcmcchain)
end

_workflowrngs(::Nothing) = default_rng()
_workflowrngs(x::Integer) = Xoshiro(x)
_workflowrngs(rng::AbstractRNG) = rng
