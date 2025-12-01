
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
    map_df, map_estimate = maximumaposterioriworkflow(
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

function priorsworkflow(
    priorsrng::AbstractRNG, model; 
    chain, name, npriors=10_000, kwargs...
)
    return _priorsworkflow(priorsrng, model, chain, name, npriors; kwargs...)
end

function _priorsworkflow(
    priorsrng, model, chain::Integer, name::AbstractString, npriors::Integer; 
    savepriors=(chain == 1), kwargs...
)
    priorschain = sample(priorsrng, model, Prior(), npriors)
    priorsdf = DataFrame(priorschain)
    return (priorsdf, priorschain)
end

function indexesformap(df, chain::Integer)
    return findall(
        x -> x == chain, ordinalrank([isnan(x) ? -Inf : x for x in df.loglikelihood]; 
        rev=true)
    )
end

function indexesformap(df, chain::AbstractVector)
    return findall(
        x -> x ∈ chain, ordinalrank([isnan(x) ? -Inf : x for x in df.loglikelihood]; 
        rev=true)
    )
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
    indexformap = indexesformap(priorsdf, chain)[1]
    initparamslastindex = size(priorsdf, 2) - 3
    initparamsformap = [values(priorsdf[indexformap, 3:initparamslastindex])...]
    map_estimate = maximum_likelihood(
        model; 
        adtype=AutoMooncake(), 
        initial_params=initparamsformap, 
        maxtime=mapmaxtime,
    )
    return __maximumlikelihoodworkflow(map_estimate, chain, name)
end

function _maximumlikelihoodworkflow(
    model, ::Nothing, chain::Integer, name::AbstractString, mapmaxtime::Integer
)
    map_estimate = maximum_likelihood(
        model; 
        adtype=AutoMooncake(), maxtime=mapmaxtime,
    )
    return __maximumlikelihoodworkflow(map_estimate, chain, name)
end

function __maximumlikelihoodworkflow(map_estimate, chain, name)
    map_df = map_DataFrame(map_estimate)
    return (map_df, map_estimate)
end

function maximumaposterioriworkflow(
    model, priorsdf=nothing; 
    chain, name, mapmaxtime=600, kwargs...
)
    return _maximumaposterioriworkflow(model, priorsdf, chain, name, mapmaxtime)
end

function _maximumaposterioriworkflow(
    model, priorsdf::DataFrame, chain::Integer, name::AbstractString, mapmaxtime::Integer
)
    indexformap = indexesformap(priorsdf, chain)[1]
    initparamslastindex = size(priorsdf, 2) - 3
    initparamsformap = [values(priorsdf[indexformap, 3:initparamslastindex])...]
    map_estimate = maximum_a_posteriori(
        model; 
        adtype=AutoMooncake(), 
        initial_params=initparamsformap, 
        maxtime=mapmaxtime,
    )
    return __maximumaposterioriworkflow(map_estimate, chain, name)
end

function _maximumaposterioriworkflow(
    model, ::Nothing, chain::Integer, name::AbstractString, mapmaxtime::Integer
)
    map_estimate = maximum_a_posteriori(
        model; 
        adtype=AutoMooncake(), maxtime=mapmaxtime,
    )
    return __maximumaposterioriworkflow(map_estimate, chain, name)
end

function __maximumaposterioriworkflow(map_estimate, chain, name)
    map_df = map_DataFrame(map_estimate)
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
        samplerng, 
        model, 
        NUTS(acceptancedelta; adtype=AutoMooncake()), 
        nsamples; 
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
        samplerng, 
        model, 
        NUTS(acceptancedelta; adtype=AutoMooncake()), 
        nsamples; 
    ) 
    return __mcmcworkflow(mcmcchain, nsamples, chain, name)
end

function __mcmcworkflow(mcmcchain, nsamples, chain, name)
    mcmcdf = DataFrame(mcmcchain)
    for i in axes(mcmcdf, 1)
        mcmcdf.chain[i] = chain 
    end
    return (mcmcdf, mcmcchain)
end

_workflowrngs(::Nothing) = default_rng()
_workflowrngs(x::Integer) = Xoshiro(x)
_workflowrngs(rng::AbstractRNG) = rng

## load saved outputs 

function loadsamples(id; analysisname="analysis$(id)", kwargs...)
    return _loadsamples(id, analysisname; kwargs...)
end

function _loadsamples(
    id, analysisname; 
    data=load(simulationdir("sim$(id).jld2"))["sim"], 
    nchains=8, 
    nsamples=[25, 100, 1000, 2000],
    initchain=1,
)
    initnsamples = _initnsamples(analysisname, nsamples; chain=initchain)
    initnsamples == 0 && return nothing 
    samples1 = load(
        datadir("sims", "$(analysisname)_results_$(initchain)_$(initnsamples)samples.jld2")
    )
    priorsdf = samples1["priorsdf"]
    mapdf = samples1["mapdf"]
    mcmcdf = samples1["mcmcdf"]
    for i in 2:nchains 
        if isfile(datadir("sims", "$(analysisname)_results_$(i)_$(initnsamples)samples.jld2"))
            additionalrows = load(
                datadir("sims", "$(analysisname)_results_$(i)_$(initnsamples)samples.jld2")
            )["mcmcdf"] 
            mcmcdf = vcat(mcmcdf, additionalrows)
        end
    end 
    return @ntuple priorsdf mapdf mcmcdf data
end

function _initnsamples(analysisname, nsamples; chain=1)
    reversensamples = sort(nsamples; rev=true)
    initnsamples = 0 
    for s in reversensamples
        initnsamples > 0 && continue
        if isfile(datadir("sims", "$(analysisname)_results_$(chain)_$(s)samples.jld2"))
            initnsamples += s
        end
    end
    return initnsamples
end
