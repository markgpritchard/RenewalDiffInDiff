
module RenewalDiffInDiff

using DrWatson
using CubicSplines, DataFrames, Distributions, ForwardDiff, Memoize, NaNMath, Random, Turing
using StatsBase: sample, Weights 
import Base: getindex, length, size, sum
import Random: default_rng

include("types.jl")
include("consts.jl")
include("stochasticsimulation.jl")
include("renewalequation.jl")
include("parameterfitting.jl")
include("processoutputs.jl")
include("extras.jl")

export 
    ## types.jl
    InterventionsMatrix, 
    InterventionsVector, 
    SEIRParameters,
    ## consts.jl
    COVIDSERIALINTERVAL, 
    POPULATION2020,
    ## stochasticsimulation.jl
    stochasticiteration, 
    stochasticmodel,
    ## renewalequation.jl
    testf,
    ## parameterfitting.jl
    generatew_gt, 
    keyvalues, 
    loadanalysisdictsasdf,
    renewaldiffindiffparameters,
    ## processoutputs.jl
    insertcumulativeeffects!,
    logbasicreproductionratios,
    logeffectivereproductionratios,
    phiquantiles,
    predictcases,
    quantilelogeffectivereproductionratios,
    quantilepredictcases,
    quantilepredictcumulativedifferenceincases,
    quantilepredictdifferenceincases,
    ## extras.jl
    getindex, 
    interventionsoffset, 
    lagleadinterventionsmatrix, 
    length, 
    seir_deterministic, 
    size
    
end  # module RenewalDiffInDiff
