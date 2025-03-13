
module RenewalDiffInDiff

using DrWatson
using CubicSplines, DataFrames, Distributions, ForwardDiff, NaNMath, Random, Turing
import Base: getindex, length, size, sum

include("types.jl")
include("consts.jl")
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
    ## renewalequation.jl
    testf,
    ## parameterfitting.jl
    diffindiffparameters, 
    diffindiffparameters_discretetimes, 
    diffindiffparameters_fittocurve_splinetimes, 
    diffindiffparameters_polytimes, 
    diffindiffparameters_twodiscretetimes, 
    diffindiffparameters_splinetimes, 
    generatew_gt, 
    generatew_gtrow, 
    generatez_gtminus1,
    keyvalues, 
    loadanalysisdictsasdf,
    ## processoutputs.jl
    insertcumulativeeffects!,
    logbasicreproductionratios,
    logeffectivereproductionratios,
    predictcases,
    quantilelogeffectivereproductionratios,
    quantilepredictcases,
    quantilepredictcumulativedifferenceincases,
    quantilepredictdifferenceincases,
    tauquantiles,
    ## extras.jl
    getindex, 
    interventionsoffset, 
    lagleadinterventionsmatrix, 
    length, 
    seir_deterministic, 
    size
    
end  # module RenewalDiffInDiff
