
module RenewalDiffInDiff

using DrWatson
using CubicSplines, DataFrames, Distributions, ForwardDiff, NaNMath, Random, Turing
import Base: getindex, length, size, sum

include("types.jl")
include("consts.jl")
include("renewalequation.jl")
include("parameterfitting.jl")
include("extras.jl")

 
## types.jl
export InterventionsMatrix, InterventionsVector, SEIRParameters
## consts.jl/
export COVIDSERIALINTERVAL, POPULATION2020

## renewalequation.jl
#export calculatesumfi, expectedinfections, poissoninfections, runrenewalequation, 
#    runrenewalequation!, runrenewalequationsample!, runrenewalequationsamples, testf
export renewalequation_expectedcases, renewalequation_expectedcases!,
    renewalequation_expectedcases_rt, renewalequation_expectedcases_rt!, 
    renewalequation_poissoncases, renewalequation_poissoncases!, 
    renewalequation_poissoncases_rt, renewalequation_poissoncases_rt!, 
    samplerenewalequation, samplerenewalequation!, samplerenewalequation_2sets, 
    samplerenewalequation_counterfactual, samplerenewalequation_polytimes_2sets,
    samplerhomatrix, samplerhomatrix!, testf

## parameterfitting.jl
export diffindiffparameters, diffindiffparameters_discretetimes, 
    diffindiffparameters_fittocurve_splinetimes, diffindiffparameters_polytimes, 
    diffindiffparameters_twodiscretetimes, diffindiffparameters_splinetimes, generatew_gt, 
    generatew_gtrow, generatez_gtminus1, keyvalues, loadanalysisdictsasdf

## extras.jl
export getindex, interventionsoffset, length, seir_deterministic, size
    
end  # module RenewalDiffInDiff
