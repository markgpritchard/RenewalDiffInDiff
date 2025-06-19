
module RenewalDiffInDiff

#using DataFrames
using Distributions
using DrWatson
#using ForwardDiff
#using NaNMath
using PrettyTables: pretty_table
using Random
using Turing
using StatsBase: sample, Weights 
import Random: default_rng

include("consts.jl")
include("interventionsarrays.jl")
include("generationinterval.jl")
include("renewalequation.jl")
include("extras.jl")

include("simulations.jl")

export
    ## types.jl
    SEIRParameters,
    ## consts.jl
    COVIDSERIALINTERVAL, 
    POPULATION2020,
    ## interventionsarrays.jl
    InterventionsMatrix, 
    InterventionsVector, 
    duration,
    offsetinterventionsmatrix,
    ## generationinterval.jl
    covidvectorg,
    vectorg,
    ## renewalequation.jl
    expectedinfections,
    initialsusceptible,
    proportionwaned,
    R0_did,
    ## extras.jl
    interventionsoffset, 
    seir_deterministic,
    ## simulations.jl
    SEIRParameters
    
end  # module RenewalDiffInDiff
