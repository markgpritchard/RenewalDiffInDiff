
module RenewalDiffInDiff

#using DataFrames
using Distributions
using DrWatson
#using ForwardDiff
#using NaNMath
using PrettyTables: pretty_table
using Random
#using Turing
using StatsBase: sample, Weights 
import Random: default_rng

include("consts.jl")
include("interventionsarrays.jl")
include("generationinterval.jl")
include("renewalequation.jl")
include("parameterfitting.jl")
include("extras.jl")

include("simulations.jl")

export
    ## types.jl
    SEIRParameters,
    ## consts.jl
    POPULATION2020,
    ## interventionsarrays.jl
    InterventionsMatrix, 
    InterventionsVector, 
    duration,
    offsetinterventionsmatrix,
    ## generationinterval.jl
    generationproportion,
    gseir,
    g_covid,
    ## renewalequation.jl
    expectedinfections,
    initialsusceptible,
    proportionwaned,
    R0_did,
    ## parameterfitting.jl
    dataforstanmodel,
    ## extras.jl
    interventionsoffset, 
    seir_deterministic,
    ## simulations.jl
    SEIRParameters
    
end  # module RenewalDiffInDiff
