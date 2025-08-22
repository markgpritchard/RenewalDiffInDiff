# functions for this analysis 

module RenewalDiffInDiff

import ReverseDiff

using CairoMakie: Axis, Colorbar, Label, Legend, TopLeft, hidespines!
using Distributions: Binomial, Distribution, Uniform
using DrWatson: @ntuple, datadir, load, safesave
using Random: AbstractRNG, Xoshiro, default_rng
using RenewalDiD: DataFrame, InterventionArray, InterventionMatrix
using RenewalDiD: map_DataFrame, simulationu0
using StatsBase: ordinalrank, sample
using Turing: AutoReverseDiff, MCMCThreads, NUTS, Prior, maximum_likelihood

export simulationdir
# simulationhelperfunctions.jl
export largepop, simu0, smallpop
# workflow.jl
export analysisworkflow, maximumlikelihoodworkflow, mcmcworkflow, priorsworkflow
export loadsamples 
# plotformatting.jl
export formataxis!, labelplots!, setvalue!
# setoffsets.jl
export addoffsetstointerventionarray

simulationdir(args...) = datadir("simulations", args...)

include("simulationhelperfunctions.jl")
include("workflow.jl")
include("plotformatting.jl")
include("setoffsets.jl")

end  # module RenewalDiffInDiff
