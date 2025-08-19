# functions for this analysis 

module RenewalDiffInDiff

import ReverseDiff

using CairoMakie: Axis, Colorbar, Label, Legend, TopLeft, hidespines!
using DrWatson: datadir, safesave
using Random: AbstractRNG, Xoshiro, default_rng
using RenewalDiD: DataFrame, map_DataFrame, simulationu0
using StatsBase: ordinalrank, sample
using Turing: AutoReverseDiff, Binomial, Distribution, MCMCThreads, NUTS, Prior, Uniform
using Turing: maximum_likelihood

# simulationhelperfunctions.jl
export largepop, simu0, smallpop
# workflow.jl
export analysisworkflow, maximumlikelihoodworkflow, mcmcworkflow, priorsworkflow
# plotformatting.jl
export formataxis!, labelplots!, setvalue!

include("simulationhelperfunctions.jl")
include("workflow.jl")
include("plotformatting.jl")

end  # module RenewalDiffInDiff
