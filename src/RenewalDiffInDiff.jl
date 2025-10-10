# functions for this analysis 

module RenewalDiffInDiff

import CSV
import ReverseDiff

using CairoMakie: Axis, Colorbar, Label, Legend, TopLeft
using CairoMakie: hidespines!, hidexdecorations!, hideydecorations!, scatter!
using DataFrames: DataFrame, insertcols!, rename!, select!
using Dates: Date
using Distributions: Binomial, Distribution, Uniform
using DrWatson: @ntuple, datadir, load, safesave
using Random: AbstractRNG, Xoshiro, default_rng
using RenewalDiD: RenewalDiD, InterventionArray, InterventionMatrix
using RenewalDiD: map_DataFrame, simulationu0
using StatsBase: ordinalrank, sample
using Turing: AutoReverseDiff, MCMCThreads, NUTS, Prior
using Turing: maximum_a_posteriori, maximum_likelihood 

export simulationdir
# consts.jl
export UKNATIONS, UKPOPULATION2020
# simulationhelperfunctions.jl
export largepop, simu0, smallpop
# loaddata.jl
export loadukmaskdata, ukinterventions, ukobservedcasesmatrix
# workflow.jl
export analysisworkflow
export indexesformap
export loadsamples 
export maximumaposterioriworkflow
export maximumlikelihoodworkflow
export mcmcworkflow
export priorsworkflow
# plotformatting.jl
export formataxis!, labelplots!, setvalue!
# setoffsets.jl
export addoffsetstointerventionarray

simulationdir(args...) = datadir("simulations", args...)

include("consts.jl")
include("simulationhelperfunctions.jl")
include("loaddata.jl")
include("workflow.jl")
include("plotformatting.jl")
include("setoffsets.jl")

end  # module RenewalDiffInDiff
