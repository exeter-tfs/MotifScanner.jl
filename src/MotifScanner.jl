module MotifScanner

using BioSequences
using Plots
using StatsBase
using DataFrames

export loadmeme, scanmotif, scanmotstats, consensus, plotseq, plotletter!, seqlogo!, seqlogo

include("loadmeme.jl")
include("scanning.jl")
include("seqlogodata.jl")
include("plotseqlogo.jl")
end
