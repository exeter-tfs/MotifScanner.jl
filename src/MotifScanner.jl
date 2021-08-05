module MotifScanner

using BioSequences
using Plots
using StatsBase
using DataFrames
import PyPlot

## This is a package for scanning motifs in DNA sequences and plotting sequence motifs.

export loadmeme, loadmemelibrary, loadhomer, scanmotif, scanmotstats, consensus, plotseq, plotletter!, seqlogo!, seqlogo

include("loadmeme.jl")
include("scanning.jl")
include("plotseqlogo.jl")
include("seqlogodata.jl")
end
