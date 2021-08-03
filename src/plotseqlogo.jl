
### definitely need to update this away from something so quick and dirty, ideally use plot recipes

### A letter is a Vector{Vector{Tuple{Float64, Float64}}}.

### Plots.Shape requires a Vector{Tuple} functions to convert this to and from a matrix
lettermatrix(L) =  [mapreduce(t -> [t[1] t[2]], vcat, t) for t in L]
lettertuple(M)  =  [map(r -> (r[1], r[2]), eachrow(T)) for T in M]


### Scale letter to unit square
function normletter(l)
    TM = lettermatrix(l)
    tmax = maximum(mapreduce(t -> maximum(t, dims=1), vcat, TM), dims=1)
    TN = [T./tmax for T in TM]
    lettertuple(TN) 
end

### offset a letter in x and y
function offsetletter(L, x=0, y=0)
    M = lettermatrix(L)
    MO = [s .+ [x y] for s in M]
    lettertuple(MO)
end

### scale a letter in x and y
function scaleletter(L, x=1, y=1)
    M = lettermatrix(L)
    MO = [s .* [x y] for s in M]
    lettertuple(MO)
end


### plot a letters scaling in x and y
function plotletter!(off, scale, letter; kwargs...)
    plotletter!(offsetletter(scaleletter(letter, scale[1], scale[2]), off[1], off[2]) ; kwargs...)
end

### plot the letters
function plotletter!(letter ; bgc=:white, kwargs...)
    p = plot!(Shape(letter[1]), lab="", line=stroke(0); kwargs...)
    for l in letter[2:end]
        plot!(Shape(l), lab="", c=bgc, line=stroke(0))
    end
    p
end

### plot seq logo
function seqlogo(mot; rc=false, xo=0.0, yo=0.0, kwargs...)
    plot()
    seqlogo!(mot, rc=rc, xo=xo, yo=yo; kwargs...)
end


### plot seq logo at specified offset
function seqlogo!(mot; weight="normal", s=10 ; rc=false, xo=0.0, yo=0.0, label=:none, pseudo=1e-6, kwargs...)
    acgt = motif_letter_data(weight=weight, s=s)
    pwm = ifelse(rc, rcm(mot.pwm), mot.pwm) .+ pseudo
    H = log2(4) .+ sum(pwm.*log2.(pwm), dims=1)
    bases = (DNA_A, DNA_C, DNA_G, DNA_T)
    S = pwm.*H
    p = plot!(; kwargs...)
    off = xo .+ 0.5
    cc = palette(:default)[1:4]
    cc = [:green, :steelblue, :orange, :red]
    for c in eachcol(S)
        si = sortperm(c)
        yoff = yo
        for s in si
            v = c[s]
            b = bases[s]
            cl = cc[s]
           plotletter!([off, yoff], [1.0, v], acgt[b], c=cl) 
           yoff += v
        end
        
        off += 1
    end
    
    if label == :left
        annotate!((xo, yo .+ 1, text(mot.name, font(8, :right, :middle))))
    elseif label == :right
        annotate!((off, yo .+ 1, text(mot.name, font(8, :left, :middle))))
        
    elseif label == :top
        annotate!(((off .+ xo)/2, yo + 2, text(mot.name, font(8, :top, :center))))
    end
    p
end


### for plotting sequences
function plotseq(seq; weight="normal", s=10, xo = 0.5, yo = 0, kwargs...)
    acgt=motif_letter_data(weight=weight, s=s)
    p = plot(; kwargs...)
    
    offx = xo
    offy = yo
    for (i, s) in enumerate(seq)
        plotletter!([offx, offy], [1.0, 0.25], acgt[s], c=:black)
        offx += 1
    end
    p
end