
## reverse complement
rcm(M) = reverse(reverse(M, dims=1), dims=2)


### loop for forward and reverse motif scanning
## no doubt there are efficiency savings available here
function scanmotif(seq, mot)
    n = length(seq)
    m = size(mot, 2)
    fscores = zeros(Float64, n)
    rscores = zeros(Float64, n)
    rot = rcm(mot)
    for i = 1:(n-m)
        for j = 1:m
            if seq[i + j - 1] == DNA_A
                fscores[i] += mot[1, j]
                rscores[i] += rot[1, j]
            elseif seq[i + j - 1] == DNA_C
                fscores[i] += mot[2, j]
                rscores[i] += rot[2, j]
            elseif seq[i + j - 1] == DNA_G
                fscores[i] += mot[3, j]
                rscores[i] += rot[3, j]
            elseif seq[i + j - 1] == DNA_T
                fscores[i] += mot[4, j]
                rscores[i] += rot[4, j]
            end
        end
        
    end
    
    
    fscores, rscores
end

### motif scann stats
function scanmotstats(mot, seq::T, thr=5) where{T}
        
    fs, rs = scanmotif(seq, mot.pbg)
    
    ## neg ecdf
    nec = ecdf([-fs ; -rs])
    
    fsi = findall(fs .> thr)
    rsi = findall(rs .> thr)
    
    res = DataFrame(Motif=String[], start=Int[], stop=Int[], score=Float64[], strand=String[], emp_p=Float64[], prmax=Float64[], seq=Vector{T}())
    
    maxscore = sum(maximum(mot.pbg, dims=1))
    
    n = size(mot.pbg, 2)
    for f in fsi
       push!(res, (mot.name, f, f + n - 1, fs[f], "+", nec(-fs[f]), fs[f]/maxscore, seq[f:f+n-1]))
    end
    
    for r in rsi
       push!(res, (mot.name, r, r + n - 1, rs[r], "-", nec(-rs[r]), rs[r]/maxscore, reverse_complement(seq[r:r+n-1])))
    end
    res
end

function consensus(mot)
    bases = (DNA_A, DNA_C, DNA_G, DNA_T)
    LongDNASeq(mapreduce(c -> bases[argmax(c)], vcat, eachcol(mot.pwm)))
end

function motifind(motname, mots)
    ms = findall(m -> occursin(motname, m.name), mots) 
    if isempty(ms)
       error("$motname not found") 
    elseif length(ms) > 1
        error("Multiple matches for $motname\n$([mots[m].name for m in ms])")
    else
       first(ms)
    end 
    
end

function consensus(motname, mots)
    ms = filter(m -> occursin(motname, m.name), mots) 
    if isempty(ms)
       error("$motname not found") 
    elseif length(ms) > 1
        error("Multiple matches for $motname\n$([m.name for m in ms])")
    else
        consensus(first(ms))
    end
end
