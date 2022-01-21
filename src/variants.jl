


### variant table
# chrom start stop id type = {var, ins, del} index refseq altseq

loadrefseqs(vartable, fastafile, mml=40) = loadrefseqs(vartable.chrom, vartable.start, vartable.stop, vartable.ref, vartable.alt, fastafile, mml)
#### currently only works for variants not indels
function loadrefseqs(chroms, starts, stops, refs, alts, file, mml=40)
    reader = open(FASTA.Reader, file, index=string(file, ".fai"))
    
    cr = 0
    ca = 0
    seqs = DataFrame(seqstart=Int[], seqstop=Int[], refind=UnitRange{Int}[], refseq=LongSequence{DNAAlphabet{4}}[], altseq=LongSequence{DNAAlphabet{4}}[])
    for (c, s, e, r, a) in zip(chroms, starts, stops, refs, alts)
        start = s - mml
        stop  = e + mml
        refseq = FASTA.extract(reader, DNAAlphabet{4}(), c, start:stop) 
        ind = (s:e) .- start .+ 1
        leftind = 1:(first(ind) - 1)
        rightind = (last(ind) + 1):length(refseq)
        altseq = refseq[leftind]*LongDNASeq(a)*refseq[rightind]
        @assert length(refseq) == length(altseq)
        if string(refseq[ind]) == r
            cr += 1
        elseif string(refseq[ind]) == a
            ca += 1
        else
            @show c, s, e, r, a, ind, refseq
            break
        end
        push!(seqs, (start, stop, ind, refseq, altseq))
    end
    
    close(reader)
    seqs
end

function motifscanall(seqtable, motifs)
    
    dfs = DataFrame[]
    for row in eachrow(seqtable)
        df = scanmots(row.refseq, row.altseq, row.refind, motifs)
        df[!, :ID] .= row.ID
        push!(dfs, df)
    end
    reduce(vcat, dfs)
end
## for variants only

function motifscanseq(seq, ind, motif)
    n = size(motif.pwm, 2)
    start = max(first(ind) - n + 1, 1)
    stop  = min(last(ind)  + n - 1, length(seq))
    seq[start:stop], start
end


function scanmots(refseq, altseq, ind, motifs)
    
    
    df = DataFrame(MotifName=String[], MotifID=String[], MotSeqStart=Int[], RefMaxScore=Float64[], RefStart=Int[], RefStop=Int[], RefStrand=String[],
                                                                            AltMaxScore=Float64[], AltStart=Int[], AltStop=Int[], AltStrand=String[],
                                                                            RefPrMax=Float64[], AltPrMax=Float64[], LR_RefAlt=Float64[], PR_RefAlt=Float64[])
    @showprogress for m in motifs
        maxscore = sum(maximum(m.pbg, dims=1))
        refmseq, refstart = motifscanseq(refseq, ind, m)
        altmseq, altstart = motifscanseq(altseq, ind, m)
        
        refres = scanmax(refmseq, m)
        altres = scanmax(altmseq, m)

        ref_prmax = first(refres)/maxscore
        alt_prmax = first(altres)/maxscore
        lr_refalt = first(refres) - first(altres)
        pr_refalt = ref_prmax - alt_prmax

        
        push!(df, (m.name, m.id, refstart, refres..., altres..., ref_prmax, alt_prmax, lr_refalt, pr_refalt))
    end
    df
    
end