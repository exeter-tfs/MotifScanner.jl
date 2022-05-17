


### variant table
# chrom start stop id type = {var, ins, del} index refseq altseq

loadrefseqs(vartable, fastafile; mml=40, reffield=:ref, altfield=:alt) = loadrefseqs(vartable.chrom, vartable.start, vartable.stop, vartable[!, reffield], vartable[!, altfield], fastafile, mml)
#### currently only works for variants where Ref is single base
function loadrefseqs(chroms, starts, stops, refs, alts, file, mml=40)
    reader = open(FASTA.Reader, file, index=string(file, ".fai"))
    
    cr = 0
    ca = 0
    seqs = DataFrame(seqstart=Int[], seqstop=Int[], Kind=String[], refind=UnitRange{Int}[], altind=UnitRange{Int}[], refseq=LongSequence{DNAAlphabet{4}}[], altseq=LongSequence{DNAAlphabet{4}}[])
    for (c, s, e, r, a) in zip(chroms, starts, stops, refs, alts)
        start = s - mml
        stop  = e + mml
        refseq = FASTA.extract(reader, DNAAlphabet{4}(), c, start:stop) 
        refind = (s:e) .- start .+ 1
        altind = first(refind) .+ (1:length(a)) .- 1
        leftind = 1:(first(refind) - 1)
        rightind = (last(refind) + 1):length(refseq)
        seqa = LongDNASeq(a)
        altseq = refseq[leftind]*seqa*refseq[rightind]
        #@assert length(refseq) == length(altseq)


        @assert seqa == altseq[altind]
        if string(refseq[refind]) == r
            cr += 1
        elseif string(refseq[refind]) == a
            ca += 1
        else
            @show c, s, e, r, a, refind, altind refseq
            break
        end

        if length(r) == length(a) == 1
            kind = "Var"
        else
            kind = "Indel"
        end
        push!(seqs, (start, stop, kind, refind, altind, refseq, altseq))
    end
    
    close(reader)
    seqs
end

function motifscanall(seqtable, motifs)
    
    dfs = DataFrame[]
    for row in eachrow(seqtable)
        df = scanmots(row.refseq, row.altseq, row.refind, row.refind, motifs)
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



function scanmots(refseq, altseq, refind, altind, motifs)
    
    
    df = DataFrame(MotifName=String[], MotifID=String[], RefMotSeqStart=Int[], AltMotSeqStart=Int[], RefMaxScore=Float64[], RefStart=Int[], RefStop=Int[], RefStrand=String[],
                                                                            AltMaxScore=Float64[], AltStart=Int[], AltStop=Int[], AltStrand=String[],
                                                                            RefPrMax=Float64[], AltPrMax=Float64[], LR_RefAlt=Float64[], PR_RefAlt=Float64[])
    @showprogress for m in motifs
        maxscore = sum(maximum(m.pbg, dims=1))
        refmseq, refstart = motifscanseq(refseq, refind, m)
        altmseq, altstart = motifscanseq(altseq, altind, m)
        
        refres = scanmax(refmseq, m)
        altres = scanmax(altmseq, m)

        ref_prmax = first(refres)/maxscore
        alt_prmax = first(altres)/maxscore
        lr_refalt = first(refres) - first(altres)
        pr_refalt = ref_prmax - alt_prmax

        
        push!(df, (m.name, m.id, refstart, altstart, refres..., altres..., ref_prmax, alt_prmax, lr_refalt, pr_refalt))
    end
    df
    
end