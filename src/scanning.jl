
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

##