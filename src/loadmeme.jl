
background_human_zero_markov()   = [0.2953, 0.2047, 0.2047, 0.2953]
background_xenopus_zero_markov() = [0.2996, 0.2004, 0.2004, 0.2996]

### Load meme into named tuple
function loadmeme(file, w = 1e-2, background=background_human_zero_markov())
   
    motifname = ""
    motifid   = ""
    mode = 1
    data = Vector{Vector{Float64}}()
    for line in eachline(file)
        if mode == 1
            if startswith(line, "MOTIF")
                fields = split(line)
                motifid   = fields[2]
                motifname = fields[3]
            elseif startswith(line, "letter-probability")
                mode = 2
            end
        else
            if startswith(line, "URL")
                mode = 1
                break
            end
            v = parse.(Float64, split(line))
            push!(data, v)
        end
    end
    pwm = reduce(hcat, data)
    pwm .+= w
    pwm ./= sum(pwm, dims=1)
    
    pbg = log2.(pwm) .- log2.(background)
    (name=motifname, id=motifid, pwm=pwm, pbg=pbg)
end

function loadhomer(file, w = 1e-2, background=background_human_zero_markov())

    io = open(file)
    line = readline(io)

    fields = split(line)
    consensus = fields[1][2:end]
    longname = fields[2]
    name = first(split(longname, "/"))
    sct = parse(Float64, fields[3])


    data = Vector{Vector{Float64}}()
    for line in eachline(io)
        v = parse.(Float64, split(line))
        push!(data, v)
    end



    pwm = reduce(hcat, data)
    pwm .+= w
    pwm ./= sum(pwm, dims=1)
    
    pbg = log2.(pwm) .- log2.(background)
    (name=name, id=longname, pwm=pwm, pbg=pbg, sct=sct, consensus=consensus)
end
