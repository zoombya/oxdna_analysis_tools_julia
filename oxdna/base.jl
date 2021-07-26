using Base: Float64, Int128, File, readuntil_vector!

struct TopInfo
    base_count::Int128
    strand_count::Int128
    TopInfo(base_count,strand_count) = new(base_count,strand_count)
end

mutable struct Conf
    t::Int128
    b::Vector{Float64}
    E::Vector{Float64}
    positions::Matrix{Float64}
    a1s::Matrix{Float64}
    a3s::Matrix{Float64}
    Conf(t,b,E,positions,a1s,a3s) = new(t,b,E,positions,a1s,a3s)
end

function read_top(f::IO)
    """
      Figure out the number of bases and strands in a provided topology file.
    """
    base_count, strand_count =  map(x->parse(Int128, x),split(readline(f)))
    return TopInfo(base_count, strand_count)
end

function read_line_count_of(f::IO)
    """
      Figure out the number of lines in a given file.
    """
    linecounter = 0
    for l in eachline(f)
        linecounter += 1
    end
    return linecounter
end

function read_conf_count(f::IO, ti::TopInfo)
    """
        Figure out the number of configurations in a given trajectory. 
    """
    line_count = read_line_count_of(f)
    seek(f,0)
    return line_count ÷ (ti.base_count + 3)
end


function binary_index(f::IO, buff_size = 52428800)
    indices = open("test/trajectory_md.dat") do f
        val = 0x74; #t
        buff_size = 52428800
        offset = 0
        cur_offset = 0
        indices = []
        count = 0
        firstRead = true
        file_size = stat(f).size
        while !eof(f)  
            buff  = Array{UInt8,1}(undef, buff_size)
            readbytes!(f, buff)
            buff_offsets = findall(isequal(val), buff)
            # to keep up with normal world we have to use 0 based indexing 
            buff_offsets .-= 1 # add the offset !!! 
            n = length(buff_offsets)
            if(n>=0)
                a = [[buff_offsets[i], buff_offsets[i+1] - buff_offsets[i], i-1] for i in 1:(n-1)]
            
                indices = vcat(
                    indices, a
                )
                #offset = offset buff_offsets[end]
                print(offset)
            end
            count += 1
            if(count == 2)
                break
            end
        end
        indices
    end
    print(indices)
end

function read_conf(f::IO, ti::TopInfo)
    """
        Read in a configuration from a dat file.
    """
    # parse header
    t = parse(Int128,split(readline(f),"=")[2])
    b = map(x->parse(Float64, x),split(split(readline(f),"=")[2]))
    E = map(x->parse(Float64, x),split(split(readline(f),"=")[2]))
    # create a conf to store parsed values
    conf = Conf(t, b, E,
                 fill(0.0, ti.base_count,3),
                 fill(0.0, ti.base_count,3), 
                 fill(0.0, ti.base_count,3))
    # parse conf
    for i = 1:ti.base_count
        # parse the line and convert values to Float64
        conf_line = map(x->parse(Float64, x), split(readline(f)))
        # extract the positions
        conf.positions[i,1] = conf_line[1]
        conf.positions[i,2] = conf_line[2]
        conf.positions[i,3] = conf_line[3]
        # extract the a1s 
        conf.a1s[i,1] = conf_line[4]
        conf.a1s[i,2] = conf_line[5]
        conf.a1s[i,3] = conf_line[6]
        # extract the a3s
        conf.a3s[i,1] = conf_line[7]
        conf.a3s[i,2] = conf_line[8]
        conf.a3s[i,3] = conf_line[9]
    end
    return conf
end

function write_conf(f::IO,conf::Conf,ti::TopInfo)
    """
      Output a configuration to a dat file.
    """ 
    write(f,conf_to_str(conf,ti))    
end

function  conf_to_str(conf::Conf,ti::TopInfo)
    out = []
    push!(out, string("t = ", conf.t))
    push!(out, string("b = ", conf.b[1], " ", conf.b[2], " ", conf.b[3]))
    push!(out, string("E = ", conf.E[1], " ", conf.E[2], " ", conf.E[3]))
    for i = 1:ti.base_count
        push!(out, string(conf.positions[i,1], " ", conf.positions[i,2], " ", conf.positions[i,3], " ",
                          conf.a1s[i,1], " ", conf.a1s[i,2], " ", conf.a1s[i,3], " ",
                          conf.a3s[i,1], " ", conf.a3s[i,2], " ", conf.a3s[i,3], " ",
                          0, " ", 0, " ", 0," ",
                          0, " ", 0, " ", 0))
    end
    return join(out,"\n")   
end

function inbox!(conf::Conf,ti::TopInfo)
    """
        Bring a configuration into the simulation box.
    """
        real_mod = (n,m) -> begin
        for i = 1:first(size(n))
            n[i,1] = ((n[i,1] % m[1]) + m[1]) % m[1]
            n[i,2] = ((n[i,2] % m[2]) + m[2]) % m[2]
            n[i,3] = ((n[i,3] % m[3]) + m[3]) % m[3]
        end
        return n
    end 


    calc_PBC_COM = (conf, ti) -> begin
        π2 = 2.0 * π    
        angle = conf.positions * π2
        for i = 1:ti.base_count
            angle[i,1] /= conf.b[1]
            angle[i,2] /= conf.b[2]
            angle[i,3] /= conf.b[3]
        end
        n = ti.base_count
        cm = [[sum(map(cos,angle[:,1])) / n, sum(map(sin,angle[:,1])) / n],
              [sum(map(cos,angle[:,2])) / n, sum(map(sin,angle[:,2])) / n],
              [sum(map(cos,angle[:,3])) / n, sum(map(sin,angle[:,3])) / n]]

        pbc_com = copy(conf.b)
        pbc_com /= π2
        for i = 1:3
            pbc_com[i] *=  (atan(-cm[i][2] , -cm[i][1]) + π)
        end 

        return pbc_com
    end

    target = [conf.b[1] / 2, conf.b[2] / 2, conf.b[3] / 2]
    center = calc_PBC_COM(conf,ti)

    factor = target - center
    for i=1:ti.base_count
        conf.positions[i,1] += factor[1]
        conf.positions[i,2] += factor[2]
        conf.positions[i,3] += factor[3]
    end
    old_poses = copy(conf.positions)
    new_poses = real_mod(copy(conf.positions), 
                              conf.b)

    conf.positions += new_poses - old_poses
end