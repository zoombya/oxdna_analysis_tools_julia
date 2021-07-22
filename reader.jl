using Base: Float64, Int128

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
    write(f,join(out,"\n"))    
end