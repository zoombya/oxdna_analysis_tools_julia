#@time main(["test\\init.top","test\\init.dat","test\\trajectory_md.dat","test\\out.dat"])
include("./align.jl")
include("./reader.jl")

function main(ARGS)
  top_file , ref_conf, alignment_conf, out_conf = ARGS 
  # figure out how many bases and strands we work with
  top_info = open(top_file) do f
      read_top(f)
  end
  # read refference conf 
  reference_conf = open(ref_conf) do f
    read_conf(f, top_info)
  end

  # retrieve number of threads which will be the number of confs we align at the same time 
  n = Threads.nthreads()
  # place to store the output configurations before write to file
  out_confs = fill("", n) 

  if(isfile(out_conf))
    rm(out_conf)
  end



  open(out_conf,"a") do out
    # go through the trajectory file
    open(alignment_conf) do trajectory
        conf_count = read_conf_count(trajectory, top_info)
        for i = 1:n:conf_count
          println("i: ",i,"/",conf_count)
          #fix end case
          if(i + n > conf_count)
            n = conf_count - i 
          end
          #preread the confs 
          confs = [read_conf(trajectory, top_info) for j in 1:n]
          
          #perform alignment 
          Threads.@threads for j=1:n
            confs[j].positions = align(confs[j].positions,reference_conf.positions)
            confs[j].a1s = align(confs[j].a1s,reference_conf.a1s)
            confs[j].a3s = align(confs[j].a3s,reference_conf.a3s)
            out_confs[j] = conf_to_str(confs[j], top_info)
          end
          
          #write output 
          write(out,join(out_confs))
        end
    end
  end
end 