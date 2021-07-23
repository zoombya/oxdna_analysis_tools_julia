#julia -t 10 .\thread_main.jl init.top init.dat .\trajectory_md.dat out.dat
include("./align.jl")
include("./reader.jl")

function main(args)
  top_file , ref_conf, alignment_conf, out_conf = args 
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
          end
          #write output 
          for j=1:n
            write_conf(out, confs[j] , top_info)
          end
        end
    end
  end
end 

@time main(ARGS)