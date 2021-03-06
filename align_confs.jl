#@time main(["test\\init.top","test\\init.dat","test\\trajectory_md.dat","test\\out.dat"])
include("oxdna/base.jl")
include("oxdna/align.jl")

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
  inbox!(reference_conf,top_info)

  if(isfile(out_conf))
    rm(out_conf)
  end

  open(out_conf,"a") do out
    # go through the trajectory file
    open(alignment_conf) do trajectory
        conf_count = read_conf_count(trajectory, top_info)
        print(time)
        for i = 1:conf_count
          println("i: ",i,"/",conf_count)
          conf = read_conf(trajectory, top_info)
          inbox!(conf,top_info)
          # do alignment 
          align!(conf,reference_conf)
          write_conf(out,conf , top_info)
        end
    end
  end
end 
