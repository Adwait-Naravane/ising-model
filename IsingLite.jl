module IsingLite

using PyPlot

export neighbors,
       nspins,
       spingrid,
       namefunc,
       wolff!,
       diagram, 
       metropolis!

include("utils.jl")
include("wolff.jl")
include("metropolis.jl")

end # module