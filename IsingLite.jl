module IsingLite

using PyPlot

export neighbors,
       nspins,
       spingrid,
       namefunc,
       wolff!,
       diagram

include("utils.jl")
include("wolff.jl")

end # module