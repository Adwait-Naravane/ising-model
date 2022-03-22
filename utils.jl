include("wolff.jl")
using  Statistics
using PyPlot
using KernelDensity, Distributions,Random,StatsPlots
randspin()                               = [1,-1][rand(1:2)] # Generate a random spin
spingrid(n::Int)                         = [randspin() for i in 1:n, j in 1:n] # Generate a random spin array
magnetization(a::Array{Int, 2})          = mean(a) |> abs # Get magnetizations of the grid
susceptibility(a::Array{Int, 2})         = var(a)
namefunc(f::Function)                    = "$f"[1:end-1]  # Get the simplified name of a function
nspins(a::Array{Int, 2}, i::Int, j::Int) = [a[x,y] for (x,y) in neighbors(a,i,j)] # Get surrounding spins
mode(a::Array{Float64}) = findmax(kde(a).density)


# Flip a spin
function flip!(grid::Array{Int, 2}, cluster::BitArray{2})
    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        if cluster[i, j] == true grid[i, j] *= -1 end
    end
end

# Return the sum of spins on a cluster
function clusterspin(grid::Array{Int, 2}, cluster::BitArray{2})
    spin = 0
    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        if cluster[i,j] == true spin += grid[i,j] end
    end
    return spin
end

# Return an array of neighbors coordinates (each one in a tuple of the kind (x,y))
function neighbors(grid::Array{Int, 2}, i::Integer, j::Int)
    n = Tuple{Int, Int}[]

    if i > 1             push!(n, (i-1, j  )) end
    if j > 1             push!(n, (i,   j-1)) end
    if i < size(grid, 1) push!(n, (i+1, j  )) end
    if j < size(grid, 2) push!(n, (i,   j+1)) end

    return n
end

# Phase diagram (magnetization by temperature) using given algorithm
function diagram(func::Function;
                 size::Integer      = 50,    # Size of the grid
                 ensembles::Integer = 300,    # Number of ensembles
                 h::Float64         = 0.0,   # External field
                 mintemp::Float64   = 2.0,   # Starting temperature
                 step::Float64      = 0.01,   # Step of temperatures
                 maxtemp::Float64   = 2.5,   # Final temperature
                 iters::Integer     = 150, # Number of the iterations
                 plot::Bool         = true,  # Plot flag
                 verbose::Bool      = true)  # Verbose flag

    name  = namefunc(func)
    temps = Float64[]
    mags  = Float64[]
    for t in mintemp:step:maxtemp
        m = mean([func(spingrid(size), h=h, temp=t, iters=iters, plot=false, verbose=false)[end] for i in 1:ensembles])
        if verbose println("(T=$t) Avg. magnetization after $name: $m") end

        push!(temps, t)
        push!(mags,  m)
    end

    if plot
        PyPlot.scatter(temps, mags, color="blue", s = 10)
        PyPlot.plot(temps, mags, "-", color="blue")
        PyPlot.plot([2.269, 2.269], [0, 0.5], "-" )
        PyPlot.grid()
        PyPlot.title("Phase Diagram ($name for size $size and h=$h)")
        PyPlot.xlabel("Temperature")
        PyPlot.ylabel("Mean magnetization")
        PyPlot.savefig("diag_magnetisation_$(name)_H=$(h)_size=$(size).png")
        PyPlot.close()
    end

    return temps, mags
end



diagram(wolff!)