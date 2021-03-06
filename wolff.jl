
function recursion_wolff!(grid::Array{Int, 2},         # Spin grig
                          cluster::BitArray{2},        # Cluster grid
                          i::Int,                      # Row coordinate
                          j::Int;                      # Column coordinate
                          h::Float64            = 0.0, # External field
                          temp::Float64         = 1.0) # Temperature

    cluster[i, j] = true
    for (x, y) in neighbors(grid, i, j)
        if cluster[x,y] == false && grid[x,y] == grid[i, j]
# Add the spin to cluster
            if rand() < 1 - exp(-2/temp) 
                recursion_wolff!(grid, cluster, x, y, h=h, temp=temp) 
            end
        end
    end
end

function wolff!(grid::Array{Int, 2};
                h::Float64           = 0.0,  # External field
                temp::Float64        = 1.0,  # Temperature
                iters::Integer       = 30,   # Number of iterations
                plot::Bool           = false, # Plot flag
                verbose::Bool        = false) # Verbose flag
        
        m = Float64[]
        for i in 1:iters
# Randomly pick a position within the grid
                x   = rand(1:size(grid, 1))
                y   = rand(1:size(grid, 2))

# Create a cluster grid full of "falses"
                cluster = falses(size(grid))

# Run the algorithm and calculate the final magnetization
                recursion_wolff!(grid, cluster, x, y, h=h, temp=temp)
# Calculate ΔE
           #     ΔE = -2h*clusterspin(grid, cluster)
# Change spin accordingly
           #     if ΔE <= 0 || rand() < exp(-ΔE) flip!(grid, cluster) end
                flip!(grid, cluster)
                push!(m, susceptibility(grid, temp))
        end

        if verbose println("(T=$temp) Wolff ended with magnetization $(m[end]) after $iters iterations") end
        if plot
                PyPlot.scatter(1:length(m), m, color="blue", s = 10)
                PyPlot.plot(1:length(m), m, "-", color="blue")
                PyPlot.title("Wolff for T=$temp for size 100 ($(length(m)) iterations and H=$h)")
                PyPlot.xlabel("Number of Iterations")
                PyPlot.ylabel("Magnetization")
                PyPlot.ylim(0,1.1)
                PyPlot.savefig("wolff_T=$temp.png")
                PyPlot.close()
        end

        return m
end