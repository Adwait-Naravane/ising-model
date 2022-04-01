function stepmetropolis!(grid::Array{Int, 2}; # Spin grid
    h::Float64    = 0.0, # External field
    temp::Float64 = 1.0) # Temperature

# Randomly pick a position within the grid
x   = rand(1:size(grid, 1))
y   = rand(1:size(grid, 2))

# Calculate the ΔE for switching the spin
m      = nspins(grid, x, y) |> sum
eplus  = -m - h
ΔE     = -2eplus * grid[x,y]

# Change spin accordingly
if ΔE <= 0 || rand() < exp(-ΔE/temp)
grid[x,y] *= -1
end
end

# Several steps of the heat bath algorithm on a Ising's spin graph
function metropolis!(grid::Array{Int, 2};  # Spin grid
                     h::Float64=0.0,       # External field
                     temp::Float64=1.0,    # Temperature
                     iters::Integer=50000, # Number of iterations
                     plot::Bool=false,      # Plot flag
                     verbose::Bool=false)   # Verbose flag

m = Float64[]
for i in 1:iters
stepmetropolis!(grid, h=h, temp=temp)
push!(m, susceptibility(grid, temp))
# Must find a better way to decide convergence
end

if verbose println("(T=$temp) Metropolis ended with susceptibility $(m[end]) after $(length(m)) iterations") end
if plot
PyPlot.plot(1:length(m), m, "o", color="blue")
PyPlot.plot(1:length(m), m, "-", color="blue")
PyPlot.title("Metropolis for T=$temp for size $size ($(length(m)) iterations and H=$h)")
PyPlot.xlabel("Number of Iterations")
PyPlot.ylabel("Magnetization")
PyPlot.ylim(0,1.1)
PyPlot.savefig("metropolis_T=$temp.png")
PyPlot.close()
end

return m
end