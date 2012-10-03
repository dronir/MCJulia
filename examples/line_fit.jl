# Simple data fitting example: Fitting a line to noisy points. 
# The example will readily generalize to any simple data fitting problem.

# First generate the data. We will make a data set consisting of (x,y)
# pairs of the form y = a*x + b + noise.
# These are the parameters used to generate the data:
a = 2.0
b = -4.0
sigma = 0.5 # noise sigma

# Generate ten data points evenly spaced in x:
n_data = 10
X = linspace(0.0, 10.0, n_data)
Y = a*X + b + randn(n_data)*sigma

# In a real application, everything above would come from the data.
# Now start the actual estimation process.

load("mcjulia.jl")
import MCJulia.*

# When fitting data, our likelihood function has the form exp(-ChiSq), where
# ChiSq = sum(model - data)^2 / sigma. We use a wide normal prior distribution
# for a and b, and an exponential prior p(x) ~= exp(-x) for sigma.
# Our log-probability function has two extra arguments after the parameter vector,
# giving the x and y values of the data points.
function log_probability(parameters::Array{Float64}, X_data::Array{Float64}, Y_data::Array{Float64})
	a = parameters[1]
	b = parameters[2]
	sigma = parameters[3]
	if sigma <= 0
		return -Inf
	end
	Y_model = a*X_data + b
	ChiSq = sum(((Y_model - Y_data)/sigma).^2)
	return -ChiSq - (a/100.0)^2 - (b/100.0)^2 - sigma
end

# We give our data points to the sampler in the extra arguments tuple.
args = (X, Y)

# Set up the sampler. It is good to use a large number of walkers.
dim = 3
walkers = 100
S = Sampler(walkers, dim, log_probability, args)

# Generate starting positions for the walkers from a N(0,10)
# distribution, squaring it for sigma. The ideal starting positions
# will be in a spherical shell around the high-probability region, but
# as long as they are close enough, the walkers will relax towards
# the region fairly quickly.
p0 = randn((walkers, 3)) * 10
p0[:,3] = p0[:,3].^2

# Do a burn-in of 100 steps (10000 samples). Discard the samples
# but keep the final position of the walkers.
println("Burn-in...")
p = sample(S, p0, 100, 1, false)
println("acceptance ratio: $(S.accepted / S.iterations)")

# Start sample for another 100 steps, saving the chain every 10
# steps, starting at the end position of the burn-in.
println("Sampling...")
sample(S, p, 100, 10, true)
println("acceptance ratio: $(S.accepted / S.iterations)")

# Save and plot the chain.
println("Plotting...")
save_chain(S, "chain.txt")
run(`python plot_chains.py chain.txt`)
