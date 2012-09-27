load("emcee.jl")
import MCJulia.*

f(X) = -(X[1]-1.0)^2

walkers = 100
S = Sampler(walkers, 1, f)
p0 = rand((walkers,1)) * 10 - 5

tic()
sample(S, p0, 100, 1, true)
toc()

flat = flat_chain(S)
save_chain(S, "chain.txt")