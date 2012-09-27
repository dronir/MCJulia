load("emcee.jl")
import MCJulia.*
S = Sampler(8, 2, log_rosenbrock)
println(S)
p0 = rand((8,2)) * 2
sample(S, p0, 1000, 1, true)
flat = flat_chain(S)
save_chain(S, "chain.txt")
println(flat)