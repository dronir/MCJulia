# emcee Julia

# Ensemble sampling MCMC based on the emcee package for Python
# http://danfm.ca/emcee/

# ref.:
# Goodman & Weare, Ensemble Samplers With Affine Invariance
#   Comm. App. Math. Comp. Sci., Vol. 5 (2010), No. 1, 65–80

abstract Blob

# Random generator for the Z distribution of Goodman&Weare, where
# p(x) = 1/sqrt(x) when 1/a <= x <= a.
randZ(a::Float64) = ((a - 1.0) * rand() + 1.0)^2 / a

# The Sampler type is the interface between the user and the machinery.
type Sampler
	n_walkers::Integer
	dim::Integer
	probfn::Function
	a::Float64
	chain::Array{Float64, 3}
	ln_posterior::Array{Float64,2}
	blobs::Array{Blob, 2}
	acceptance_fraction::Float64
	threads::Int64
	args::(Any...)
end

# Minimal constructor for Sampler
function Sampler(n_walkers::Integer, dim::Integer, probfn::Function)
	n_walkers % 2 == 0 ? true : error("n_walkers must be even!")
	n_walkers > 2*dim ? true : error("n_walkers must be greater than 2*dim!")

	chain = zeros(Float64, (n_walkers, dim, 0))
	ln_post = zeros(Float64, (n_walkers, 0))
	blobs = Array(Blob, (n_walkers, 0))

	a = 2.0
	accpt = 0.0
	threads = 1
	Sampler(n_walkers,dim,probfn,a,chain,ln_post,blobs,accpt,threads,())
end

type Ensemble
	# Type for the ensemble
end

function reset(Sampler)
	# This will reset a Sampler.
end