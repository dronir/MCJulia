# emcee Julia

# Ensemble sampling MCMC inspired by the emcee package for Python
# http://danfm.ca/emcee/

# ref.:
# Goodman & Weare, Ensemble Samplers With Affine Invariance
#   Comm. App. Math. Comp. Sci., Vol. 5 (2010), No. 1, 65–80

abstract Blob

# Random generator for the Z distribution of Goodman & Weare, where
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

# Return lnprobs and blobs at given position
function get_lnprob(S::Sampler, p::Array{Float64,2})
	return 0,0
end


function sample(S::Sampler, p0::Array{Float64,2}, N::Int64, thin::Int64, storechain::Bool)
	halfk = fld(S.n_walkers, 2)

	p = p0
	lnprob, blobs = get_lnprob(S, p)

	# Add N/thin columns of zeroes to the Sampler's chain and ln_posterior
	if storechain
		S.chain = hcat(S.chain, zeros(Float64, (S.n_walkers, S.dim, fld(N,thin))))
		S.ln_posterior = hcat(S.ln_posterior, zeros(Float64, (S.n_walkers, fld(N,thin))))
	end

	randZ() = randZ(S.a)

	first = 1 : halfk
	second = halfk+1 : S.n_walkers

	# FIXME
	#i0 = S.iteration
	i0 = 0
	for i = i0+1 : i0+N
		for t in [(first, second), (second, first)]
			active, passive = t
			l_act = length(active)
			l_pas = length(passive)
			for k in active
				X0 = p[k,:]
				choice = passive[randi(l_pas)]
				X1 = p[choice,:]
				z = randZ()
				Y = X1 + z * (X0 - X1)
				new_lnprob = get_lnprob(S, Y)
				q = (S.dim - 1) * log(z) + new_lnprob - lnprob
				if q > log(rand())
					lnprob = new_lnprob
					p[k,:] = Y
				end
				if (i-i0) % thin == 0
					S.ln_posterior[k, i/thin] = lnprob
					S.chain[k, :, i/thin] = p[k,:]
				end
			end
		end
	end
end

# Reset a Sampler to state after construction.
function reset(Sampler)
end