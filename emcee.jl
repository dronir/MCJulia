# MC Julia

# Ensemble sampling MCMC inspired by the emcee package for Python
# http://danfm.ca/emcee/

# ref.:
# Goodman & Weare, Ensemble Samplers With Affine Invariance
#   Comm. App. Math. Comp. Sci., Vol. 5 (2010), No. 1, 65–80

require("options.jl")

module MCJulia

export Blob, Sampler, sample, reset

import Base.*
import OptionsMod.*


abstract Blob

# Random generator for the Z distribution of Goodman & Weare, where
# p(x) = 1/sqrt(x) when 1/a <= x <= a.
randZ(a::Float64) = ((a - 1.0) * rand() + 1.0)^2 / a

# The Sampler type is the interface between the user and the machinery.
type Sampler
	n_walkers::Int64
	dim::Int64
	probfn::Function
	a::Float64
	chain::Array{Float64, 3}
	ln_posterior::Array{Float64,2}
	blobs::Array{Blob, 2}
	acceptance_fraction::Float64
	threads::Int64
	args::(Any...)
end

function Sampler(k, dim, f, opts::Options)
	@defaults opts a=2.0 threads=1 accpt=0.0
	@defaults opts chain = zeros(Float64, (k, dim, 0))
	@defaults opts ln_p = zeros(Float64, (k, 0))
	@defaults opts blobs = Array(Blob, (k, 0))
	@defaults opts args = ()
	S = Sampler(k,dim,f,a,chain,ln_p,blobs,accpt,threads,args)
	@check_used opts
	return S
end

# Return lnprobs and blobs at given position
function get_lnprob(S::Sampler, pos::Array{Float64,2})
	lnprob = zeros(Float64, S.n_walkers)
#	blobs = Array(Blob, S.n_walkers)
#	b = false
	for k = 1:S.n_walkers
		p = S.probfn(vec(pos[k,:]))
		if length(p) == 1
			lnprob[k] = p
		elseif length(p) == 2
			b = true
			lnprob[k] = p[1]
			blobs[k] = p[2]
		else
			error("something weird came out of S.lnprobf()")
		end
	end
	# FIXME
	return lnprob, None
end


function sample(S::Sampler, p0::Array{Float64,2}, N::Int64, thin::Int64, storechain::Bool)
	k = S.n_walkers
	halfk = fld(k, 2)

	p = copy(p0)
	lnprob, blobs = get_lnprob(S, p)

	i0 = size(S.chain, 3)

	# Add N/thin columns of zeroes to the Sampler's chain and ln_posterior
	if storechain
		S.chain = cat(3, S.chain, zeros(Float64, (k, S.dim, fld(N,thin))))
		S.ln_posterior = cat(2, S.ln_posterior, zeros(Float64, (k, fld(N,thin))))
	end

	first = 1 : halfk
	second = halfk+1 : k
	subsets = [(first, second), (second, first)]

	# FIXME
	#i0 = S.iteration
#	i0 = 0
	for i = i0+1 : i0+N
		for t in subsets
			active, passive = t
#			l_act = length(active)
			l_pas = length(passive)
			for k in active
				X0 = p[k,:]
				choice = passive[randi(l_pas)]
				X1 = p[choice,:]
				z = randZ(S.a)
				Y = X1 + z * (X0 - X1)
				new_lnprob = S.probfn(Y)
				log_ratio = (S.dim - 1) * log(z) + new_lnprob - lnprob[k] # FIXME: lnprob
				if log_ratio < log(rand())
					lnprob[k] = new_lnprob
					p[k,:] = Y
				end
				if (i-i0) % thin == 0
					S.ln_posterior[k, i/thin] = lnprob[k]
					S.chain[k, :, fld(i,thin)] = vec(p[k,:])
				end
			end
		end
	end
end

# Reset a Sampler to state after construction.
function reset(S::Sampler)
	k = S.n_walkers
	S.chain = zeros(Float64, (k, S.dim, 0))
	S.ln_posterior = zeros(Float64, (k, 0))
	S.blobs = Array(Blob, (k, 0))
	S.acceptance_fraction = 0.0
	return S
end

end