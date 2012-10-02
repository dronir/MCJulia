# MC Julia

# Ensemble sampling MCMC inspired by the emcee package for Python
# http://danfm.ca/emcee/

# ref.:
# Goodman & Weare, Ensemble Samplers With Affine Invariance
#   Comm. App. Math. Comp. Sci., Vol. 5 (2010), No. 1, 65–80

# Start module
module MCJulia

import Base.*
export Blob, Sampler, sample, reset, flat_chain, save_chain

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
	iterations::Int64
	accepted::Int64
	args::(Any...)
end

# Constructor
function Sampler(k::Integer, dim::Integer, f::Function, a::Real, args::(Any...))
	accpt=0
	iter=0
	chain = zeros(Float64, (k, dim, 0))
	ln_p = zeros(Float64, (k, 0))
	blobs = Array(Blob, (k, 0))
	S = Sampler(k,dim,f,a,chain,ln_p,blobs,iter,accpt,args)
	return S
end

# Minimal constructors
Sampler(k::Integer, dim::Integer, f::Function, a::Real) = Sampler(k, dim, f, a, ())
Sampler(k::Integer, dim::Integer, f::Function, args::(Any...)) = Sampler(k, dim, f, 2.0, args)
Sampler(k::Integer, dim::Integer, f::Function) = Sampler(k, dim, f, 2.0, ())

call_lnprob(S::Sampler, pos::Array{Float64}) = S.probfn(pos, S.args...)

# Return lnprobs and blobs at given position
function get_lnprob(S::Sampler, pos::Array{Float64})
	lnprob = zeros(Float64, S.n_walkers)
#	blobs = Array(Blob, S.n_walkers)
#	b = false
	for k = 1:S.n_walkers
		p = call_lnprob(S, vec(pos[k,:]))
		if length(p) == 1
			lnprob[k] = p
		elseif length(p) == 2
			b = true
			lnprob[k] = p[1]
			blobs[k] = p[2]
		else
			error("Something weird came out of S.lnprobf()")
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
	divisions = [(first, second), (second, first)]

	for i = i0+1 : i0+N
		for ensembles in divisions
			active, passive = ensembles
			l_pas = length(passive)
			for k in active
				X_active = p[k,:]
				choice = passive[randi(l_pas)]
				X_passive = p[choice,:]
				z = randZ(S.a)
				proposal = X_passive + z * (X_active - X_passive)
				new_lnprob = call_lnprob(S, proposal)
				log_ratio = (S.dim - 1) * log(z) + new_lnprob - lnprob[k]
				if log(rand()) <= log_ratio
					lnprob[k] = new_lnprob
					p[k,:] = proposal
					S.accepted += 1
				end
				S.iterations += 1
				if storechain && (i-i0) % thin == 0
					S.ln_posterior[k, fld(i,thin)] = lnprob[k]
					S.chain[k, :, fld(i,thin)] = vec(p[k,:])
				end
			end
		end
	end
	return p
end

function sample(S::Sampler, N::Int64, thin::Int64, storechain::Bool)
	N = size(S.chain, 3)
	if N == 0
		error("No initial position for chain!")
	end
	p0 = S.chain[:,:,N]
	sample(S, p0, N, thin, storechain)
end

sample(S::Sampler, N::Int64) = sample(S, N, 1, true)
sample(S::Sampler, N::Int64, storechain::Bool) = sample(S, N, 1, storechain)


# Reset a Sampler to state after construction.
function reset(S::Sampler)
	k = S.n_walkers
	S.chain = zeros(Float64, (k, S.dim, 0))
	S.ln_posterior = zeros(Float64, (k, 0))
	S.blobs = Array(Blob, (k, 0))
	S.accepted = 0
	S.iterations = 0
	return S
end

# Flatten the chain along the walker axis
function flat_chain(S::Sampler)
	walkers,dimensions,steps = size(S.chain)
	flatchain = zeros((dimensions, walkers*steps))
	for step = 1:steps
		k = (step-1)*walkers + 1
		for dim = 1:dimensions
			flatchain[dim, k:(k+walkers-1)] = S.chain[:,dim,step]
		end
	end
	return flatchain
end

# Squash the chains and save them in a csv file
function save_chain(S::Sampler, filename::String)
	flatchain = flat_chain(S)
	csvwrite(filename, flatchain)
end


end # module ends
