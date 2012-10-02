# MC Julia

MC Julia is an affine-invariant ensemble MCMC sampler, written in the [Julia programming language](http://julialang.org). It implements a highly efficient MCMC algorithm by [Goodman & Weare (2010)](http://msp.org/camcos/2010/5-1/p04.xhtml), and is heavily inspired by the [`emcee`](http://danfm.ca/emcee/) Python package.

The algorithm is very efficient at sampling probability distributions which are highly skewed but more or less linearly correlated.

# Usage


The algorithm is implemented through the `Sampler` type. It keeps track of the sampled function and some optional parameters, and stores the chain of samples and other sampling output.

See the `examples` directory for some example programs.

## Constructors

#### Sampler(k, dim, f)
This is the minimal constructor for a new Sampler object. The parameter `k` gives the number of walkers in the ensemble, `dim` gives the dimension of the sampling space and `f` is the log-probability function. The number of walkers must be even and greater than twice the dimension.

The log-probability function must be defined to take the sampled position as an Array{Float64}, even when the sampling space is one-dimensional. 

```
log_probability(X::Array{Float64}) = -(X[1] - X[2])^2 / 0.5 - (X[1] + X[2])^2 / 1.5

S = Sampler(50, 2, log_probability)
```

#### Sampler(k, dim, f, a, args)
#### Sampler(k, dim, f, args)
#### Sampler(k, dim, f, a)

It is also possible to give some additional parameters to the sampler. The scale parameter `a::Float64` defines the stretch move scale (see Goodman & Weare, 2010). Its default value is 2.0, which is good for most cases, but tuning the parameter may give a performance increase in some situations. The parameter `args` is an arbitrary tuple which is passed to the log-probability function along with the position vector. The args tuple is expanded in the function call, so a tuple `(a,b,c)` will lead to a function call `f(X, a, b, c)`.

The following is a simple example of constant parameters with a one-dimensional Gaussian:

```
log_probability(X::Array{Float64}, mu::Float64, sigma::Float64) = -(X[1] - mu)^2 / (2 * sigma^2)

args = (1.0, 0.5)
S = Sampler(50, 2, log_probability, args)
```


## Functions

#### p = sample(S, p0, N, thin, storechain)
#### p = sample(S, p0, N, storechain)
#### p = sample(S, p0, N)

This function does the actual sampling. It takes as parameters a `Sampler` object `S`, an initial position `p0`, the number of samples `N`, the thinning factor `thin` and a boolean `storechain`.

The initial position `p0` must be an `Array{Float64, 2}` with the shape `(k, dim)`, and give the initial position in the `dim`-dimensional sampling space for each of the `k` walkers in the ensemble.

The `Int64` parameter `thin` gives thinning factor, the interval by which the positions of the walkers are saved to `S.chain`. That is, with `thin=10`, the chain is updated once every ten samping steps.

The boolean parameter `storechain` tells whether to save the steps of the walkers to the chain. Setting it to false will produce a `burn-in`, where the walkers move around, settling towards their steady state around the high-probability region of the sampled density, without any of the steps being saved.

The return value `p` is an `Array{Float64, 2}` giving the positions of the walkers at the last step of the sampling. It can directly be used as the initial position of another call to `sample`. It is useful to first do a burn-in with `storechains` set to `false`, then perform the actual sampling starting from the positions returned by the first call.

#### reset(S)

Reset a `Sampler` object to the state it was in right after construction. This will set `S.iterations` and `S.accepted` to zero and also empty the `S.chain`, `S.ln_posterior` and `S.blobs` arrays.

#### flat_chain(S)

Return the `S.chain` array flatted along the walker dimension.

#### save_chain(S, filename)

Save the flattened chain into a plaintext file in CSV format.


## Fields

These are the fields defined for the `Sampler` type. The values of some of them will be of interest to the user. The user should not need to change the values of any fields manually.

#### Sampler.n_walkers
The number of walkers in the ensemble. Set by the user through the constructor.

#### Sampler.a
The stretch move scale parameter. Set by the user through the constructor.

#### Sampler.dim
The dimension of the sampling space. Set by the user through the constructor.

#### Sampler.probfn
The log-probability function to be sampled. Set by the user through the constructor.

#### Sampler.args
The constant arguments for the log-probability function. Set by the user through the constructor.

#### Sampler.chain
An `Array{Float64, 3}` with shape `(k, dim, N)`, where `k` is the number of walkers in the ensemble, `dim` is the dimension of the sampling space and `N` is the number of samples taken (see the `sample` function, below). This is what you are mostly interested in acquiring.

#### Sampler.lnprob

An `Array{Flaot64, 2}` with shape `(k, N)`, giving the log-probability value of each walker at each point in the chain.

#### Sampler.iterations
An integer giving the number of moves proposed since the last reset (or creation) of the `Sampler` object. Includes also sampling steps done with `storechain == false`. One sampling step with `k` walkers will add `k` to this variable.

#### Sampler.accepted
An integer giving the number of accepted moves since the last reset (or creation) of the `Sampler` object. Includes also sampling steps done with `storechain == false`. Every accepted move of a walker will add one to this counter. With `Sampler.iterations` can be used to estimate the acceptance fraction of the full ensemble.




# Future work
### Parallelization
This is definitely planned as the next major improvement.

### Blobs
The Python package `emcee` allows the log-probability function to return arbitrary metadata "Blobs" which are also saved for each step in the chain. This may be a useful feature and will be implemented after parallelization.