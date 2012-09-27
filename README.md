# MC Julia

MC Julia is an affine-invariant ensemble MCMC sampler, written in the [Julia programming language](http://julialang.org). It implements a highly efficient MCMC algorithm by [Goodman & Weare (2010)](http://msp.org/camcos/2010/5-1/p04.xhtml), and is heavily inspired by the [`emcee`](http://danfm.ca/emcee/) Python package.

The algorithm is very efficient at sampling probability distributions which are highly skewed but more or less linearly correlated.

# Usage

TODO: detailed instructions. In the meantime, look at the `test.jl` file for a simple example.