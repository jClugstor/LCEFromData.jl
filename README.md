# LCEFromData

[![Build Status](https://github.com/jClugstor/LCEFromData.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jClugstor/LCEFromData.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a package with methods for finding the [Lyapunov Characteristic Exponents (LCEs)](https://en.wikipedia.org/wiki/Lyapunov_exponent) of a dynamical system given state data. 

This package has the goal of combining disparate software packages and methods for computing LCEs from data in to one place. It uses an API similar to that of the [SciML Common Interface](https://docs.sciml.ai/Overview/stable/highlevels/interfaces/) so that different methods can be used with the same interface.

# API
This package is organized around the idea of `LCEProblem`s and `LCEAlgorithm`s. An `LCEProblem` object holds the time series, the delay embedded data set to use, and the dimension and lag value of the dataset. `LCEAlgorithm` objects determine what method will be used, and holds solver parameters.

To use this package, first construct an `LCEProblem` using `LCEProblem(timeseries, dim, tau)`. This will construct an `LCEProblem` that holds a delay embedding of `timeseries` using the `dim` and `tau` parameters. Then construct an `LCEAlgorithm` using the desired parameters, e.g. `DivergenceAlgorithm(ks; w = 1, distance = FirstElement(), ntype = NeighborNumber(1), dxi = 1, tol = 0.25, ignore_saturation = true)`. This will result in an `LCESolution` which holds the solution information.  


# Methods

So far only two methods are available.
## Divergence Method
Measures the logarithm of the divergence of two nearby trajectories of an lag embedded data set. Uses methods from [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/dev/)

## Wolf Method
Based on Alan Wolf's paper

Alan Wolf, Jack B. Swift, Harry L. Swinney, John A. Vastano,
Determining Lyapunov exponents from a time series,
Physica D: Nonlinear Phenomena,
Volume 16, Issue 3,
1985,
Pages 285-317,
ISSN 0167-2789,
https://doi.org/10.1016/0167-2789(85)90011-9.
(https://www.sciencedirect.com/science/article/pii/0167278985900119)