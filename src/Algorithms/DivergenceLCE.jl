struct DivergenceAlgorithm <: LCEAlgorithm
    ks
    #same as the keyword args in `lyapunov_from_data` in DynamicalSystems.jl, minus refstates
    w 
    distance
    ntype

    #same as the keyword args in `linear_region` in DynamicalSystems.jl
    dxi 
    tol
    ignore_saturation 

end

"""
    DivergenceAlgorithm(kwargs...)
    
Constructor for an object that holds hyperparameters for finding an LCE using the Divergence method.
Uses functions `lyapunov_from_data` and `linear_region` from [ChaosTools.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/chaostools/stable/lyapunovs/)

Keyword arguments: 

    - lyapspan: 

    - w = 1: passed to `lyapunov_from_data`, Theiler window, minimum time seperation that neighbors should have

    - distance = FirstElement(): passed to `lyapunov_from_data`, this specifies what distance function to use for the logarithmic distance

    - ntype = NeighborNumber(1): passed to `lyapunov_from_data`, the type of [neighborhood](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/chaostools/stable/lyapunovs/#Neighborhood.NeighborNumber) used

    - dxi = 1: passed to `LargestLinearRegion`

    - tol = 0.25: passed to `LargestLinearRegion`

    - ignore_saturation = true: passed to `slopefit`

"""
function DivergenceAlgorithm(ks; w = 1, distance = FirstElement(), ntype = NeighborNumber(1), dxi = 1, tol = 0.25, ignore_saturation = true)
    DivergenceAlgorithm(ks,w,distance,ntype,dxi,tol,ignore_saturation)
end

export DivergenceAlgorithm

function solve(prob::LCEProblem, alg::DivergenceAlgorithm)
    Y = prob.embedded_data.embedded_data
    ks = alg.ks
    stepsize = prob.dt
    
    e = lyapunov_from_data(Y,ks, w = alg.w, distance = alg.distance, ntype = alg.ntype)
    slope = slopefit(ks .* stepsize,e, LargestLinearRegion(tol=alg.tol,dxi = alg.dxi), ignore_saturation = alg.ignore_saturation)[1]
    LCEMaxSolution(slope,e,alg)
end

export solve
