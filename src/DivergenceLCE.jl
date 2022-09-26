include("LCEProblem.jl")

Base.@kwdef struct DivergenceAlgorithm
    lyapspan
    #same as the keyword args in `lyapunov_from_data` in DynamicalSystems.jl, minus refstates
    w = 1
    distance = Euclidean()
    ntype = NeighborNumber(1)

    #same as the keyword args in `linear_region` in DynamicalSystems.jl
    dxi = 1
    tol = 0.25
    ignore_saturation = true

end


function solve(prob::LCEProblem, alg::DivergenceAlgorithm)

    Y = prob.embedded_data.embedded_data
    lyapspan = alg.lyapspan
    stepsize = prob.timestep
    
    e = lyapunov_from_data(Y,lyapspan, w = eval(alg.w), distance = eval(alg.distance), ntype = eval(alg.ntype))
    slope = linear_region(lyapspan .* stepsize,e,tol=alg.tol,dxi = alg.dxi, ignore_saturation = alg.ignore_saturation)[2]
    LCEMaxSolution(slope,e)
end

