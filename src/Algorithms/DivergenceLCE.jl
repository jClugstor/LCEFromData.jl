struct DivergenceAlgorithm <: LCEAlgorithm
    lyapspan
    #same as the keyword args in `lyapunov_from_data` in DynamicalSystems.jl, minus refstates
    w 
    distance
    ntype

    #same as the keyword args in `linear_region` in DynamicalSystems.jl
    dxi 
    tol
    ignore_saturation 

end

function DivergenceAlgorithm(;lyapspan, w = 1, distance = FirstElement(), ntype = NeighborNumber(1), dxi = 1, tol = 0.25, ignore_saturation = true)
    DivergenceAlgorithm(lyapspan,w,distance,ntype,dxi,tol,ignore_saturation)
end

export DivergenceAlgorithm

function solve(prob::LCEProblem, alg::DivergenceAlgorithm)
    Y = prob.embedded_data.embedded_data
    lyapspan = alg.lyapspan
    stepsize = prob.timestep
    
    e = lyapunov_from_data(Y,lyapspan, w = alg.w, distance = alg.distance, ntype = alg.ntype)
    slope = linear_region(lyapspan .* stepsize,e,tol=alg.tol,dxi = alg.dxi, ignore_saturation = alg.ignore_saturation)[2]
    LCEMaxSolution(slope,e,alg)
end

export solve
