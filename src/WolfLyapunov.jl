using DynamicalSystems, LinearAlgebra


function wolf_lyapunov(data::Vector, SCALMX::Int, SCALMN::Int, EVOLV::Int; ANGLMX = 0.2)
    Y, τ_vals, ts_vals, Ls, ϵs = pecuzal_embedding(data)

    #big loop
    for point in Y
        
    end
end