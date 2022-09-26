include("../src/DivergenceLCE.jl")
include("../src/LCEProblem.jl")
using Plots

Δt = 0.05
x = trajectory(Systems.lorenz(),1000.0; Δt )[:,1]
lyapunovspectrum(Systems.lorenz(), 8000)

dim = 4
tau = 7
em = embed(x,dim,tau)

em = EmbeddedData(x, em,dim,tau)
prob = LCEProblem(em, Δt)
alg = DivergenceAlgorithm(lyapspan = 0:4:200, ntype = NeighborNumber(5))
sol = solve(prob,alg)