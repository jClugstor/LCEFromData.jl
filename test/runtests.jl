using LCEFromData
using Test

using DynamicalSystems

@testset "LCEFromData.jl" begin
    # Write your tests here.
    function LorenzSystem!(du, u, p, t)
        σ, β, ρ = p
        x, y, z = u
        du[1] = σ * (y - x)
        du[2] = x * (ρ - z) - y
        du[3] = x * y - β * z
    end

    pars = [16.0, 4.0, 45.92]
    tspan = (0, 500)
    u0 = [0.1, 0.1, 0.1]

    lorenz_prob = CoupledODEs(LorenzSystem!,u0, pars)
    Y,t = trajectory(lorenz_prob, 500.0, Δt = 0.1)

    x_dat = Y[:,1]
    prob = LCEProblem(x_dat,0.1, 3, 3)
    sol = solve(prob, DivergenceAlgorithm(0:1:100))
    @test sol isa LCEMaxSolution

    sol = solve(prob, WolfAlgorithm())
    @test sol isa LCEMaxSolution
end

function LorenzSystem!(du, u, p, t)
    σ, β, ρ = p
    x, y, z = u
    du[1] = σ * (y - x)
    du[2] = x * (ρ - z) - y
    du[3] = x * y - β * z
end

pars = [16.0, 4.0, 45.92]
tspan = (0, 500)
u0 = [0.1, 0.1, 0.1]

lorenz_prob = CoupledODEs(LorenzSystem!, u0, pars)
Y, t = trajectory(lorenz_prob, 500.0, Δt=0.1)

x_dat = Y[:, 1]
prob = LCEProblem(x_dat, 0.1, 3, 3)
sol = solve(prob, DivergenceAlgorithm(0:1:100))
@test sol isa LCEMaxSolution

sol = solve(prob, WolfAlgorithm())
@test sol isa LCEMaxSolution