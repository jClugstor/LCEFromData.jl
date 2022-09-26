includet("../src/LCEProblem.jl")
includet("../src/Wolf Lyapunov/WolfLyapunov.jl")

x = readlines("src/Wolf Lyapunov/Data2.lor")
#data = zeros(1,datcnt)
data = parse.(Float64,x)
#this will change

tau = 10;
ndim = 3;
em = embed(data,ndim,tau)
em = EmbeddedData(data,em,ndim,tau)
prob = LCEProblem(em,0.01)

ires = 10;
maxbox = 6000;
dt = 0.01;
evolve = 20;
dismin = 0.001;
dismax = 0.3;
thmax = 30;

w = WolfAlgorithm(ires = ires, maxbox = maxbox, dt = dt, evolve = evolve, dismin = dismin, dismax = 0.3, thmax = thmax)
sol = solve(prob,w)
plot(sol.LCEprogression)