include("../src/Wolf Lyapunov/WolfLyapunov.jl")

x = readlines("src/Wolf Lyapunov/Data2.lor")
#data = zeros(1,datcnt)
data = parse.(Float64,x)
#this will change
prob = LCEProblem(data,nothing, nothing)





fname = "./src/Wolf Lyapunov/Data2.lor"
datcnt = 16384;
tau = 10;
ndim = 3;
ires = 10;
maxbox = 6000;
dt = 0.01;
evolve = 20;
dismin = 0.001;
dismax = 0.3;
thmax = 30;

w = WolfAlgorithm(fname,tau,ndim,ires,datcnt,maxbox,dt,evolve,dismin,dismax,thmax)
solve(prob,w)