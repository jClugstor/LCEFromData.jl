prob = LCEProblem()





fname = "./src/Wolf Lyapunov/Data2.lor"
datcnt = 16384;
tau = 10;
ndim = 3;
ires = 10;
maxbox = 6000;
db = basgen(fname, tau, ndim, ires, datcnt, maxbox);
dt = .01;
evolve = 20;
dismin = 0.001;
dismax = 0.3;
thmax = 30;