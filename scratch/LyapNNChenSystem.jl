include("../src/LyapWithNN.jl")

plotlyjs()

function ChenSystem(du,u,p,t)
    α, b, γ = p

    x = u[1] 
    y = u[2]
    z = u[3]

    dx = α*(y - x)
    dy = (γ - α)*x - x*z + γ*y
    dz = x*y - b*z

    du[1] = dx
    du[2] = dy
    du[3] = dz
end

Span = (0,100)
y0 = [1.0, 1.0, 1.0]
pars = [35,3,27]
q = [1,1,1]
stepsize = 0.01
prob = FODESystem(ChenSystem,q,y0,Span,pars)
Ys = FractionalDiffEq.solve(prob, stepsize, PECE())
x = Ys.u[1,:]
Y, τ_vals, ts_vals, Ls, ϵs = pecuzal_embedding(x)
indat, outdat = create_training_data(Y,τ_vals[end])

indat = indat[begin:10:end]
outdat = outdat[begin:10:end]

model = train(indat, outdat, 500)

jac = Flux.jacobian(model, [1,2,3,4,5])

plot(reduce(vcat,model.(indat)), label = "Prediction", seriestype = :scatter)
plot!(reduce(vcat,outdat), label = "Actual", seriestype = :scatter)



indat = collect(0:0.4:20)
outdat = sin.(indat)

plot(indat,reduce(vcat,model.([[dat] for dat in indat])), label = "Prediction", seriestype = :scatter)
plot!(indat,outdat, label = "Actual", seriestype = :scatter)
