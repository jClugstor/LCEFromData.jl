include("../src/LyapWithNN.jl")
using Plots
plotlyjs()

A = trajectory(Systems.lorenz(),1000, Δt = 0.01)
lyapunovspectrum(Systems.lorenz(), 8000)


Y = embed(A[3000:end,1],3, 10)
plot(Y[:,1], Y[:,2], seriestype = :scatter)


indat, outdat = create_training_data(A,10)
train_indat = indat[begin:5:end]
train_outdat = outdat[begin:5:end]

model = train(train_indat, train_outdat,900)

plot(reduce(vcat,model.(train_indat)), label = "Prediction", seriestype = :scatter)
plot!(reduce(vcat,train_outdat), label = "Actual", seriestype = :scatter)

lyap_spectrum(indat,model)