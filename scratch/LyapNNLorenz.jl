include("../src/LyapWithNN.jl")
using Plots
plotlyjs()

A = trajectory(Systems.lorenz(),100)
lyapunovspectrum(Systems.lorenz(), 8000)


Y = embed(A[5000:end,1],3, 10)
plot(Y[:,1], Y[:,2], seriestype = :scatter)


indat, outdat = create_training_data(Y,10)
train_indat = indat[begin:10:end]
train_outdat = outdat[begin:10:end]

model = train(train_indat, train_outdat,200)

plot(reduce(vcat,model.(train_indat)), label = "Prediction", seriestype = :scatter)
plot!(reduce(vcat,train_outdat), label = "Actual", seriestype = :scatter)




lyap_spectrum(indat,model)