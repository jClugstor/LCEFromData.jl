include("../src/LyapWithNN.jl")
using Plots


A = trajectory(Systems.henon(), 100000)
dat = A[:,1]
#plot(A[:,1], A[:,2], seriestype = :scatter)

Y = embed(A[:,1],2, 1)
#plot(Y[:,1], Y[:,2], seriestype = :scatter)


indat, outdat = create_training_data(Y,1)
train_indat = indat[begin:100:end]
train_outdat = outdat[begin:100:end]


model = train(train_indat, train_outdat)


#plot(reduce(vcat,model.(indat)), label = "Prediction", seriestype = :scatter)
#plot!(reduce(vcat,outdat), label = "Actual", seriestype = :scatter)
lyapunovspectrum(Systems.henon(),10000)
lyap_spectrum(indat,model)


lce_error(nnres, dynsysres) = abs.((nnres .- dynsysres) ./ dynsysres)
lce_error(lyap_spectrum(indat,model), lyapunovspectrum(Systems.henon(),10000))