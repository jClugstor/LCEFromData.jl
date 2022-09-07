#Jadon Clugston
#8/5/2022
#Implementing the method for calculating lyapunov spectrum from time series using Neural Networks
#From Lyudmila A. Dimtrieva et. al

using FractionalDiffEq, Flux, DynamicalSystems, LinearAlgebra

#need to collect data and set up the training set 
function create_training_data(data::Dataset,τ)
    input_set = []
    output_set = []
    T = length(data)

    for i in 1:(T - τ)
        push!(input_set,data[i])
        push!(output_set,[data[i+1][end]])
    end
    input_set,output_set
end


function build_model(input_length, output_length, hidden)
    return Flux.Chain(
        Flux.Dense(input_length => hidden,tanh),
        Flux.Dense(hidden => hidden, tanh),
        Flux.Dense(hidden => output_length))
end

function train(in_vals, out_vals, iters)
    hidden = 50

    input_length = length(in_vals[1])
    output_length = length(out_vals[1])

    model = build_model(input_length, output_length, hidden)
    loss(x,y) = Flux.Losses.mse(model(x), y)
    dat = zip(in_vals,out_vals)
    opt = Descent(0.01)
    

    Flux.@epochs iters my_custom_train!(loss,Flux.params(model),dat, opt) 
    model
end


function my_custom_train!(loss, ps, data, opt)
    # training_loss is declared local so it will be available for logging outside the gradient calculation.
    local training_loss
    ps = Flux.Params(ps)
    for d in data
      gs = gradient(ps) do
        training_loss = loss(d...)
        # Code inserted here will be differentiated, unless you need that gradient information
        # it is better to do the work outside this block.
        return training_loss
      end
      # Insert whatever code you want here that needs training_loss, e.g. logging.
      # logging_callback(training_loss)
      # Insert whatever code you want here that needs gradients.
      # e.g. logging histograms with TensorBoardLogger.jl to check for exploding gradients.
      Flux.update!(opt, ps, gs)
      # Here you might like to check validation set accuracy, and break out to do early stopping.
    end
end

function Ghat(invec, embedded_dat, model)
    inlen = length(invec)
    ind = nothing
    for (index,point) in enumerate(embedded_dat)
        if point == invec
           ind = index 
        end
    end 
    if isnothing(ind)
        error("Supplied vector is not on the trajectory")
    end
    out = embedded_dat[ind+1][1:(inlen-1)]
    out = vcat(out, model(invec))
end

function DGhat(invec, model)
    inlen = length(invec)
    mat = [Int(y == (x+1)) for x = 1:(inlen-1), y = 1:inlen]
    jac = Flux.jacobian(model,invec)[1]
    mat = vcat(mat, jac)
end


function lyap_spectrum(invecs, model)
    inlen = length(invecs)
    finalR = nothing
    qrs = []
    

    push!(qrs,qr(DGhat(invecs[1],model)))
    for l in 2:(inlen - 1)
        push!(qrs,qr(DGhat(invecs[l],model)*qrs[l-1].Q))
    end


    rs = [qrpair.R for qrpair in qrs]
    #finalR = foldl(*,[mat.R for mat in qrs])
    (1/inlen)*sum(x-> log.(abs.(diag(x))), rs)
end





