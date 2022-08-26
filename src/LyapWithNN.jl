#Jadon Clugston
#8/5/2022
#Implementing the method for calculating lyapunov spectrum from time series using Neural Networks
#From Lyudmila A. Dimtrieva et. al

using FractionalDiffEq, Plots, Flux, DynamicalSystems
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
        Flux.Dense(hidden => hidden, tanh),
        Flux.Dense(hidden => output_length))
end

function train(in_vals, out_vals)
    hidden = 50

    input_length = length(in_vals[1])
    output_length = length(out_vals[1])

    model = build_model(input_length, output_length, hidden)
    loss(x,y) = Flux.Losses.mse(model(x), y)
    dat = zip(in_vals,out_vals)
    opt = Descent(0.0001)
    

    Flux.@epochs 2000 my_custom_train!(loss,Flux.params(model),dat, opt) 
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



Span = (0,25)
y0 = [1.0, 1.0, 1.0]
pars = [35,3,27]
q = [1,1,1]
stepsize = 0.01
prob = FODESystem(ChenSystem,q,y0,(0,5),pars)
Ys = FractionalDiffEq.solve(prob, stepsize, PECE())
x = Ys.u[1,:]
Y, τ_vals, ts_vals, Ls, ϵs = pecuzal_embedding(x)
indat, outdat = create_training_data(Y,τ_vals[end])

indat = indat[begin:2:end]
outdat = outdat[begin:2:end]

model = train(indat, outdat)

plot(reduce(vcat,model.(indat)), label = "Prediction", seriestype = :scatter)
plot!(reduce(vcat,outdat), label = "Actual", seriestype = :scatter)




indat = collect(0:0.4:20)
outdat = sin.(indat)

plot(indat,reduce(vcat,model.([[dat] for dat in indat])), label = "Prediction", seriestype = :scatter)
plot!(indat,outdat, label = "Actual", seriestype = :scatter)



