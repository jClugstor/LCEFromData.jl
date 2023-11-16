struct LyapNNAlgorithm <: LCEAlgorithm
    hidden
    loss
end

function solve(prob::LCEProblem, alg::NNLyapAlgorithm; plot = true)
    timestep = prob.timestep
    Y = prob.embedded_data.embedded_data
    τ = prob.embedded_data.tau
    embedm = prob.embedded_data.dim

    input_set = []
    output_set = []
    T = length(Y)
    
    for i in 1:(T - τ)
        push!(input_set,Y[i])
        push!(output_set,[Y[i+1][begin]])
    end

    dat = [(x,y) for (x,y) in zip(input_set, _output_set)]

    dat = shuffle(dat)


    model = Flux.Chain(
            Flux.Dense(input_length => hidden,tanh,init=Flux.glorot_normal),
            Flux.Dense(hidden => output_length,init=Flux.glorot_normal))

    opt_state = Flux.setup(AMSGrad(), model)

    loss(mod,x,y) = Flux.Losses.mse(mod(x), y)
            
            
    meaner = 1
    i = 0
    while meaner > alg.loss
        i = i+1
        Flux.train!(loss, model, dat, opt_state)
        meaner = mean([loss(model,x...) for x in dat])
        if i%10 == 0
            println(i)
            println(meaner)
        end
    end
            
    dg = only.([Flux.jacobian(x -> model(x),invec) for invec in scaled_input_set])
    dg = reduce(vcat,dg)

    mat = [Int(y == x) for x = 1:(embedm-1), y = 1:(embedm)]

    rs = []
    nmats = []
    Q0 = I
    N = length(dg[:,1])

    LLs = []
    for i = 1:embedm:N
        nmat = vcat(dg[i,:]',mat)
        push!(nmats,nmat)
        Q,R = qr(nmat*Q0)
        push!(rs,R)
        Q0 = Q
        push!(LLs,1/(timestep*τ)*mean(x-> log.(abs.(diag(x))), rs))
    end

    
    LCESpectrumSolution((1/(timestep*τ))*mean(x-> log.(abs.(diag(x))), rs))
end






