using UnicodePlots

abstract type LCEAlgorithm end
export LCEAlgorithm

struct EmbeddedData
    timeseries::AbstractArray
    embedded_data::StateSpaceSet
    dim::Number
    tau::Number
end

function EmbeddedData(timeseries, dim, tau)
    embedded_data = DynamicalSystems.embed(timeseries,dim,tau)
    EmbeddedData(timeseries,embedded_data,dim,tau)
end

export EmbeddedData

struct LCEProblem
    embedded_data::EmbeddedData
end

function LCEProblem(timeseries, tau, dim)
    em = EmbeddedData(timeseries, tau, dim)
    LCEProblem(em)
end

export LCEProblem



export LCEProblem


struct LCEMaxSolution
    maxLCE::Number
    LCEprogression::Array
    algorithm::LCEAlgorithm
end

export LCEMaxSolution

function Base.show(io::IO,sol::LCEMaxSolution)
    println(io, "Solved using $(typeof(sol.algorithm))")
    println(io, "Max LCE: $(sol.maxLCE)")
    display(lineplot(sol.LCEprogression))
end


struct LCESpectrumSolution
    LCEs::Array
    LCEprogression::Array
    algorithm::LCEAlgorithm
end

export LCESpectrumSolution

function solve end
export solve