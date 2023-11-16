using UnicodePlots

abstract type LCEAlgorithm end
export LCEAlgorithm

struct EmbeddedData
    timeseries::AbstractArray
    embedded_data::StateSpaceSet
    dim::Number
    tau::Number
end
export EmbeddedData

struct LCEProblem
    timeseries::Vector
    timestep::Number
    embedded_data::EmbeddedData
end
export LCEProblem

function LCEProblem(em::EmbeddedData, step::Number)
    ts = em.timeseries
    LCEProblem(ts,step,em)
end

function LCEProblem(ts::AbstractArray, step::Number)

end

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