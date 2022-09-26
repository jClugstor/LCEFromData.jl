using UnicodePlots

abstract type LCEAlgorithm end

struct EmbeddedData
    timeseries::Vector
    embedded_data::Dataset
    dim::Number
    tau::Number
end

struct LCEProblem
    timeseries::Vector
    timestep::Number
    embedded_data::EmbeddedData
end


function LCEProblem(em::EmbeddedData, step::Number)
    ts = em.timeseries
    LCEProblem(ts,step,em)
end



struct LCEMaxSolution
    maxLCE::Number
    LCEprogression::Array
    algorithm::LCEAlgorithm
end

function Base.show(io::IO,sol::LCEMaxSolution)
    println(io, "Solved using $(typeof(sol.algorithm))")
    println(io, "Max LCE: $(sol.maxLCE)")
    display(lineplot(sol.LCEprogression))
    sol
end


struct LCESpectrumSolution
    LCEs::Array
end


function solve end