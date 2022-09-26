struct EmbeddedData
    timeseries::Vector
    embedded_data::Dataset
    dim::Number
    tau::Number
end

struct LCEProblem
    timeseries::Array
    timestep::Union{Number,Nothing}
    embedded_data::Union{EmbeddedData,Nothing}
end

function LCEProblem(em::EmbeddedData, step::Number)
    ts = em.timeseries
    LCEProblem(ts,step,em)
end


struct LCEMaxSolution
    maxLCE::Number
    LCEprogression::Array
end

struct LCESpectrumSolution
    LCEs::Array
end