struct LCEProblem
    timeseries::Array
    timestep::Number
    embedded_data::EmbeddedData
end

struct EmbeddedData
    timeseries::Dataset
    embedded_data::Dataset
    dim::Number
    tau::Number
end

struct LCEMaxSolution
    maxLCE::Number
    LCEprogression::Array
end

struct LCESpectrumSolution
    LCEs::Array
end