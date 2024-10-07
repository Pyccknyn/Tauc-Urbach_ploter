using DataFrames
using CSV
using GLMakie
using GLM
using StatsBase

@kwdef struct DATA
    R::DataFrame
    T::DataFrame
end

@kwdef struct CALCULATIONS
    sample_name::String
    wavelengths::Vector
    absorbance::Vector{Float64}
    T100_minus_R::Vector{Float64}
    minus_ad::Vector{Float64}
    hv::Vector{Float64}
    hv_absorbance_y::Vector{Float64}
    urbach::Vector{Float64}
end

include("functions.jl")

begin
    config = read_config("config")
    γ = config["γ"]
    if config["onlyTdata"] == false
        dataT_name = config["dataT_name"]
        dataR_name = config["dataR_name"]
    else
        dataT_name = config["dataT_name"]
    end
end

begin

    nrowsD = Dict{String, Integer}()
    ncolsD = Dict{String, Integer}()

    if config["onlyTdata"] == false
        data = import_data(dataT_name, dataR_name)
        nrowsD = Dict("R" => nrow(data.R),"T" => nrow(data.T))
        ncolsD = Dict("R" => ncol(data.R),"T" => ncol(data.T))
    else
        data = import_data(dataT_name)
    end

    #nrowsD = Dict("R" => nrow(data.R),"T" => nrow(data.T))
    #ncolsD = Dict("R" => ncol(data.R),"T" => ncol(data.T))

    names_of_columns = list_names(data.T, ncolsD["T"])

    namesD = sort(Dict(i => (names_of_columns[i]) for i in 1:length(names_of_columns)))

    sample = CALCULATIONS[]

    for i in 1:length(namesD)

        wavelengths = get_wavelengths(data.R, i, nrowsD["R"])
        absorbance = calculate_absorbance(data.R, data.T, i, nrowsD["R"], nrowsD["T"])
        T100_minus_R = calculate_T100_minus_R(data.R, data.T, i, nrowsD["R"], nrowsD["T"])
        minus_ad = calculate_minus_ad(T100_minus_R)
        hv = calculate_hv(wavelengths, 1239.841)

        sample_calculations = CALCULATIONS(
            namesD[i],
            wavelengths,
            absorbance,
            T100_minus_R,
            minus_ad,
            hv,
            calculate_hv_absorbance_1_r(hv, absorbance, γ),
            calculate_urbach_energy(minus_ad)
        )
        push!(sample, sample_calculations)
    end
end

begin
    
    tauc = Figure[] 

    for i in 1:length(sample)

        tauc_plot = create_tauc_plot(sample[i].sample_name, sample[i].hv, sample[i].hv_absorbance_y, nrowsD, γ)

        push!(tauc, tauc_plot)

    end

    urbach = Figure[]
    
    for i in 1:length(sample)

        urbach_plot = create_urbach_plot(sample[i].sample_name, sample[i].hv, sample[i].urbach, nrowsD)

        push!(urbach, urbach_plot)
    end

end

