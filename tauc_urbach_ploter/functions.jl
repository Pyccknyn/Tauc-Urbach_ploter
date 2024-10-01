function read_config(config::String)::Dict{String, Any}
    lines = readlines("$config.txt")
    configD = Dict{String, Any}()

    for line in lines
        key, value = split(line, "=")
        key = strip(key)
        value = strip(value)
        
        if value == "true" || value == "false"
            parsed_value = parse(Bool, value)
        elseif occursin(r"^-?\d+\.?\d*$", value)  # Regex to match Float64
            parsed_value = parse(Float64, value)
        else
            parsed_value = value  # Treat as String
        end

        configD[key] = parsed_value
    end

    return configD
end

function get_wavelengths(data::DataFrame, sample::Integer, nrows::Integer)::AbstractVector
    wavelengths = []
    if sample%2 == 0
        sample += 1
    end
    for i in 1:nrows
        push!(wavelengths, data[i, sample])
    end
    return wavelengths
end

function calculate_absorbance(dataR::DataFrame, dataT::DataFrame, sample::Integer, nrowsR::Integer, nrowsT::Integer)::Vector{Float64}
    absorbance = Float64[]  

    for i in 1:(nrowsR)  
        push!(absorbance, 100.0 - dataR[i, 2 * sample]) 
    end
    
    for i in 1:(nrowsT)
        absorbance[i] -= dataT[i, 2 * sample] 
    end
    
    return absorbance
end

function calculate_absorbance_onlyT(dataT::DataFrame, sample::Integer, nrowsT::Integer)::Vector{Float64}
    absorbance = Float64[]
    for i in 1:(nrowsT)
        if dataT[i, 2 * sample] <= 0
            push!(absorbance, NaN)
        else
            push!(absorbance, -log10(T))
        end
    end
    return absorbance
end

function calculate_T100_minus_R(dataR::DataFrame, dataT::DataFrame, sample::Integer, nrowsR::Integer, nrowsT::Integer)::Vector{Float64}
    T100_minus_R = Float64[] 

    for i in 1:(nrowsR)  
        push!(T100_minus_R, dataT[i, 2 * sample]) 
    end
    
    for i in 1:(nrowsT)
        T100_minus_R[i] = T100_minus_R[i]./(100.0 - dataR[i, 2 * sample]) 
    end
    
    return T100_minus_R
end

function calculate_minus_ad(T100_minus_R::Vector{Float64})::Vector{Float64}
    minus_ad = Float64[]
    for i in 1:length(T100_minus_R)
        if T100_minus_R[i] <= 0
            push!(minus_ad, NaN)
        else
            push!(minus_ad,log(ℯ,T100_minus_R[i]))
        end
    end
    return minus_ad
end

function calculate_hv(wavelengths::AbstractVector, factor::Float64)::Vector{Float64}
    hv = Float64[]  
    for i in 1:length(wavelengths)
        push!(hv, factor / wavelengths[i])  
    end
    return hv  
end

calculate_hv_absorbance_1_r(hv::Vector{Float64}, absorbance::Vector{Float64}, y::Number)::Vector{Float64} = (hv.*absorbance).^y

function calculate_urbach_energy(minus_ad::Vector{Float64})::Vector{Float64}
    urbach_energy = Float64[]
    for i in 1:length(minus_ad)
        if isnan(minus_ad[i])
            push!(urbach_energy, 0)
        else
            push!(urbach_energy, log(ℯ, -minus_ad[i]))
        end
    end
    return urbach_energy
end

function list_names(data::DataFrame, ncols::Integer)::Vector{String} 
    colNames = names(data)[1:2:ncols]
    colNames = map(x -> replace(x, ".asc [nm]" => ""), colNames)
    colNames = map(x -> replace(x, ".Cycle1" => ""), colNames)
    colNames = map(x -> replace(x, r"^[^0-9]*" => ""), colNames)
    return colNames
end

function create_tauc_plot(sample_name::String, hv::Vector{Float64}, hv_absorbance_y::Vector{Float64}, nrowsD::Dict{String, Int64}, γ::Number)::Figure
    fig = Figure()  # Create a new figure for the sample

    # Create a DataFrame for the current sample's data
    fig_data = DataFrame(X = hv, Y = hv_absorbance_y, Visible = trues(nrowsD["T"]))

    # Create axis for the figure
    ax = Axis(fig[1, 1], title = sample_name, xlabel = "hv [eV]", ylabel = "(αhv)^$γ [(eVcm⁻¹)^$γ]")
    
    # Set the x-axis limits based on the hv range of the current sample
    xlims!(ax, hv[begin], hv[end])

    #slider for selecting a range in hv values
    interrange = IntervalSlider(fig[2, 1], range = LinRange(hv[begin], hv[end], 10000), startvalues = (hv[begin], hv[end]))

    #variable for color update
    points = Point2(hv, hv_absorbance_y)

    ols = lm(@formula(Y ~ X), fig_data)

    #update slider text with the selected range values
    slidertext = lift(interrange.interval) do int
        string(round.(int, digits = 3))
    end

    Label(fig[2, 1], slidertext, tellwidth = false)

    #observables vertical line positions
    vl1_pos = Observable(hv[end])
    vl2_pos = Observable(hv[end])

    vl1 = vlines!(ax, vl1_pos, color = :blue)
    vl2 = vlines!(ax, vl2_pos, color = :blue)

    #observables for linear regression coefficients
    a_reglin = Observable(0.0)
    b_reglin = Observable(0.0)

    #update the vertical lines and perform regression on visible data
    lift(interrange.interval) do slider_dot
        vl1_pos[] = slider_dot[1]  
        vl2_pos[] = slider_dot[2]
    end

    # Observable for the x-intercept of the linear regression (Eg)
    urbach_energy = Observable(0.0)

    # Observable for R² value (goodness of fit)
    r_squared = Observable(0.0)

    visible_points_helperX = Observable([0.0])
    visible_points_helperY = Observable([0.0])

    lift(interrange.interval) do slider_dot 

        fig_data.Visible .= map(x -> slider_dot[1] <= x <= slider_dot[2], fig_data.X)

        #filter data used for linear regression
        visible_data = fig_data[fig_data.Visible, :]

        visible_points_helperX[] = visible_data[!, :X]
        visible_points_helperY[] = visible_data[!, :Y]

        if nrow(visible_data) > 1
            ols = lm(@formula(Y ~ X), visible_data)
        end
        
        a_reglin[] = coef(ols)[2]
        b_reglin[] = coef(ols)[1]
    
        r_squared[] = round(r2(ols), digits = 5)
        urbach_energy[] = round((-b_reglin[] / a_reglin[]), digits = 3)

    end
    
    #update colors of scattered points
    colors = lift(interrange.interval) do coloring
        map(points[1]) do p
            coloring[1] < p[1] < coloring[2]
        end
    end

    eg = lift(urbach_energy) do x
        Point2f(x,0)
    end

    #create scatter plot with color mapping
    s = scatter!(ax, points[1], points[2], color = colors, colormap = [(:red, 0.2), :red], strokewidth = 0, markersize = 7)

    Label(fig[1, 2], lift(urbach_energy, r_squared) do x,r2 "Eg: $x eV\nR²: $r2" end, width = 100, tellheight = false)

    #linear regression line
    reglin = ablines!(ax, b_reglin, a_reglin, color = :green)

    eg_val = scatter!(ax, eg, color = :black, markersize = 6)

    # Add a button to save the plot as a PNG
    save_button = Button(fig[2, 2], label = "Save")
    on(save_button.clicks) do _
        save("$(sample_name)_Tauc.png", fig, update = false)
        write("$(sample_name)_Tauc.txt", "x:$(visible_points_helperX[])\ny:$(visible_points_helperY[])
        \na_reglin:$(a_reglin[])\nb_reglin:$(b_reglin[])\nEg:$(urbach_energy[])\nR²:$(r_squared[])")
        println("Saved sample: $(sample_name)")
    end

    fig
end

function create_urbach_plot(sample_name::String, hv::Vector{Float64}, urbach::Vector{Float64}, nrowsD::Dict{String, Int64})::Figure
    fig = Figure()  # Create a new figure for the sample

    # Create a DataFrame for the current sample's data
    fig_data = DataFrame(X = hv, Y = urbach, Visible = trues(nrowsD["T"]))

    # Create axis for the figure
    ax = Axis(fig[1, 1], title = sample_name, xlabel = "hv [eV]", ylabel = "ln(-αd)")
    
    # Set the x-axis limits based on the hv range of the current sample
    xlims!(ax, hv[begin], hv[end])

    #slider for selecting a range in hv values
    interrange = IntervalSlider(fig[2, 1], range = LinRange(hv[begin], hv[end], 10000), startvalues = (hv[begin], hv[end]))

    #variable for color update
    points = Point2(hv, urbach)

    ols = lm(@formula(Y ~ X), fig_data)

    #update slider text with the selected range values
    slidertext = lift(interrange.interval) do int
        string(round.(int, digits = 3))
    end

    Label(fig[2, 1], slidertext, tellwidth = false)

    #observables vertical line positions
    vl1_pos = Observable(hv[end])
    vl2_pos = Observable(hv[end])

    vl1 = vlines!(ax, vl1_pos, color = :blue)
    vl2 = vlines!(ax, vl2_pos, color = :blue)

    #observables for linear regression coefficients
    a_reglin = Observable(0.0)
    b_reglin = Observable(0.0)

    #update the vertical lines and perform regression on visible data
    lift(interrange.interval) do slider_dot
        vl1_pos[] = slider_dot[1]  
        vl2_pos[] = slider_dot[2]
    end

    # Observable for the x-intercept of the linear regression (Eg)
    urbach_energy = Observable(0.0)

    # Observable for R² value (goodness of fit)
    r_squared = Observable(0.0)

    visible_points_helperX = Observable([0.0])
    visible_points_helperY = Observable([0.0])

    lift(interrange.interval) do slider_dot 

        fig_data.Visible .= map(x -> slider_dot[1] <= x <= slider_dot[2], fig_data.X)

        #filter data used for linear regression
        visible_data = fig_data[fig_data.Visible, :]

        visible_points_helperX[] = visible_data[!, :X]
        visible_points_helperY[] = visible_data[!, :Y]

        if nrow(visible_data) > 1
            ols = lm(@formula(Y ~ X), visible_data)
        end
        
        a_reglin[] = coef(ols)[2]
        b_reglin[] = coef(ols)[1]
    
        r_squared[] = round(r2(ols), digits = 5)
        urbach_energy[] = round((1000/a_reglin[]), digits = 3)

    end
    
    #update colors of scattered points
    colors = lift(interrange.interval) do coloring
        map(points[1]) do p
            coloring[1] < p[1] < coloring[2]
        end
    end

    eg = lift(urbach_energy) do x
        Point2f(x,0)
    end

    #create scatter plot with color mapping
    s = scatter!(ax, points[1], points[2], color = colors, colormap = [(:red, 0.2), :red], strokewidth = 0, markersize = 7)

    Label(fig[1, 2], lift(urbach_energy, r_squared) do e,r2 "Urbach_E: $e meV\nR²: $r2" end, width = 150, tellheight = false)

    #linear regression line
    reglin = ablines!(ax, b_reglin, a_reglin, color = :green)

    eg_val = scatter!(ax, eg, color = :black, markersize = 6)

    # Add a button to save the plot as a PNG
    save_button = Button(fig[2, 2], label = "Save")
    on(save_button.clicks) do _
        save("$(sample_name)_Urbach.png", fig, update = false)
        write("$(sample_name)_Urbach.txt", "x:$(visible_points_helperX[])\ny:$(visible_points_helperY[])
        \na_reglin:$(a_reglin[])\nb_reglin:$(b_reglin[])\nUrbach_E [meV]:$(urbach_energy[])\nR²:$(r_squared[])")
        println("Saved sample: $(sample_name)")
    end

    fig
end

function export_data(sample)
    write("$(sample.sample_name)_data.txt",
    "wavelengths: $(sample.wavelengths)\n
    absorbance: $(sample.absorbance)\n
    T100_minus_R: $(sample.T100_minus_R)\n
    minus_ad: $(sample.minus_ad)\n
    hv: $(sample.hv)\n
    hv_absorbance_y: $(sample.hv_absorbance_y)\n
    urbach: $(sample.urbach)\n"
    )
end