"""Background correction for a 2D data array
Args:
    data (2D array): channel data
    type (Background): type of background correction, i.e. plane_linear_fit, subtract_minimum line_average, vline_average, line_linear_fit, vline_linear_fit.
    offset (bool): if true (default), then the data-array will be shifted such that its minimum is 0.
Returns:
    data (2D array): background corrected data.
"""
function correct_background(data::Array{<:Number,2}, type::Background, offset::Bool=true)
    #returns Array{<:Number,2}

    if type == no_correction
        return data
    end

    data = copy(data)
    if type == plane_linear_fit  # subtract plane
        # see: https://math.stackexchange.com/a/2306029
        ci = vec(CartesianIndices(data))
        ci1 = [c[1] for c in ci]
        ci2 = [c[2] for c in ci]
        X = [ones(length(ci)) ci1 ci2]
        y = vec(data)

        not_nan = findall(!isnan, y)  # we need to skip the NaN values, otherwise we get NaN results (using "missing" didn't help either)
        @views @inbounds p = X[not_nan,:] \ y[not_nan]

        data -= reshape(X * p, size(data))
    elseif type == line_average  # subtract line by line average
        data = data .- mean(data, dims=2)
    elseif type == vline_average  # subtract line by line average for vertical lines (i.e. slow scan direction)
        # we want to account for NaNs here
        # return data .- mean(data, dims=1)
        for y in eachcol(data)
            not_nan = findall(!isnan, y)
            if length(not_nan) > 0
                @views @inbounds y .-= mean(y[not_nan])
            end
        end
    elseif type == line_linear_fit_legacy
        # efficient linear fit: https://discourse.julialang.org/t/efficient-way-of-doing-linear-regression/31232/26
        x = 1:size(data)[2]
        for y in eachrow(data)
            X = [ones(size(x)) x]

            α, β = X \ y
            y .= y .- α - x .* β
        end
    elseif type == vline_linear_fit_legacy
        # as line_linear_fit_, but for columns (and here we have to account for NaN values)
        x = 1:size(data)[1]
        for y in eachcol(data)
            X = [ones(size(x)) x]

            not_nan = findall(!isnan, y)
            if length(not_nan) > 0
                @views @inbounds α, β = X[not_nan,:] \ y[not_nan]
                y .= y .- α - x .* β
            end
        end
    elseif type == line_linear_fit
        # https://en.wikipedia.org/wiki/Ordinary_least_squares#Simple_linear_regression_model
        # this seems to be faster than the method above
        x = 1:size(data)[2]
        varx = var(x)  # same for each row
        meanx = mean(x)
        for y in eachrow(data)
            β = cov(x, y) / varx
            α = mean(y) - β * meanx
            y .= y .- α - x .* β
        end
    elseif type == vline_linear_fit
        # as line_linear_fit, but for columns (and here we have to account for NaN values)
        x = 1:size(data)[1]
        for y in eachcol(data)
            not_nan = findall(!isnan, y)
            if length(not_nan) > 0
                @views @inbounds β = cov(x[not_nan], y[not_nan]) / var(x[not_nan])
                @views @inbounds α = mean(y[not_nan]) - β * mean(x[not_nan])
                y .= y .- α - x .* β
            end
        end
    end
    if type == subtract_minimum || offset
        filtered = filter(!isnan, data)
        if length(filtered) > 0
            m = minimum(filtered)
            data .-= m
        end
    end
    return data
end
