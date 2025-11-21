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
    elseif type == line_diff_mean
        # line correction using differences between adjacent lines, using the mean of the differences
        prev = @view data[begin, :]
        for y in eachrow(@view data[begin+1:end, :])
            diff = filter(!isnan, y .- prev)
            if length(diff) > 0
                y .-= mean(diff)
            end
            prev = y
        end
    elseif type == line_diff_median
        # line correction using differences between adjacent lines, using the median of the differences
        prev = @view data[begin, :]
        for y in eachrow(@view data[begin+1:end, :])
            diff = filter(!isnan, y .- prev)
            if length(diff) > 0
                y .-= median(diff)
            end
            prev = y
        end
    elseif type == plane_facets
        data, _, _ = plane_correction_normals(data)
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


"""plane correction using the calculated normals for each subsquare (with size `subsize` and distance `substep`), 
then subtract the plane corresponding to a weighted average of the normals"""
function plane_correction_normals(data::Array{<:Number,2}; sigma=0, subsize=9, step=7)
    data_filtered = imfilter(data, Kernel.gaussian(sigma))

    rows, cols = size(data_filtered)

    rrows = 1:step:(rows - subsize + 1)
    rcols = 1:step:(cols - subsize + 1)
    normals = fill([0.0, 0.0, 0.0], length(rrows), length(rcols))
    for (idx_i, i) in enumerate(rrows)
        for (idx_j, j) in enumerate(rcols)
            subdata = data_filtered[i:(i + subsize - 1), j:(j + subsize - 1)]
            subdata .-= mean(subdata)  # remove mean to improve numerical stability

            # skip patches with non-finite values
            if any(!isfinite, subdata)
                continue
            end

            # Fit a plane to the subimage
            # Make columns: ones, x, y so coeffs[2] is x-coeff, coeffs[3] is y-coeff.
            A = [ones(subsize^2) vec(repeat((1:subsize)', subsize)) vec(repeat(1:subsize, subsize))]
            b = vec(subdata)
            coeffs = A \ b  # Plane: z = coeffs[1] + coeffs[2]*x + coeffs[3]*y

            normal = [-coeffs[2], -coeffs[3], 1]
            normal /= sqrt(coeffs[2]^2 + coeffs[3]^2 + 1) # normalize

            # force normals to same hemisphere (nz >= 0)
            if normal[3] < 0
                normal .= -normal
            end

            normals[idx_i, idx_j] = normal
        end
    end

    if isempty(normals)
        @warn "plane_correction_normals: no valid subsquares; returning original image"
        return data, Matrix{Float64}(), Matrix{Float64}()
    end

    # Average the normals, weighted by their z-component
    # weights = [exp((-1.0 - n[3]) * 0.51) for n in normals]
    avg_normal = zeros(3)
    sum_weights = 0.0
    weights = zeros(size(normals, 1), size(normals, 2))
    for i in 2:size(normals, 1)-1, j in 2:size(normals, 2)-1
        n = normals[i, j]
        # stddev = std(normals[(i-1):(i+1), (j-1):(j+1)][3])
        # w = 1 / (stddev + 1) * exp(-abs((n[3] - 1.0)) * 0.5)  # custom weight: prefer normals with nz close to 1 and low stddev
        w = exp(-abs((n[3] - 1.0)) * 5)  # custom weight: prefer normals with nz close to 1 and low stddev
        weights[i, j] = w
        sum_weights += w
        avg_normal .+= (w .* n)
    end
    avg_normal ./= sum_weights

    # subtract the plane defined by the average normal
    cx = (cols + 1) / 2
    cy = (rows + 1) / 2
    d = -avg_normal[1] * cx - avg_normal[2] * cy

    # make X and Y arrays (Float64) that broadcast to the image shape (rows x cols)
    X = reshape(Float64.(1:cols), 1, cols)   # 1 x cols (x = column index)
    Y = reshape(Float64.(1:rows), rows, 1)   # rows x 1 (y = row index)

    plane = -(avg_normal[1] .* X .+ avg_normal[2] .* Y .+ d) ./ avg_normal[3]

    return data .- plane, normals, weights
end


"""calculates histogram"""
function hist(v, nbins; return_indices=false)
    mn, mx= extrema(v)
    edges = collect(range(mn, stop=mx, length=nbins+1))
    counts = zeros(Int, nbins)
    denom = mx - mn
    if return_indices
        indices = fill(Int[], nbins)
    else
        indices = nothing
    end 
    for val in v
        t = (val - mn) / denom  # in [0,1]
        # map t to bin index 1..nbins, ensure mx falls into last bin
        idx = clamp(floor(Int, t * nbins) + 1, 1, nbins)
        counts[idx] += 1
        if return_indices
            push!(indices[idx], val)
        end
    end
    return counts, edges, indices
end


"""Sets the 0-level to the largest plane in the image"""
function set_baselevel(data::Array{<:Number,2}, binwidth=0.2)
    v = vec(data)
    # drop NaNs/Infs
    v = v[isfinite.(v)]
    if isempty(v)
        @warn "set_baselevel: no finite pixels found; returning original image"
        return data
    end

    mn, mx = extrema(v)

    # if nearly constant, subtract that constant
    if isapprox(mn, mx; atol=eps(Float64))
        data .-= mn
        return data
    end

    # if few pixels, just use the mean
    if length(v) < 192
        return data .- mean(v)
    end

    # build histogram manually (we could use the clamp function, but this might lead to InexacrErrors as the values can be larger than typemax(Int64))
    if mx - mn < binwidth
        nbins = 16
    elseif mx - mn > binwidth * 512
        nbins = 512
    else
        nbins = ceil(Int, (mx - mn) / binwidth)
    end

    counts, edges, _ = hist(v, nbins)

    bin_centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    max_idx = findmax(counts)[2]
    base_level = bin_centers[max_idx]

    data .-= base_level
    return data
end