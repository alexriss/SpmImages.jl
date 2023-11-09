
"""
    crosscorr_norm(img::Matrix{Float32}, template::Matrix{Float32})::Matrix{Float32}

Calculates a normalized cross-correlation between two images.
See: https://discourse.julialang.org/t/cross-correlation-function-in-julia-vs-python/7573/7
"""
function crosscorr_norm(img::Matrix{Float32}, template::Matrix{Float32})::AbstractArray{Float32}
    # z = imfilter(Float64, centered(img), centered(template))
    # q = sqrt.(imfilter(Float64, centered(img.^2), centered(ones(size(template)))))
    z = imfilter(Float32, centered(img), centered(template))
    q = sqrt.(imfilter(Float32, centered(img.^2), centered(ones(size(template)))))
    m = sqrt(sum(template.^2))
    corr = z ./ (m * max.(eps(), q))
    return corr
end


"""
    calc_drift(filename1::AbstractString, filename2::AbstractString, channel="Z")::Vector{AbstractFloat}

Calculates xy drift between `image1` and `image2` via cross-correlations of `channel` of the two images.
The images do not need to be equal in size, but rotation is not accounted for explicitely.
Returns the x and y drift in `nm/s`.

*Experimental and not well tested yet.*
"""
function calc_drift_xy(filename1::AbstractString, filename2::AbstractString, channel="Z")::Vector{Float64}
    image1 = load_image(filename1)
    image2 = load_image(filename2)
    return calc_drift_xy(image1, image2, channel)
end


"""
    calc_drift(image1::SpmImage, image2::SpmImage, channel="Z")::Vector{AbstractFloat}

Calculates xy drift between `image1` and `image2` via cross-correlations of `channel` of the two images.
The images do not need to be equal in size, but the rotation of both should be the same.
Returns the x and y drift in `nm/s` with respect to the angle of the images.

*Experimental and not well tested yet.*
"""
function calc_drift_xy(image1::SpmImage, image2::SpmImage, channel="Z")::Vector{Float64}
    if image1.angle != image2.angle
        throw(ArgumentError("Both images should have the same rotation. Rotations found: $(image1.angle) and $(image2.angle) degrees."))
    end

    ch1 = get_channel(image1, channel).data
    ch2 = get_channel(image2, channel).data

    density1 = image1.pixelsize ./ image1.scansize
    density2 = image2.pixelsize ./ image2.scansize

    # adapt densities, use lower density
    image_base = image1
    if any(density2 .> density1)
        newsize = round.(Int, Tuple(density1 .* image2.scansize))
        ch2 = imresize(ch2, reverse(newsize))  # imresize wants row/col order
    elseif any(density1 .> density2)
        newsize = round.(Int, Tuple(density2 .* image2.scansize))
        ch1 = imresize(ch1, reverse(newsize))  # imresize wants row/col order
        image_base = image2
    end

    # remove NaN values - they cause problems in the crosscorrelation
    ch1[isnan.(ch1)] .= 0f0
    ch2[isnan.(ch2)] .= 0f0
    corr = SpmImages.crosscorr_norm(ch1, ch2)
    # calculate cross-correlation
    # corr = imfilter(corr, Kernel.gaussian(5))

    # find maximum
    drift_px = -collect(reverse(Tuple(findmax(corr)[2])))
    # @show drift_px
    diff_secs = datetime2unix(image2.start_time) - datetime2unix(image1.start_time) +
         (image2.acquisition_time - image1.acquisition_time) / 2
    drift_nm = pixels_to_nm(image_base, drift_px) - image1.center + image2.center
    # @show diff_secs, drift_nm, drift_nm ./ diff_secs

    return drift_nm ./ diff_secs
end


"""
    correct_drift!(image::SpmImage, drift::Vector{Float64};
        ref_time::Union{DateTime,SpmImage,Nothing}=nothing, skew::Bool=true)::Nothing

Applies a drift correction to `image`. It re-calculates the image center and z-value
assuming that the `drift` started at `ref_time-time`. For `ref_time` a `DateTime` can be
given, or an `SpmImage`. If `nothing` is given, then no change of the origin is made.
`drift` should be a 2-element or 3-element vector for x,y or x,y,z drift, respecitvely.
If `skew` is `true` (default), then the image will be
transformed and the image dimensions updated accodingly.

*Experimental and not well tested yet.*
"""
function correct_drift!(image::SpmImage, drift::Vector{Float64};
    ref_time::Union{DateTime,SpmImage,Nothing}=nothing, skew::Bool=true)::Nothing

    @assert length(drift) >= 2

    image.drift = drift

    if isnothing(ref_time)
        diff_secs = 0.
    elseif isa(ref_time, SpmImage)
        start_time_secs = datetime2unix(ref_time.start_time) + ref_time.acquisition_time / 2
        diff_secs = datetime2unix(image.start_time) + image.acquisition_time / 2 - start_time_secs
    else
        start_time_secs = datetime2unix(start_time)
        diff_secs = datetime2unix(image.start_time) + image.acquisition_time / 2 - start_time_secs
    end

    driftxy = @view drift[1:2]
    # @show image.center
    image.center -= driftxy .* diff_secs
    # @show image.center, driftxy .* diff_secs
    if length(drift) > 2
        image.z -= drift[3] .* diff_secs
    end

    if skew
        density = image.pixelsize ./ image.scansize
        # test transformation on dummy matrix
        dummy = zeros(Float32, reverse(image.pixelsize)...)
        dummyw = drift_corr_data(image, dummy)
        image.pixelsize = collect(reverse(size(dummyw)))
        image.scansize = image.pixelsize ./ density
        image.drift_correction = drift_full
    else
        image.drift_correction = drift_translation
    end

    return nothing
end


"""
    correct_drift(image::SpmImage, drift::Vector{Float64};
    ref_time::Union{DateTime,SpmImage,Nothing}=nothing, skew::Bool=true)::SpmImage

Applies a drift correction to a copy of `image`. It re-calculates the image center and z-value
assuming that the `drift` started at `start-time`. For `ref_time` a `DateTime` can be
given, or an `SpmImage`. If `nothing` is given, then the drift is assumed to start at
the ref_time of the current image.
`drift` should be a 2-element or 3-element vector for x,y or x,y,z drift, respecitvely.
If `skew` is `true` (default), then the image will be
transformed and the image dimensions updated accodingly.

*Experimental and not well tested yet.*
"""
function correct_drift(image::SpmImage, drift::Vector{Float64};
    ref_time::Union{DateTime,SpmImage,Nothing}=nothing, skew::Bool=true)::SpmImage
    
    image_copy = deepcopy(image)
    correct_drift!(image_copy, drift, ref_time=ref_time, skew=skew)
    return image_copy
end


"""
    linear_map(image::SpmImage, drift::Union{Nothing,Vector{Float64}}=nothing)::LinearMap{Matrix{Float64}}

Calculates the linear map that corrects the xy `drift` for an `image`.
"""
function linear_map(image::SpmImage, drift::Union{Nothing,Vector{Float64}}=nothing)::LinearMap{Matrix{Float64}}
    if isnothing(drift)
        drift = image.drift
    end
    drift = -drift[1:2] ./ image.scansize  # xy drift relative to imagesize
    cols, rows = image.pixelsize

    # x-position at end of first row
    x_resize = 1 + drift[1] * image.acquisition_time / rows
    # y position at end of first row
    y_shear = drift[2] * image.acquisition_time / rows *
        rows / cols  # in multiples of cols
    # y position at beginning of last row
    y_resize = 1 + drift[2] * image.acquisition_time / rows * (rows - 1)
    # x position at beginning of last row
    x_shear = drift[1] * image.acquisition_time /
        rows * (rows - 1) *  # x posistion in the last row
        cols / rows  # in multiples of rows

    return LinearMap([y_resize y_shear; x_shear x_resize])  # image is in row/col, i.e. y/x
end


"""
    drift_corr_data(data::Matrix{Float32}, transform::LinearMap{Matrix{Float64}}, fill::Float32=NaN32)::Matrix{Float32}

Skewes `data` according to `transform`. Raw data should be given, i.e. before accounting for
the scan direction.
"""
function drift_corr_data(data::Matrix{Float32}, transform::LinearMap{Matrix{Float64}}, fill::Float32=NaN32)::Matrix{Float32}
    return parent(warp(data, transform, fill=fill))
end


"""
    function drift_corr_data(image::SpmImage, data::AbstractArray{Float32}, fill::Float32=NaN32)::Matrix{Float32}

Skewes `data` according to `image.drift`. Raw data should be given, i.e. before accounting for
the scan direction.
"""
function drift_corr_data(image::SpmImage, data::AbstractArray{Float32}, fill::Float32=NaN32)::Matrix{Float32}
    # sometimes the dimensions are off by 1, so we use the indices here to force the dimensions
    indices = (1:image.pixelsize[2], 1:image.pixelsize[1])
    return parent(warp(data, linear_map(image), indices, fill=fill))
end


