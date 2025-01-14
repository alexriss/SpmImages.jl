module SpmImages

using CoordinateTransformations
using DataStructures:OrderedDict
using Dates
using ImageFiltering
using ImageTransformations
using NetCDF
using Printf
using StaticArrays
using Statistics
using StructIO
using TOML

export SpmImage, load_image, get_channel, plot_channel, plot_data, correct_background, line_profile
export calc_drift_xy, correct_drift!, correct_drift
export pixels_to_nm, nm_to_pixels
export Background
export no_correction, subtract_minimum, plane_linear_fit, line_average, vline_average, line_linear_fit, vline_linear_fit

const VERSION = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "../Project.toml"))["version"]) 


@enum FileType sxm nc ibw bru
@enum ScanDirection up down
@enum Direction bwd fwd
@enum Background no_correction subtract_minimum plane_linear_fit line_average vline_average line_linear_fit vline_linear_fit line_linear_fit_legacy vline_linear_fit_legacy
@enum DriftCorrection drift_none drift_translation drift_full
const DriftCorrection_display = Dict(
    drift_none => "none",
    drift_translation => "translation",
    drift_full => "full"
)

mutable struct SpmImage
    filename::Union{String,Vector{String}}
    type::FileType
    header::AbstractDict
    data::Array{Float32}

    channel_names::Vector{String}
    channel_units::Vector{String}
    channel_indices_fwd::Vector{Int}  # indices of fwd channels in `data`. 0-index means that this direction is missing
    channel_indices_bwd::Vector{Int}  # indices of fwd channels in `data`. 0-index means that this direction is missing
    
    scansize::Vector{Float64}
    scansize_unit::String
    center::Vector{Float64}
    angle::Float64
    pixelsize::Vector{Int}
    scan_direction::ScanDirection

    bias::Float64
    z_feedback::Bool
    z_feedback_setpoint::Float64
    z_feedback_setpoint_unit::String
    z::Float64
    
    start_time::DateTime
    acquisition_time::Float64

    drift_correction::DriftCorrection
    drift::Vector{Float64}   # drift in units of [scansize] / [acquisition_time]
end
SpmImage(filename::Union{String,Vector{String}}, type::FileType) = SpmImage(filename, type, OrderedDict(), Float32[],
    String[], String[], String[], String[],
    Float64[], "", Float64[], 0., Float64[], up,
    0., true, NaN, "", 0.,
    Date(-2), 0.,
    drift_none, Float64[])

mutable struct SpmImageChannel
    name::String
    unit::String
    direction::Direction
    data::Array{Float32}
end
SpmImageChannel() = SpmImageChannel("", "", fwd, [])


include("nanonis_functions.jl")
include("netCDF_functions.jl")
include("ibw_functions.jl")
include("bru_functions.jl")
include("plot_functions.jl")
include("drift_functions.jl")


function Base.show(io::IO, i::SpmImage)
    if get(io, :compact, false)
        fname = basename.(i.filename)
        fname_str = typeof(fname) == String ? fname : fname * "... (" * string(length(i.filename)) * " files)"
        print(io, "SpmImage(\"", fname_str, "\")")
    else
        fname = basename.(i.filename)
        fname_str = typeof(fname) == String ? fname : fname[1] * "... (" * string(length(i.filename)) * " files)"
        b = @sprintf "%0.2g" i.bias
        scansize = join([@sprintf("%0.2g", x) for x in i.scansize], " x ")
        pixelsize = join([string(x) for x in i.pixelsize], " x ")
        if scansize == ""
            scansize = "-"
        end
        if pixelsize == ""
            pixelsize = "-"
        end
        dc = DriftCorrection_display[i.drift_correction]
        print(io, "SpmImage(\"", fname_str, "\", ",
        "bias: ", b, " V, ",
        length(i.channel_names), " channels, ",
        scansize, " nm, ", pixelsize, " pixels, ",
        "drift correction: ", dc, ")")
    end
end


"""
    load_image(fname::Union{String,Vector{String}}; output_info::Int=0, header_only::Bool=false)::SpmImage

Loads an SPM image from a file, and returns an `SpmImage` object. The file extension is used to determine the file type.
Nanonis .sxm files are supported, as well as GXSM netCDF files. For Nanonis .sxm files, only one file can be read - if many files are specified, only the first one is read.
For netCDF files, the `fname` argument can be a vector of strings, specifying all the file names of the netCDF files to load.
Most of the header information is then read from the first file.

`output_info` specifies the amount of output info to print to stdout when reading the files. 0 for no output, 1 for limited output, 2 for detailed output.
If `header_only` is `true` then only header data is read, image data is ignored.
"""
function load_image(fname::Union{String,Vector{String}}; output_info::Int=0, header_only::Bool=false)
    #returns SpmImage

    if isa(fname, String)
        fname = [fname]
    end

    ext = splitext(fname[1])[2]
    if ext == ".sxm"
        image = load_image_nanonis(fname[1], output_info, header_only)
    elseif ext == ".nc"
        image = load_image_netCDF(fname, output_info, header_only)
    elseif ext == ".ibw"
        image = load_image_ibw(fname[1], output_info, header_only)
    elseif ext == ".spm"
        image = load_image_bru_spm(fname[1], output_info, header_only)
    else
        throw(ErrorException("Error: Unknown file type: $ext"))
    end
    return image
end


"""
    get_i_channel(image::SpmImage, channel_name::String)

Returns the channel index and scan direction corresponding to `channel_name`. This is used for indexing `image.data.`
"""
function get_channel_index(image::SpmImage, channel_name::String)
    backward = false
    if channel_name in image.channel_names
        i_channel = findfirst(x -> x == channel_name, image.channel_names)
    elseif endswith(channel_name, " bwd") && channel_name[1:end-4] in image.channel_names
        i_channel = findfirst(x -> x == channel_name[1:end-4], image.channel_names)
        backward = true
    elseif endswith(channel_name, " fwd") && channel_name[1:end-4] in image.channel_names
        i_channel = findfirst(x -> x == channel_name[1:end-4], image.channel_names)
    else  # try lower case, as well as forward and backward suffixes
        channel_name_ = string_simplify(channel_name)
        channel_names = string_simplify.(image.channel_names)
        reg = r"[^a-zA-Z0-9](bwd|b|back|backward|backwards|retrace|\[bwd\]|\[b\]|\[back\]|\[backward\]|\[backwards\])_*$"
        if occursin(reg, channel_name_)
            backward = true
        end
        channel_name_ = replace(channel_name_, reg => "")
        if !(channel_name_ in channel_names)
            channel_name_ = replace(channel_name_, r"[^a-zA-Z0-9]" => "")
            channel_names = replace.(channel_names, r"[^a-zA-Z0-9]" => "")
        end
        if !(channel_name_ in channel_names)
            println("Error: Channel $channel_name does not exist in $(image.filename).")
            println("Available channels are: " * join(image.channel_names, ", ") * ".")
            error("Channel $channel_name not found.")
            return SpmImageChannel()
        end
        i_channel = findfirst(x -> x == channel_name_, channel_names)
    end

    direction = backward ? bwd : fwd
    return i_channel, direction
end


"""
    get_channel(image::SpmImage, channel_name::String; origin::String="lower")

Gets the channel object for `channel_name` in `image`.
 `origin` specifies the origin of the image, either "upper" or "lower".
"""
function get_channel(image::SpmImage, channel_name::String; origin::String="lower")
    i_channel, direction = get_channel_index(image, channel_name)

    # if image.type == sxm
    #     return get_channel_nanonis(image, i_channel, direction; origin=origin)
    # elseif image.type == nc
    #     return get_channel_netCDF(image, i_channel, direction; origin=origin)
    # elseif image.type == ibw
    #     return get_channel_ibw(image, i_channel, direction; origin=origin)
    # end

    empty = false
    if direction == bwd
        i = image.channel_indices_bwd[i_channel]
        if i == 0 
            data = fill(NaN32, reverse(image.pixelsize)...)
            empty = true
        else
            if image.type == ibw || image.type == bru
                @views data = transpose(image.data[:, :, i])
            else
                @views data = transpose(reverse(image.data[:, :, i], dims=1))
            end
            
        end
    else
        i = image.channel_indices_fwd[i_channel]
        if i == 0 
            data = fill(NaN32, reverse(image.pixelsize)...)
            empty = true
        else
            @views data = transpose(image.data[:, :, i])
        end
    end

    # drift correction
    if !empty && image.drift_correction === drift_full
        data = drift_corr_data(image, data)
    end
    
    if image.type == sxm || image.type == ibw
        if !empty && origin == "upper" && image.scan_direction == up
            data = reverse(data, dims=1)
        elseif !empty && origin == "lower" && image.scan_direction == down
            data = reverse(data, dims=1)
        end
    elseif image.type == nc
        if origin == "lower" 
            data = reverse(data, dims=1)
        end
    elseif image.type == bru
        if origin == "upper" 
            data = reverse(data, dims=1)
        end
    end
        
    return SpmImageChannel(image.channel_names[i_channel], image.channel_units[i_channel], direction, data)
end


"""Simplifies a string (i.e. replaces space for "_", and makes it lowercase

Args:
    str (str): Input string.
Returns:
    str: Simplified output string.
"""
function string_simplify(str::String)
    return lowercase(replace(str,' ' => '_'))
end


"""Prettifies (~un-simplifies) a string (i.e. replaces "_" for space

Args:
    str (str): Input string.
Returns:
    str: Prettified output string.
"""
function string_prettify(str::String)
    if uppercase(str) == str
        str = uppercasefirst(lowercase(str))
    end
    return replace(str,'_' => ' ')
end


"""converts nm units to pixel units

Args:
    image (SpmImage): SpmImage object.
    point (Vector of Number): x and y coordinates in scan units (nm)
    origin (str): specifying the pixel-data origin (same as for get_channel)
Returns:
    point (Vector of Number): x and y (column and row) coordinates in pixel units
"""
function nm_to_pixels(image::SpmImage, point::Vector{<:Number}, origin::String="lower")
    #returns Vector{<:Number}

    pixelsize = image.pixelsize .- 1
    if origin == "upper"  # pixel-origin is top left, scan-units origin is bottom left
        return 1 .+ [point[1] * pixelsize[1] / image.scansize[1], pixelsize[2] - point[2] * pixelsize[2] / image.scansize[2]]  # 1.+ because of 1-based indexing
    else
        return 1 .+ [point[1] * pixelsize[1] / image.scansize[1], point[2] * pixelsize[2] / image.scansize[2]]  # 1.+ because of 1-based indexing
    end
end


"""converts pixel units to nm units

Args:
    image (SpmImage): SpmImage object.
    point (Vector of Number): x and y (column and row) coordinates in pixel units
    origin (str): specifying the pixel-data origin (same as for get_channel)
Returns:
    point (Vector of Number): x and y coordinates in nm units
"""
function pixels_to_nm(image::SpmImage, point::Vector{<:Number}, origin::String="lower")
    #returns Vector{<:Number}

    pixelsize = image.pixelsize .- 1
    x, y = point[1] - 1, point[2] - 1  # -1 because of 1-based indexing
    if origin == "upper"  # pixel-origin is top left, scan-units origin is bottom left
        return [x / pixelsize[1] * image.scansize[1], image.scansize[2] - y / pixelsize[2] * image.scansize[2]]
    else
        return [x / pixelsize[1] * image.scansize[1], y / pixelsize[2] * image.scansize[2]]
    end
end


"""Returns value at pixel coordinates. Float coordinates are rounded to the nearest integer.

Args:
    data: 2d array
    point (Vector of Number): x and y (column and row) coordinates in pixel units
Returns:
    value (Number): value at the coordinates (or missing)
"""
function get_value(data::Array{<:Number,2}, point::Vector{<:Number})
    #returns Union{Number,Missing}

    col, row = round.(Int, point)
    if (row < 1 || col < 1)
        return missing
    elseif (row > size(data)[1] || col > size(data)[2])
        return missing
    else
        return data[row,col]
    end
end


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


"""Calculates a line profile
Args:
    image (SpmImage): SpmImage object
    channel_name (str): string specifying the channel name to return, backward channels are generated by a suffix " bwd".
    background (Background): type of background correction.
    start_point (1D array of Number): coordinates of the start point (in scan-units)
    end_point (1D array of Number): coordinates of the start point (in scan-units)
    width (Number):: width of the line (in scan-units)
Returns:
    A tuple of three arrays (x/y coordinates, lengths, data values) and two numbers (start point value and end point value).
    The line profile values are subject to the line width, whereas the start and end point values are values at the exact point.
"""
function line_profile(image::SpmImage, channel_name::String, start_point::Vector{<:Number}, end_point::Vector{<:Number}, width::Number=0; background::Background=no_correction)
    #returns Tuple{Vector{Vector{<:Number}}, Vector{<:Number}, Vector{Union{<:Number,Missing}}, Union{<:Number,Missing}, Union{<:Number,Missing}}

    # check dimensions of start and end points
    length(start_point) == 2 || throw(ArgumentError("Invalid length of start_point ($start_point): should be of length 2"))
    length(end_point) == 2 || throw(ArgumentError("Invalid length of end_point ($end_point): should be of length 2"))

    origin = "lower"
    channel = get_channel(image, channel_name, origin=origin)
    if background != no_correction
        data = correct_background(channel.data, background)
    else
        data = channel.data
    end

    return line_profile(image, data, start_point, end_point, width, origin=origin)
end


"""Calculates a line profile for 2d data.
Args:
    image (SpmImage): SpmImage object
    data (2D array): channel data.
    start_point (1D array of Number): coordinates of the start point (in scan-units)
    end_point (1D array of Number): coordinates of the start point (in scan-units)
    width (Number):: width of the line (in scan-units)
    origin (str): origin of the image ("lower" or "upper")
Returns:
    A tuple of three arrays (x/y coordinates, lengths, data values) and two numbers (start point value and end point value).
    The line profile values are subject to the line width, whereas the start and end point values are values at the exact point.
"""
function line_profile(image::SpmImage, data::Array{<:Number,2},
    start_point::Vector{<:Number}, end_point::Vector{<:Number}, width::Number=0;
    origin::String="lower")
    #returns Tuple{Vector{Vector{<:Number}}, Vector{<:Number}, Vector{Union{<:Number,Missing}}, Union{<:Number,Missing}, Union{<:Number,Missing}}

    # x and y width (nm units)
    d_x, d_y = end_point - start_point
    theta = atan(d_y, d_x)  # arctan2
    width_x = width * sin(-theta) / 2
    width_y = width * cos(theta) / 2

    # now switch to pixel units
    start_point_px = nm_to_pixels(image, start_point, origin)
    end_point_px = nm_to_pixels(image, end_point, origin)
    width_col = width_x * image.pixelsize[1] / image.scansize[1]
    width_row = width_y * image.pixelsize[2] / image.scansize[2]

    # similar to skimage.measure_profile_line
    @views src_row, src_col = start_point_px[2], start_point_px[1]  # row is y, col is x
    @views dst_row, dst_col = end_point_px[2], end_point_px[1]  # row is y, col is x
    d_col, d_row = end_point_px - start_point_px
    length_px = ceil(Int, hypot(d_row, d_col) + 1)  # +1 to include end-point
    line_col = collect(range(src_col, dst_col, length=length_px))
    line_row = collect(range(src_row, dst_row, length=length_px))
    width_length_px = ceil(Int, hypot(width_row, width_col) + 1)  # +1 to include end-point

    t_data = Union{typeof(first(data)),Missing}
    distances = Vector{Float64}(undef, length_px)  # vector of distances
    coords = fill([0.0, 0.0], length_px)  # vector of points
    values = Vector{t_data}(undef, length_px)  # vector of values

    for i in 1:length_px
        row_i = line_row[i]
        col_i = line_col[i]
        coords[i] = pixels_to_nm(image, [col_i, row_i], origin)
        distances[i] = (i == 1) ? 0 : hypot((coords[i] - coords[1])...)
        perp_rows = collect(range(row_i - width_row, row_i + width_row, length=width_length_px))
        perp_cols = collect(range(col_i - width_col, col_i + width_col, length=width_length_px))
        values_perp = Vector{t_data}(undef, width_length_px)
        for j in 1:width_length_px
            values_perp[j] = get_value(data, [perp_cols[j],perp_rows[j]])
        end
        values_perp_no_missing = collect(skipmissing(values_perp))
        if (length(values_perp_no_missing) == 0)
            values[i] = missing
        else
            values[i] = mean(values_perp_no_missing)
        end
    end

    start_point_value = get_value(data, start_point_px)
    end_point_value = get_value(data, end_point_px)

    return coords, distances, values, start_point_value, end_point_value
end


end