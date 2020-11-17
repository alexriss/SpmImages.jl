module SpmImages

using Dates
using Plots

export load_image, get_channel, plot_channel, plot_data

@enum ScanDirection up down
@enum Direction bwd fwd

mutable struct SpmImage
    filename::String
    header::Dict
    # channels::Vector{Dict}
    data::Array{Float32}
    channel_names::Vector{String}
    channel_units::Vector{String}
    
    scansize::Vector{Float64}
    scansize_unit::String
    pixelsize::Vector{Int64}
    scan_direction::ScanDirection
    
    start_time::DateTime
    acquisition_time::Float64
end
SpmImage(filename::String) = SpmImage(filename, Dict(), [], [], [], [], "", [], up, Dates.now(), 0)

mutable struct Channel
    name::String
    unit::String
    direction::Direction
    data::Array{Float32}
end
Channel() = Channel("", "", fwd, [])


"""Loads SPM data from a Nanonis sxm file.
Args:
    fname (str): Filename of the data file to read.
    output_info (int): Specifies the amount of output info to print to stdout when reading the files. 0 for no output, 1 for limited output, 2 for detailed output.
    header_only (bool): If true, then only header data is read, image data is ignored

Raises:
    Error: If the file extension is not known.
"""
function load_image(fname::String; output_info::Int64=1, header_only::Bool=false)
    if !isfile(fname)
        error("Cannot find file $fname")
        return nothing
    end
    
    image = SpmImage(fname)
    ext = rsplit(fname, "."; limit=2)[2]
    if ext == "sxm"
        _load_image_nanonis!(image, output_info, header_only)
    else
        throw(ErrorException("Error: Unknown file type: $ext"))
    end
    image
end


"""Simplifies a string (i.e. replaces space for "_", and makes it lowercase

Args:
    str (str): Input string.

Returns:
    str: Simplified output string.
"""
function _string_simplify(str::String)
    return lowercase(replace(str,' ' => '_'))
end


"""Un-simplifies (beautifies) a string (i.e. replaces "_" for space

Args:
    str (str): Input string.

Returns:
    str: Un-simplified output string.
"""
function _string_unsimplify(str::String)
    return replace(str,'_' => ' ')
end


"""Gets channel information from image headers

Args:
    image (SpmImage): SpmImage object
Returns:
    names, units: a tuple of arrays of strings specifying the channel names and their respective units; or empty tuple if information cant be extracted
"""
function _get_channel_names_units(image::SpmImage)
    lines = split(image.header["data_info"], "\n")
    popfirst!(lines) # first row are the headers: Channel Name Unit Direction Calibration Offset
    names = []
    units = []
    for line in lines
        entries = split(line)
        if length(entries) > 1
            push!(names, entries[2])
            push!(units, entries[3])
            if entries[4] != "both"
                error("Error: Only one direction recorded. This is not implemented yet. ($entries)")
                return ()
            end
        end
    end
    return (names, units)
end


"""Loads header data from a Nanonis .sxm file

Args:
    image (SpmImage): SpmImage struct
    output_info (int): Specifies the amount of output info to print to stdout when reading the files. 0 for no output, 1 for limited output, 2 for detailed output.
    header_only (bool): If true, then only header data is read, image data is ignored
"""
function _load_image_nanonis!(image::SpmImage, output_info::Int64=1, header_only::Bool=false)
    if output_info > 0
        println("Reading header of $(image.filename)")
    end
    
    header_ended = false
    caption = r"^:.*:"
    key = ""
    contents = ""
    open(image.filename) do f
        # read header data
        for line in eachline(f)
            if line == ":SCANIT_END:"  # end of header
                header_ended = true
                image.header[key] = contents
            else
                if occursin(caption, line)
                    if key != ""
                        image.header[key] = strip(contents)
                    end
                    key = _string_simplify(line[2:end-1])  # set new name
                    contents = ""  # reset contents
                else  # if not caption, it is content
                    if contents != ""
                        contents *= "\n"
                    end
                    contents *= line
                end
            end
            
            if header_ended
                break
            end
        end
        
        # parse some of the extracted data
        image.scansize = map(x -> parse(Float64, x) * 1e9, split(image.header["scan_range"]))  # 1e9 is to convert to nm
        image.scansize_unit = "nm"
        image.pixelsize = map(x -> parse(Int64, x), split(image.header["scan_pixels"]))
        image.scan_direction = image.header["scan_dir"] == "up" ? up : down

        image.start_time = DateTime(image.header["rec_date"] * " " * image.header["rec_time"], dateformat"d.m.Y H:M:S")
        image.acquisition_time = parse(Float64, image.header["acq_time"])
        
        # read body
        if !header_only
            if output_info > 0
                println("Reading body of $(image.filename)")
            end
            r = _get_channel_names_units(image)
            if length(r) != 2
                return
            end
            image.channel_names, image.channel_units = r
            
            num_channels = length(image.channel_names) * 2    # the "*2" is there because of forward and backward channels
            x_pixels, y_pixels = image.pixelsize
            data = Array{Float32}(undef, x_pixels, y_pixels, num_channels)
            skip(f, 4)
            read!(f, data)
            image.data = ntoh.(data)  # big-endian to host endian
        end
    end
end


"""Gets the channel dictionary for channel_name
Args:
    image (SpmImage): SpmImage object
    channel_name (str): string specifying the channel name to return, backward channels are generated by a suffix "_bwd"
    origin (str): "upper" or "lower". heatmap function usually uses "lower" origin and imageview functions "upper" origin
Returns:
    Channel: struct containing the channel name, unit and 2D data corresponding to channel_name.
"""
function get_channel(image::SpmImage, channel_name::String; origin::String="lower")
    backward = false
    if channel_name in image.channel_names
        i_channel = findfirst(x -> x == channel_name, image.channel_names)
    elseif endswith(channel_name, "_bwd") && channel_name[1:end-4] in image.channel_names
        i_channel = findfirst(x -> x == channel_name, image.channel_names)
        backward = true
    elseif endswith(channel_name, "_fwd") && channel_name[1:end-4] in image.channel_names
        i_channel = findfirst(x -> x == channel_name, image.channel_names)
    else  # try lower case, as well as forward and backward suffixes
        channel_name_ = lowercase(channel_name)
        channel_names = lowercase.(image.channel_names)
        if occursin(r"[^a-zA-Z0-9](bwd|b|back|backward|backwards)$", channel_name_)
            backward = true
        end
        channel_name_ = replace(channel_name_, r"[^a-zA-Z0-9](bwd|b|back|backward|backwards|f|fwd|forward)$" => "")
        if !(channel_name_ in channel_names)
            channel_name_ = replace(channel_name_, r"[^a-zA-Z0-9]" => "")
            channel_names = replace.(channel_names, r"[^a-zA-Z0-9]" => "")
        end
        if !(channel_name_ in channel_names)
            println("Error: Channel $channel_name does not exist in $(image.filename).")
            println("Available channels are: " * join(image.channel_names, ", ") * ".")
            error("Channel $channel_name not found.")
            return Channel()
        end
        i_channel = findfirst(x -> x == channel_name_, channel_names)
    end

    if backward
        data = transpose(reverse(image.data[:, :, i_channel * 2], dims=1))  # *2 for forward and backward
        direction = bwd
    else
        data = transpose(image.data[:, :, i_channel * 2 - 1])  # -1 because forward scan 
        direction = fwd
    end
    
    if origin == "upper" && image.scan_direction == up
        data = reverse(data, dims=1)
    elseif origin == "lower" && image.scan_direction == down
        data = reverse(data, dims=1)
    end
        
    return Channel(image.channel_names[i_channel], image.channel_units[i_channel], direction, data)
end


"""Plots the channel image
Args:
    image (SpmImage): SpmImage object
    channel_name (str): string specifying the channel name to return, backward channels are generated by a suffix "_bwd"
    pixel_units (bool): specifies whether to use pixel units (otherwise physical units are used)
    args: extra keyword arguments that will be passed on to plot_data (and there to heatmap)
Returns:
    plot of channel
"""
function plot_channel(image::SpmImage, channel_name::String; pixel_units::Bool=false, args...)
    channel = get_channel(image, channel_name)
    title = _string_unsimplify(channel_name) * " [$(channel.unit)]"

    if pixel_units
        x_label = "px"
        y_label = "px"
        return plot_data(channel.data, title=title, x_label=x_label, y_label=y_label; args...)
    else
        x_label = image.scansize_unit
        y_label = image.scansize_unit
        return plot_data(channel.data, title=title, x_label=x_label, y_label=y_label, scansize=image.scansize; args...)
    end
end


"""Plots 2D data
Args:
    data (array): 2d data rray to be plotted
    title (str): string specifying the image title
    x_label (str): x label
    y_label (str): y label
    scansize (array of float): scansize in physical units
    args: extra keyword arguments to be passed to the heatmap function
Returns:
    plot of the data
"""
function plot_data(data::Array{<:Number,2}; title::String="", x_label::String="", y_label::String="", scansize::Vector{<:Number}=Vector{Float64}(undef, 0), args...)
    if length(scansize) == 2  # physical units
        xs = range(0, scansize[1], length=size(data)[1])
        ys = range(0, scansize[2], length=size(data)[2])
        xlim = (0, scansize[1])
        ylim = (0, scansize[2])
        p = heatmap(xs, ys, data, aspect_ratio=1, color=:grays; args...)
    else  # pixel units
        p = heatmap(data, aspect_ratio=1, color=:grays; args...)
        xlim = (0, size(data)[1])
        ylim = (0, size(data)[2])
    end
    title!(p, title)
    xlabel!(p, x_label)
    ylabel!(p, y_label)
    xlims!(p, xlim)
    ylims!(p, ylim)
    return p
end


end