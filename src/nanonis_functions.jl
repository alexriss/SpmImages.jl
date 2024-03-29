"""
    get_channel_names_units_nanonis(image::SpmImage)

Gets a list of all channel names and units in `image`.
"""
function get_channel_names_units_nanonis(image::SpmImage)
    #returns Tuple{Vector{String},Vector{String}}
    
    lines = split(image.header["Data info"], "\n")
    popfirst!(lines) # first row are the headers: SpmImageChannel Name Unit Direction Calibration Offset
    names = Vector{String}(undef, 0)
    units = Vector{String}(undef, 0)
    for line in lines
        entries = split(line)
        if length(entries) > 1
            push!(names, replace(entries[2], "_" => " "))  # this is a quick fix (maybe one should get the actual names from "Scan>channels", but they are not in the same order
            push!(units, entries[3])
            if entries[4] != "both"
                error("Error: Only one direction recorded. This is not implemented yet. ($entries)")
                return ()
            end
        end
    end
    return (names, units)
end


"""
    load_image_nanonis(fname::String, output_info::Int=1, header_only::Bool=false)

Loads data from `fname`, specifying the file name of a Nanonis .sxm file, and returns an `SpmImage` object.

`output_info` specifies the amount of output info to print to stdout when reading the files. 0 for no output, 1 for limited output, 2 for detailed output.
If `header_only` is `true` then only header data is read, image data is ignored.
"""
function load_image_nanonis(fname::String, output_info::Int=1, header_only::Bool=false)
    image = SpmImage(fname, sxm)

    output_info > 0 && println("Reading header of $(image.filename)")
    
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
                    key = string_prettify(line[2:end-1])  # set new name
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
        image.scansize = map(x -> parse(Float64, x) * 1e9, split(image.header["Scan range"]))  # 1e9 is to convert to nm
        image.scansize_unit = "nm"
        image.center = map(x -> parse(Float64, x) * 1e9, split(image.header["Scan offset"]))  # 1e9 is to convert to nm
        image.angle = parse(Float64, image.header["Scan angle"])
        image.pixelsize = map(x -> parse(Int, x), split(image.header["Scan pixels"]))
        image.scan_direction = image.header["Scan dir"] == "up" ? up : down

        image.bias = parse(Float64, image.header["Bias"])
        
        z_controller_data = split(image.header["Z-controller"], "\t")
        image.z_feedback = z_controller_data[8] == "1" ? true : false
        image.z_feedback_setpoint = parse(Float64, split(z_controller_data[9])[1])
        image.z_feedback_setpoint_unit = split(z_controller_data[9])[2]
        if haskey(image.header, "Z-Controller>Z (m)")
            image.z = parse(Float64, image.header["Z-Controller>Z (m)"])
        end

        image.start_time = DateTime(image.header["Rec date"] * " " * image.header["Rec time"], dateformat"d.m.Y H:M:S")
        image.acquisition_time = parse(Float64, image.header["Acq time"])

        r = get_channel_names_units_nanonis(image)
        if length(r) != 2
            return
        end
        image.channel_names, image.channel_units = r

        # read body
        if !header_only
            output_info > 0 && println("Reading body of $(image.filename)")

            num_channels = length(image.channel_names) * 2    # the "*2" is there because of forward and backward channels
            x_pixels, y_pixels = image.pixelsize
            data = Array{Float32}(undef, x_pixels, y_pixels, num_channels)
            image.channel_indices_fwd = collect(1:2:num_channels)
            image.channel_indices_bwd = collect(2:2:num_channels)
            skip(f, 4)
            read!(f, data)
            image.data = ntoh.(data)  # big-endian to host endian
        end
    end
    return image
end
