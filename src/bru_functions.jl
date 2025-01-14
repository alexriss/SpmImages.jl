# Loads bruker icon files (.spm)
#
# resources:
#   https://gist.github.com/kogens/03cbe570449201c4431b22ad08aeb407
#   
#   Example file from:
#   https://github.com/George-R-Heath/NanoLocz


"""
    get_channel_info_spm(image::SpmImage)

Gets a list of all channel names and units in `image`.
"""
function get_channel_info_spm(image::SpmImage)
    channel_info = Dict{String,Union{String,Int,Direction,Float64}}[]
    i = 1
    while haskey(image.header, "Ciao image $i - Data offset")
        info = Dict{String,Union{String,Int,Direction,Float64}}()
        name_data = image.header["Ciao image $i - Image Data"]
        pos = findall(isequal('"'), name_data)
        if length(pos) != 2
            println("Warning: Could not parse channel name. ($name_data)")
            i += 1
            continue
        end
        info["name"] = String(name_data[pos[1]+1:pos[2]-1])

        info["unit"] = "V"
        calibration_factor = 1.
        calibration_factor2 = 1.
        calibration_offset = 0.
        calibration_factor_entry = ""
        unit_data = image.header["Ciao image $i - Z scale"]
        # example entry: V [Sens. Zsens] (0.006713867 V/LSB) 440.0000 V
        pattern = r"\[(.+)\] \(([\d\.-]+) (.+)\/LSB\)\s+([\d.-]+) "
        unit_data = replace(unit_data, "\xba" => "°")  # replace deg symbol encoding
        m = match(pattern, unit_data)
        if isnothing(m)
            println("Warning: Could not parse channel unit. ($unit_data)")
        else
            calibration_factor_entry = "@" * m.captures[1]
            calibration_factor = parse(Float64, m.captures[2])
            info["unit"] = String(m.captures[3])
        end

        offset_data = image.header["Ciao image $i - Z offset"]
        # example: V [Sens. Zsens] (0.006713867 V/LSB)       0 V
        offset_data = replace(offset_data, "\xba" => "°")  # replace deg symbol encoding
        m = match(pattern, offset_data)
        if isnothing(m) && info["unit"] != ""
            println("Warning: Could not parse channel offset. ($offset_data)")
        else
            calibration_offset = parse(Float64, m.captures[4])
        end

        if calibration_factor_entry != ""
            if haskey(image.header, calibration_factor_entry)
                calib_data = image.header[calibration_factor_entry]
                # example: Sens. Zsens: V 6.127566 nm/V
                pattern = r"([\d\.-]+) (.+)\/(.+)"
                m = match(pattern, calib_data)
                if isnothing(m)
                    pattern = r"([\d\.-]+)"  # try without the unit
                    m = match(pattern, calib_data)
                end

                if isnothing(m)
                    println("Warning: Could not parse calibration factor. ($calib_data)")
                else
                    if length(m.captures) == 1
                        calibration_factor2 = parse(Float64, m.captures[1])
                    elseif m.captures[3] == info["unit"]
                        info["unit"] = String(m.captures[2])
                        calibration_factor2 = parse(Float64, m.captures[1])
                    elseif "m" * m.captures[3] == info["unit"]   # milli in channel unit, todo: more general
                        info["unit"] = String(m.captures[2])
                        calibration_factor2 = parse(Float64, m.captures[1]) * 1e-3
                    else
                        println("Warning: Unit of calibration factor ($(m.captures[3])) does not match channel unit ($(info["unit"])). ($calib_data)")
                    end
                end
            else
                println("Warning: Could not find calibration factor '$calibration_factor_entry'.")
            end
        end

        info["calibration_factor"] = calibration_factor
        info["calibration_factor2"] = calibration_factor2
        info["calibration_offset"] = calibration_offset
                
        # just do some sanity checks
        scan_direction = image.header["Ciao image $i - Frame direction"] == "Up" ? up : down
        pixelsize = [parse(Int, image.header["Ciao image $i - Samps/line"]), parse(Int, image.header["Ciao image $i - Number of lines"])]
        info["x_pixels"] = pixelsize[1]
        info["y_pixels"] = pixelsize[2]
        scansize_parts = split(image.header["Ciao image $i - Scan Size"])
        scansize = parse.(Float64, scansize_parts[1:2])
        if length(scansize_parts) == 3 && scansize_parts[3] != "nm"  # sometimes it is nm
            # we assume otherwise it is um
            scansize .*= 1000
        end
        
        if scan_direction != image.scan_direction
            println("Warning: Direction of channel '$(info["name"])' is different from main scan direction. This is unexpected.")
        end
        # this can happen for partial images
        # if pixelsize[1] != image.pixelsize[1] || pixelsize[2] != image.pixelsize[2]
        #     println("Warning: Pixelsize of channel '$(info["name"])' ($pixelsize) is different from main pixelsize ($(image.pixelsize)). This is unexpected.")
        # end
        if !isapprox(scansize, image.scansize)
            println("Warning: Scansize of channel '$(info["name"])' ($scansize) is different from main scansize ($(image.scansize). This is unexpected.")
        end

        info["direction"] = image.header["Ciao image $i - Line Direction"] == "Trace" ? fwd : bwd

        info["data_offset"] = parse(Int, image.header["Ciao image $i - Data offset"])
        info["data_length"] = parse(Int, image.header["Ciao image $i - Data length"])
        info["bytes_per_pixel"] = parse(Int, image.header["Ciao image $i - Bytes/pixel"])

        push!(channel_info, info)
        i += 1
    end
    
    return channel_info
end


"""
    addtoheader_icon_spm!(image::SpmImage, key::Union{String,SubString{String}}, value::Union{String,SubString{String}})

Adds a key-value pair to the header of `image`. If the key already exists, a number is appended to the key.
"""
function addtoheader_icon_spm!(image::SpmImage, key::Union{String,SubString{String}}, value::Union{String,SubString{String}})
    if haskey(image.header, key)
        i = 1
        while haskey(image.header, key * " $i")
            i += 1
        end
        key *= " $i"
    end
    image.header[key] = value
end


"""
    function load_image_bru_spm(fname::String, output_info::Int=1, header_only::Bool=false)

Loads data from `fname`, specifying the file name of a Bruker .spm file, and returns an `SpmImage` object.

`output_info` specifies the amount of output info to print to stdout when reading the files. 0 for no output, 1 for limited output, 2 for detailed output.
If `header_only` is `true` then only header data is read, image data is ignored.
"""
function load_image_bru_spm(fname::String, output_info::Int=1, header_only::Bool=false)
    image = SpmImage(fname, bru)

    output_info > 0 && println("Reading header of $(image.filename)")
    
    key = ""
    contents = ""
    num_ciao_image = 0
    key_prefix = ""
    open(image.filename) do f
        # read header data
        for line in eachline(f)
            if line == "\\*File list end"  # end of header
                break
            elseif startswith(line, "\\*Ciao image list")
                # new image starts
                num_ciao_image += 1
                key_prefix = "Ciao image $num_ciao_image - "
            elseif startswith(line, "\\*")  # ignore other section headers
                key_prefix = ""  # most likely not a ciao image anymore; but anyways these seem to be at the end of the headers
                continue
            elseif startswith(line, "\\@") && line[4] == ':'
                parts = split(line, ":")
                if length(parts) == 3
                    key = strip(parts[2])
                    contents = strip(parts[3])
                    addtoheader_icon_spm!(image, key_prefix * key, contents)
                elseif output_info > 0
                    println("Warning: Ignoring line:\n  $line")
                end
            elseif startswith(line, "\\@")
                parts = split(line, ":", limit=2)
                if length(parts) == 2
                    key = strip(parts[1][2:end])
                    contents = strip(parts[2])
                    addtoheader_icon_spm!(image, key_prefix * key, contents)
                elseif output_info > 0
                    println("Warning: Ignoring line:\n  $line")
                end
            elseif startswith(line, "\\")
                parts = split(line, ":", limit=2)
                if length(parts) == 2
                    key = strip(parts[1][2:end])
                    contents = strip(parts[2])
                    addtoheader_icon_spm!(image, key_prefix * key, contents)
                elseif output_info > 0
                    println("Warning: Ignoring line:\n  $line")
                end
            end
        end
        
        # parse some of the extracted data
        if haskey(image.header, "Scan Size") && haskey(image.header, "Slow Axis Size")
            s1 = split(image.header["Scan Size"])
            s2 = split(image.header["Slow Axis Size"])
            if length(s1) == 2 && length(s2) == 2  && s2[2] == "nm"
                s1_val = parse(Float64, s1[1])
                s2_val = parse(Float64, s2[1])
                image.scansize = [s1_val, s2_val]
                image.scansize_unit = "nm"
            else
                println("Warning: Could not parse scansize (parameters: 'Scan Size' and 'Slow Axis Size').")
            end
        elseif haskey(image.header, "Scan Size")  # we assume a square image
            s1 = split(image.header["Scan Size"])
            if length(s1) == 2 && s1[2] == "nm"
                s1_val = parse(Float64, s1[1])
                image.scansize = [s1_val, s1_val]
                image.scansize_unit = "nm"
            else
                println("Warning: Could not parse scansize (parameters: 'Scan Size'.")
            end           
        end

        image.center = [0, 0]
        if haskey(image.header, "X Offset") && haskey(image.header, "Y Offset")
            s1 = split(image.header["X Offset"])
            s2 = split(image.header["Y Offset"])
            if length(s1) == 2 && length(s2) == 2 && s1[2] == "nm" && s2[2] == "nm"
                s1_val = parse(Float64, s1[1])
                s2_val = parse(Float64, s2[1])
                image.center = [s1_val, s2_val]
            else
                println("Warning: Could not parse image center (parameters: 'X Offset' and 'Y Offset').")
            end
        end
        # add stage offset
        if haskey(image.header, "Stage X") && haskey(image.header, "Stage Y")
            # these dont have units, we assume um
            s1_val = parse(Float64, image.header["Stage X"]) * 1e3
            s2_val = parse(Float64, image.header["Stage Y"]) * 1e3
            image.center[1] += s1_val
            image.center[2] += s2_val
        end
        if haskey(image.header, "Rotate Ang")
            image.angle = parse(Float64, image.header["Rotate Ang"])
        end
        if haskey(image.header, "Samps/line") && haskey(image.header, "Lines")
            image.pixelsize = [parse(Int, image.header["Samps/line"]), parse(Int, image.header["Lines"])]
        end
        
        if haskey(image.header, "Capture direction")
            image.scan_direction = image.header["Capture direction"] == "Up" ? up : down
        elseif haskey(image.header, "Ciao image $i - Frame direction")
            image.scan_direction = image.header["Ciao image $i - Frame direction"] == "Up" ? up : down
        else
            println("Warning: Cant read scan direction.")
        end
        image.bias = 0.0  # todo
        image.z_feedback = true  # todo
        image.z_feedback_setpoint = 0.0  # todo
        image.z_feedback_setpoint_unit = ""  # todo
        image.z = 0.0  # todo

        # \Date: 02:18:33 PM Wed Sep 07 2022  (is this start or end time?)
        if haskey(image.header, "Date")
            image.start_time = DateTime(image.header["Date"], DateFormat("HH:MM:SS p eee u dd yyyy"))
        end
        
        if haskey(image.header, "Ciao image 1 - Relative frame time")
            image.acquisition_time = parse(Float64, image.header["Ciao image 1 - Relative frame time"])
        end

        channel_info = get_channel_info_spm(image)
        channel_names = [info["name"] for info in channel_info]
        channel_units = [info["unit"] for info in channel_info]
        idx_unique = unique(i -> channel_names[i], eachindex(channel_names))
        image.channel_names = channel_names[idx_unique]
        image.channel_units = channel_units[idx_unique]

        if !startswith(image.header["Start context"], "OL")
            println("Warning: Unsupported file type. Expected 'Start context' of 'OLxxx' . Skipping.")
            return
        end

        # read body
        if !header_only
            output_info > 0 && println("Reading body of $(image.filename)")

            num_channels = length(channel_info)
            x_pixels, y_pixels = image.pixelsize
            image.data = Array{Float32}(undef, x_pixels, y_pixels, num_channels)

            image.channel_indices_fwd = fill(0, length(image.channel_names))
            image.channel_indices_bwd = fill(0, length(image.channel_names))

            for (i_info, info) in enumerate(channel_info)  # backward and forward channels
                i_ch = findfirst(x -> x == info["name"], image.channel_names)
                if info["direction"] == fwd
                    image.channel_indices_fwd[i_ch] = i_info
                else
                    image.channel_indices_bwd[i_ch] = i_info
                end
                skip(f, info["data_offset"])

                x_pixels = info["x_pixels"]  # in case we have a partial image
                y_pixels = info["y_pixels"]  # in case we have a partial image
                if info["data_length"] != x_pixels * y_pixels * 4
                    println("Warning: Data length does not match expected length for channel $(info["name"]).")
                    image.data[:,:,i_info] .= NaN32
                    continue
                end
                data_int = Array{Int32}(undef, x_pixels, y_pixels)  # this seems to be always Int16, despite the parameter 'bytes_per_pixel'
                seek(f, info["data_offset"])
                read!(f, data_int)
                data_int = ltoh.(data_int)  # little-endian (LSB) to host endian
                image.data[1:x_pixels,1:y_pixels,i_info] = (data_int .* info["calibration_factor"] .+ info["calibration_offset"]) .* info["calibration_factor2"]
                if x_pixels != image.pixelsize[1] ||  y_pixels != image.pixelsize[1]   # partial image
                    image.data[x_pixels + 1:image.pixelsize[1], y_pixels + 1:image.pixelsize[2], i_info] .= NaN32
                end
    
            end
        end
    end
    return image
end
