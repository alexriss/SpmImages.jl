const GXSM_FILENAME_ATTRIB_SEPARATOR = "-"
const GXSM_FORWARD_SCAN_DIR = "Xp"
const GXSM_BACKWARD_SCAN_DIR = "Xm"
const GXSM_MAXLENGTH_HEADER = 20



"""
    get_channel_names_units_netCDF(fnames::Vector{String})

Gets a list of all channel names and units in `image`.
"""
function get_channel_names_units_netCDF(fnames::Vector{String})
    names = Vector{String}(undef, 0)
    units = Vector{String}(undef, 0)
    files_fwd = Dict{String, String}()
    files_bwd = Dict{String, String}()

    for fname in fnames
        fbase, ext = rsplit(fname, "."; limit=2)
        if ext != "nc"
            error("File ($fname) has the wrong extension.")
            continue
        end
        entries = split(fbase, GXSM_FILENAME_ATTRIB_SEPARATOR)
        if length(entries) < 3
            error("File ($fname) has the wrong naming format.")
            continue
        else
            if entries[end] ∉ names  # there are backwards and forward channels, we only add once
                push!(names, entries[end])
                push!(units, "V")  # todo: how can I know the unit?
            end

            # create list of files for fwd and bwd directions
            if entries[end-1] == GXSM_BACKWARD_SCAN_DIR
                files_bwd[entries[end]] = fname
            elseif entries[end-1] == GXSM_FORWARD_SCAN_DIR
                files_fwd[entries[end]] = fname
            elseif !haskey(files_fwd, entries[end])
                files_fwd[entries[end]] = fname
            end
        end
    end
    return names, units, files_fwd, files_bwd
end


"""
    load_image_netCDF(fnames::Vector{String}, output_info::Int=1, header_only::Bool=false)

Loads data from `fnames`, specifying a list of file names of a GSXM netCDF files, and returns an `SpmImage` object.
The first file is used to read the header data, the image data is read from all given files.

`output_info` specifies the amount of output info to print to stdout when reading the files. 0 for no output, 1 for limited output, 2 for detailed output.
If `header_only` is `true` then only header data is read, image data is ignored.
"""
function load_image_netCDF(fnames::Vector{String}, output_info::Int=1, header_only::Bool=false)
    image = length(fnames) == 1 ? SpmImage(fnames[1], nc) : SpmImage(fnames, nc)

    output_info > 0 && println("Reading header of $(image.filename)")

    files_fwd = Dict{String, String}()
    files_bwd = Dict{String, String}()
   
    NetCDF.open(fnames[1]) do nc
        # read header data

        # global attributes
        for (k,v) in nc.gatts
            image.header[k] = v
        end

        for name in keys(nc.vars)
            var = nc.vars[name]
            if haskey(var.atts, "short_name")
                key = var.atts["short_name"]
            elseif haskey(var.atts, "long_name")
                key = var.atts["long_name"]
            elseif haskey(var.atts, "label")
                key = var.atts["label"] * " (" * name * ")"
            else
                key = name
            end
            unit = get(var.atts, "var_unit", "")
            if unit == "1"
                key *= " []"
            elseif unit != ""
                key *= " [" * unit * "]"
            end

            # check if key already exists
            if haskey(image.header, key)
                i_key = 1
                while haskey(image.header, key * " $(i_key)")
                    i_key += 1
                end
                key *= " $(i_key)"
            end

            if typeof(var) <: NcVar{NetCDF.ASCIIChar}
                if var.ndim == 0
                    val = var[:]
                elseif var.ndim == 1
                    val = var[1:var.dim[1].dimlen]
                else
                    val = var[:]
                end
               image.header[key] = strip(replace(String(val), '\0' => ' '))
            elseif typeof(var) <: NcVar{<:Number}
                if var.ndim == 0
                    val = var[1]
                    image.header[key] = string(val)
                elseif var.ndim == 1
                    l = var.dim[1].dimlen
                    if l == 1
                        val = var[1]
                        image.header[key] = string(val)
                    elseif l > GXSM_MAXLENGTH_HEADER
                        val = var[1:GXSM_MAXLENGTH_HEADER]
                        image.header[key] = string(val) * " ..."
                    else
                        val = var[1:l]
                        image.header[key] = string(val)
                    end
                # we skip everything that is > 1 dim
                end
            elseif output_info > 1
                println("Warning: unknown type of variable $(name): $(typeof(var))")
            end
        end

        # parse some of the extracted data
        image.scansize = [nc.vars["rangex"][1], nc.vars["rangey"][1]] .* 0.1  # 0.1 is to convert to nm
        image.scansize_unit = "nm"
        # todo: is this the center or the start of the scan?
        image.center = [nc.vars["dx"][1], nc.vars["dy"][1]] .* 0.1  # 0.1 is to convert to nm
        # todo: is rotation clockwise or counterclockwise?
        image.angle = nc.vars["alpha"][1]
        image.pixelsize = [nc.dim["dimx"].dimlen, nc.dim["dimy"].dimlen]
        if haskey(image.header, "spm_scancontrol: scan direction")
            image.scan_direction = image.header["spm_scancontrol: scan direction"] == "TopDown" ? down : up
        end
        if "sranger_mk2_hwi_bias" in keys(nc.vars)
            image.bias = nc.vars["sranger_mk2_hwi_bias"][1]
        end
        
        # todo: z-controller/feedback
        # image.z_feedback = z_controller_data[8] == "1" ? true : false
        # header -  ncvars
        # Current Setpt. [nA] - sranger_mk2_hwi_mix0_current_set_point
        # Voltage Setpt. [nA] - sranger_mk2_hwi_mix1_voltage_set_point
        # Aux2 Setpt. [nA] - sranger_mk2_hwi_mix2_aux2_set_point
        # Aux3 Setpt. [nA] - sranger_mk2_hwi_mix3_aux3_set_point
        if haskey(image.header, "Current Setpt. [nA]")
            image.z_feedback_setpoint = parse(Float64, image.header["Current Setpt. [nA]"]) .* 1e-9  # convert to A
            # todo: check if feedback is active
            # image.z_feedback_setpoint_unit = "A"
            image.z_feedback_setpoint_unit = ""  # this means that we don't know if it is active
        end

        # todo image.z

        image.start_time = unix2datetime(nc.vars["t_start"][1])
        image.acquisition_time = nc.vars["time"][1]

        image.channel_names, image.channel_units, files_fwd, files_bwd = get_channel_names_units_netCDF(fnames)
    end

    # read data
    if !header_only
        num_channels = length(image.channel_names) * 2    # the "*2" is there because of forward and backward channels (we assume so)
        x_pixels, y_pixels = image.pixelsize
        image.data = Array{Float32}(undef, x_pixels, y_pixels, num_channels)
        i = 1
        channel_indices_fwd = fill(0, length(image.channel_names))
        channel_indices_bwd = fill(0, length(image.channel_names))
        for (i_ch, ch) in enumerate(image.channel_names)
            if ch in keys(files_fwd)
                fname = files_fwd[ch]
                output_info > 0 &&  println("Reading body of $(image.filename)")
                image.data[:,:,i] = NetCDF.ncread(fname, "FloatField")[:,:,1,1]
                channel_indices_fwd[i_ch] = i
                i += 1
            end
            if ch in keys(files_bwd)
                fname = files_bwd[ch]
                output_info > 0 &&  println("Reading body of $(image.filename)")
                image.data[:,:,i] = NetCDF.ncread(fname, "FloatField")[:,:,1,1]
                channel_indices_bwd[i_ch] = i
                i += 1
            end
        end
        image.channel_indices_fwd = channel_indices_fwd
        image.channel_indices_bwd = channel_indices_bwd
    end
    return image
end


