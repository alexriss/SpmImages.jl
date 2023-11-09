const GXSM_FILENAME_ATTRIB_SEPARATOR = "-"
const GXSM_FORWARD_SCAN_DIR = "Xp"
const GXSM_BACKWARD_SCAN_DIR = "Xm"
const GXSM_BACKWARD_MAXLENGTH_HEADER = 20



"""
    get_channel_names_units_netCDF(image::SpmImage)

Gets a list of all channel names and units in `image`.
"""
function get_channel_names_units_netCDF(image::SpmImage, fnames::Vector{String})
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

                if entries[end-1] == GXSM_BACKWARD_SCAN_DIR
                    files_bwd[entries[end]] = fname
                elseif entries[end-1] == GXSM_FORWARD_SCAN_DIR
                    files_fwd[entries[end]] = fname
                elseif !haskey(files_fwd, entries[end])
                    files_fwd[entries[end]] = fname
                end
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
    image = SpmImage(fnames[1], nc)

    if output_info > 0
        println("Reading header of $(image.filename)")
    end

    files_fwd = Dict{String, String}()
    files_bwd = Dict{String, String}()
   
    NCDataset(fnames[1]) do ds
        # read header data

        # global attributes
        for (k,v) in ds.attrib
            image.header[k] = v
        end

        for name in keys(ds)
            var = ds[name].var
            if haskey(var.attrib, "short_name")
                key = var.attrib["short_name"]
            elseif haskey(var.attrib, "long_name")
                key = var.attrib["long_name"]
            elseif haskey(var.attrib, "label")
                key = var.attrib["label"] * " (" * name * ")"
            else
                key = name
            end
            unit = get(var.attrib, "var_unit", "")
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

            if eltype(var) <: Complex
                if ndims(var) == 0
                    val = var[:]
                elseif ndims(var) == 1
                    val = var[1:size(var, 1)]
                else
                    val = var[:]
                end
               image.header[key] = strip(replace(String(val), '\0' => ' '))
            elseif eltype(var) <: Complex
                if ndims(var) == 0
                    val = var[1]
                    image.header[key] = string(val)
                elseif ndims(var) == 1
                    l = size(var, 1)
                    if l == 1
                        val = var[1]
                        image.header[key] = string(val)
                    elseif l > GXSM_BACKWARD_MAXLENGTH_HEADER
                        val = var[1:GXSM_BACKWARD_MAXLENGTH_HEADER]
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
        image.scansize = [ds["rangex"][1], ds["rangey"][1]] .* 0.1  # 0.1 is to convert to nm
        image.scansize_unit = "nm"
        # todo: is this the center or the start of the scan?
        image.center = [ds["dx"][1], ds["dy"][1]] .* 0.1  # 0.1 is to convert to nm
        # todo: is rotation clockwise or counterclockwise?
        image.angle = ds["alpha"][1]
        image.pixelsize = [ds.dim["dimx"],  ds.dim["dimy"]]
        if haskey(image.header, "spm_scancontrol")
            image.scan_direction = image.header["spm_scancontrol"] == "TopDown" ? down : up
        end
        if "sranger_mk2_hwi_bias" in keys(ds)
            image.bias = ds["sranger_mk2_hwi_bias"][1]
        end
        
        # todo: z-controller/feedback
        # z_controller_data = split(image.header["Z-controller"], "\t")
        # image.z_feedback = z_controller_data[8] == "1" ? true : false
        # image.z_feedback_setpoint = parse(Float64, split(z_controller_data[9])[1])
        # image.z_feedback_setpoint_unit = split(z_controller_data[9])[2]

        image.z = ds["dz"][1] .* 0.1  # 0.1 is to convert to nm

        image.start_time = unix2datetime(ds["t_start"][1])
        image.acquisition_time = ds["time"][1]

        image.channel_names, image.channel_units, files_fwd, files_bwd = get_channel_names_units_netCDF(image, fnames)
    end

    # read data
    if !header_only
        num_channels = length(image.channel_names) * 2    # the "*2" is there because of forward and backward channels (we assume so)
        x_pixels, y_pixels = image.pixelsize
        image.data = Array{Float32}(undef, x_pixels, y_pixels, num_channels)
        for (i,ch) in enumerate(image.channel_names)
            i_data_fwd = 2 * i - 1
            i_data_bwd = 2 * i
            if ch in keys(files_fwd)
                fname = files_fwd[ch]
                NCDatasets.load!(NCDataset(fname)["FloatField"].var, image.data[:,:,i_data_fwd], :, :, 1, 1)
            else
                image.data[:,:,i_data_fwd] .= fill(NaN32, x_pixels, y_pixels)
            end
            if ch in keys(files_bwd)
                fname = files_bwd[ch]
                output_info > 0 &&  println("Reading body of $(image.filename)")
                NCDatasets.load!(NCDataset(fname)["FloatField"].var, image.data[:,:,i_data_bwd], :, :, 1, 1)
            else
                image.data[:,:,i_data_bwd] .= fill(NaN32, x_pixels, y_pixels)
            end
        end
    end
    return image
end
