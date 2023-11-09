const GXSM_FILENAME_ATTRIB_SEPARATOR = "-"
const GXSM_FORWARD_SCAN_DIR = "Xp"
const GXSM_BACKWARD_SCAN_DIR = "Xm"
const GXSM_BACKWARD_MAXLENGTH_HEADER = 10

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
   
    NetCDF.open(fnames[1]) do nc
        # read header data

        # global attributes
        for (k,v) in nc.gatts
            break
            image.header[k] = v
        end

        for name in keys(nc.vars)
            break
            var = nc.vars[name]
            if haskey(var.atts, "long_name")
                key = var.atts["long_name"]
            elseif haskey(var.atts, "short_name")
                key = var.atts["short_name"]
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
                if haskey(nc.dim, name)
                    val = var[1:nc.dim[name].dimlen]
                elseif haskey(nc.dim, name * "_dim")
                    val = var[1:nc.dim[name * "_dim"].dimlen]
                else
                    val = var[:]
                end
                image.header[key] = strip(replace(String(val), '\0' => ' '))
            elseif typeof(var) <: NcVar{<:Number}
                if haskey(nc.dim, name)
                    if nc.dim[name].dimlen > GXSM_BACKWARD_MAXLENGTH_HEADER
                        val = var[1:GXSM_BACKWARD_MAXLENGTH_HEADER]
                        image.header[key] = string(val) * " ..."
                    else
                        val = var[1:nc.dim[name].dimlen]
                        image.header[key] = string(val)
                    end
                elseif length(var) > GXSM_BACKWARD_MAXLENGTH_HEADER
                    val = var[1:GXSM_BACKWARD_MAXLENGTH_HEADER]
                    image.header[key] = string(val) * " ..."
                elseif length(var) > 1
                    val = var[:]
                    image.header[key] = string(val)
                else
                    val = var[1]
                    image.header[key] = string(val)
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
        if haskey(image.header, "spm_scancontrol")
            image.scan_direction = image.header["spm_scancontrol"] == "TowDown" ? down : up
        end
        if "sranger_mk2_hwi_bias" in keys(nc.vars)
            image.bias = nc.vars["sranger_mk2_hwi_bias"][1]
        end
        
        # todo: z-controller/feedback
        # z_controller_data = split(image.header["Z-controller"], "\t")
        # image.z_feedback = z_controller_data[8] == "1" ? true : false
        # image.z_feedback_setpoint = parse(Float64, split(z_controller_data[9])[1])
        # image.z_feedback_setpoint_unit = split(z_controller_data[9])[2]

        image.z = nc.vars["dz"][1] .* 0.1  # 0.1 is to convert to nm

        image.start_time = unix2datetime(nc.vars["t_start"][1])
        image.acquisition_time = nc.vars["time"][1]

        return image

        # todo: get channel names and units

        # todo: read data

    end
end
