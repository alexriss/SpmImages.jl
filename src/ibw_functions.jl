# Loads Igor binary wave files (.ibw)
#
# resources:
#   https://www.wavemetrics.net/Downloads/FTP_Archive/IgorPro/Technical_Notes/index.html
#     see TN #003: Igor binary file format and PTN #003: Writing Packed Files
#   https://de.mathworks.com/matlabcentral/fileexchange/42679-igor-pro-file-format-ibw-to-matlab-variable
#   https://blog.tremily.us/posts/igor/
#   https://github.com/AFM-analysis/igor2/tree/master/igor2
#   https://www.wavemetrics.com/forum/general/igor-binary-files-image-processing-asylum-research-data


const IBW_MAXDIMS = 4
const IBW_MAX_WAVE_NAME2 = 18  # Maximum length of wave name in version 1 and 2 files. Does not include the trailing null.
const IBW_MAX_WAVE_NAME5 = 31  # Maximum length of wave name in version 5 files. Does not include the trailing null.
const IBW_MAX_UNIT_CHARS = 3
const IBW_MAXLENGTH_HEADER = 128

const TYPE_TABLE = Dict(
    0 => nothing,     # Text wave, not handled in ReadWave.c
    1 => nothing,  # NT_CMPLX, makes number complex.
    2=> Float32,  # NT_FP32, 32 bit fp numbers.
    3=> Complex{Float64},
    4=> Float64,  # NT_FP64, 64 bit fp numbers.
    5=> Complex{Float64},
    8=> Int8,    # NT_I8, 8 bit signed integer. Requires Igor Pro
                      # 2.0 or later.
    9=> Complex{Int8},
    0x10=> Int16,  # NT_I16, 16 bit integer numbers. Requires Igor
                      # Pro 2.0 or later.
    0x11=> Complex{Int16},
    0x20=> Int32,  # NT_I32, 32 bit integer numbers. Requires Igor
                      # Pro 2.0 or later.
    0x21=> Complex{Int32},
    #   0x40=>None,        # NT_UNSIGNED, Makes above signed integers
    #                     # unsigned. Requires Igor Pro 3.0 or later.
    0x48=> UInt8,
    0x49=> Complex{UInt8},
    0x50=> UInt16,
    0x51=> Complex{UInt16},
    0x60=> UInt32,
    0x61=> Complex{UInt32}
)


@io struct BinHeader1 
	version::Int16   # Version number for backwards compatibility.
	wfmSize::Int32   # The size of the WaveHeader2 data structure plus the wave data plus 16 bytes of padding.
	checksum::Int16  # Checksum over this header and the wave header.
end align_packed

@io struct BinHeader2 
	version::Int16   # Version number for backwards compatibility.
	wfmSize::Int32   # The size of the WaveHeader2 data structure plus the wave data plus 16 bytes of padding.
	noteSize::Int32  # The size of the note text.
	pictSize::Int32  # Reserved. Write zero. Ignore on read.
	checksum::Int16  # Checksum over this header and the wave header.
end align_packed

@io struct BinHeader3 
	version::Int16      # Version number for backwards compatibility.
	wfmSize::Int32      # The size of the WaveHeader2 data structure plus the wave data plus 16 bytes of padding.
	noteSize::Int32     # The size of the note text.
	formulaSize::Int32  # The size of the dependency formula, if any.
	pictSize::Int32     # Reserved. Write zero. Ignore on read.
	checksum::Int16     # Checksum over this header and the wave header.
end align_packed

@io struct BinHeader5 
	version::Int16                 # Version number for backwards compatibility.
	checksum::Int16                # Checksum over this header and the wave header.
	wfmSize::Int32                 # The size of the WaveHeader5 data structure plus the wave data.
	formulaSize::Int32             # The size of the dependency formula, if any.
	noteSize::Int32                # The size of the note text.
	dataEUnitsSize::Int32          # The size of optional extended data units.
	dimEUnitsSize::SVector{IBW_MAXDIMS,Int32}  # The size of optional extended dimension units.
	dimLabelsSize::SVector{IBW_MAXDIMS,Int32}  # The size of optional dimension labels.
	sIndicesSize::Int32            # The size of string indicies if this is a text wave.
	optionsSize1::Int32            # Reserved. Write zero. Ignore on read.
	optionsSize2::Int32            # Reserved. Write zero. Ignore on read.
end align_packed

@io struct WaveHeader2
	type::Int16							# See types (e.g. NT_FP64) above. Zero for text waves.
	WaveHeader2::Int32			# Used in memory only. Write zero. Ignore on read.

	bname::SVector{IBW_MAX_WAVE_NAME2+2,UInt8} 	# Name of wave plus trailing null.
	whVersion::Int16					# Write 0. Ignore on read.
	srcFldr::Int16						# Used in memory only. Write zero. Ignore on read.
	fileName::Int32					# Used in memory only. Write zero. Ignore on read.

	dataUnits::SVector{IBW_MAX_UNIT_CHARS+1,UInt8}	# Natural data units go here - null if none.
	xUnits::SVector{IBW_MAX_UNIT_CHARS+1,UInt8}		# Natural x-axis units go here - null if none.

	npnts::Int32							# Number of data points in wave.

	aModified::Int16					# Used in memory only. Write zero. Ignore on read.
	hsA::Float64						# X value for point p = hsA*p + hsB
    hsB::Float64						

	wModified::Int16					# Used in memory only. Write zero. Ignore on read.
	swModified::Int16					# Used in memory only. Write zero. Ignore on read.
	fsValid::Int16						# True if full scale values have meaning.
	topFullScale::Float64 	# The min full scale value for wave.
    botFullScale::Float64 
		   
	useBits::UInt8						# Used in memory only. Write zero. Ignore on read.
	kindBits::UInt8						# Reserved. Write zero. Ignore on read.
	formula::Int32						# Used in memory only. Write zero. Ignore on read.
	depID::Int32							# Used in memory only. Write zero. Ignore on read.
	creationDate::UInt32			# DateTime of creation. Not used in version 1 files.
	wUnused::SVector{2, UInt8}					# Reserved. Write zero. Ignore on read.

	modDate::UInt32				# DateTime of last modification.
	waveNoteH::Int32					# Used in memory only. Write zero. Ignore on read.

	wData::SVector{4,Float32}						# The start of the array of waveform data.
end align_packed

@io struct WaveHeader5
	WaveHeader5::Int32			# link to next wave in linked list.

	creationDate::UInt32			# DateTime of creation.
	modDate::UInt32				# DateTime of last modification.

	npnts::Int32							# Total number of points (multiply dimensions up to first zero).
	type::Int16							# See types (e.g. NT_FP64) above. Zero for text waves.
	dLock::Int16						# Reserved. Write zero. Ignore on read.

	whpad1::SVector{6,UInt8}						# Reserved. Write zero. Ignore on read.
	whVersion::Int16					# Write 1. Ignore on read.
	bname::SVector{IBW_MAX_WAVE_NAME5+1,UInt8}		# Name of wave plus trailing null.
	whpad2::Int32						# Reserved. Write zero. Ignore on read.
	DataFolder::Int32		# Used in memory only. Write zero. Ignore on read.

	# Dimensioning info. [0] == rows, [1] == cols etc
	nDim::SVector{IBW_MAXDIMS,Int32}					# Number of of items in a dimension -- 0 means no data.
	sfA::SVector{IBW_MAXDIMS,Float64}				# Index value for element e of dimension d = sfA[d]*e + sfB[d].
	sfB::SVector{IBW_MAXDIMS,Float64}

	# SI units
	dataUnits::SVector{IBW_MAX_UNIT_CHARS+1,UInt8}			# Natural data units go here - null if none.
	dimUnits::SVector{IBW_MAXDIMS * (IBW_MAX_UNIT_CHARS+1),UInt8}	# Natural dimension units go here - null if none.

	fsValid::Int16						# TRUE if full scale values have meaning.
	whpad3::Int16						# Reserved. Write zero. Ignore on read.
	topFullScale::Float64	# The max and max full scale value for wave.
    botFullScale::Float64

	dataEUnits::Int32					# Used in memory only. Write zero. Ignore on read.
	dimEUnits::SVector{IBW_MAXDIMS,Int32}			# Used in memory only. Write zero. Ignore on read.
	dimLabels::SVector{IBW_MAXDIMS,Int32}			# Used in memory only. Write zero. Ignore on read.
	
	waveNoteH::Int32					# Used in memory only. Write zero. Ignore on read.
	whUnused::SVector{16,Int32}					# Reserved. Write zero. Ignore on read.

	# The following stuff is considered private to Igor.

	aModified::Int16					# Used in memory only. Write zero. Ignore on read.
	wModified::Int16					# Used in memory only. Write zero. Ignore on read.
	swModified::Int16					# Used in memory only. Write zero. Ignore on read.
	
	useBits::UInt8						# Used in memory only. Write zero. Ignore on read.
	kindBits::UInt8						# Reserved. Write zero. Ignore on read.
	formula::Int32						# Used in memory only. Write zero. Ignore on read.
	depID::Int32							# Used in memory only. Write zero. Ignore on read.
	
	whpad4::Int16						# Reserved. Write zero. Ignore on read.
	srcFldr::Int16						# Used in memory only. Write zero. Ignore on read.
	fileName::Int32					# Used in memory only. Write zero. Ignore on read.
	
	sIndices::Int32					# Used in memory only. Write zero. Ignore on read.

	wData::Float32						# The start of the array of data. Must be 64 bit aligned.
end align_packed

BinHeaderAll = Union{BinHeader1, BinHeader2, BinHeader3, BinHeader5}
WaveHeaderAll = Union{WaveHeader2, WaveHeader5}

header_skip_fields = [
    :pictSize, :optionsSize1, :optionsSize2,
    :WaveHeader2, :whVersion, :srcFldr, :fileName, :aModified, :wModified, :swModified, :useBits, :kindBits, :formula, :depID, :wUnused, :waveNoteH,
    :dLock, :whPad1, :whPad2, :DataFolder, :whpad3, :dataEUnits, :dimEUnits, :dimLabels, :whUnused, :whpad4, :sIndices,
    :wData
]

"""
    load_image_ibw(fname::String, output_info::Int=1, header_only::Bool=false; extra_checks::Bool=false)

Loads data from `fname`, specifying the file name of a Nanonis .sxm file, and returns an `SpmImage` object.

`output_info` specifies the amount of output info to print to stdout when reading the files. 0 for no output, 1 for limited output, 2 for detailed output.
If `header_only` is `true` then only header data is read, image data is ignored.
If `extra_checks` is `true` then additional checks are performed to ensure that the file is not corrupted.
"""
function load_image_ibw(fname::String, output_info::Int=1, header_only::Bool=false; extra_checks::Bool=false)
    image = SpmImage(fname, ibw)

    if output_info > 0
        println("Reading header of $(image.filename)")
    end

    creationDate = missing
    dimLabels = String[]
    dimUnits = String[]
    open(image.filename) do f
        version = peek(f, UInt16)

        switch_endian = (version & 0xFF) == 0
        if switch_endian
            if ENDIAN_BOM == 0x01020304 # we are on big endian
                from_endian = :LittleEndian
            else
                from_endian = :BigEndian
            end
            version = bswap(version)
        else
            if ENDIAN_BOM == 0x01020304 # we are on big endian
                from_endian = :BigEndian
            else
                from_endian = :LittleEndian
            end
        end

        checksum_diff = 0
        if version == 1
            BinHeader = BinHeader1
            WaveHeader = WaveHeader2
        elseif version == 2
            BinHeader = BinHeader2
            WaveHeader = WaveHeader2
        elseif version == 3
            BinHeader = BinHeader3
            WaveHeader = WaveHeader2
        elseif version == 5
            BinHeader = BinHeader5
            WaveHeader = WaveHeader5
            checksum_diff = -4	# Version 5 checksum does not include the wData field.
        else
            error("Unknown version of $(image.filename): $version.")
        end

        binHeader = unpack(f, BinHeader, from_endian)
        waveHeader = unpack(f, WaveHeader, from_endian)

        if extra_checks
            if checksum(binHeader, waveHeader, 0, checksum_diff, from_endian, switch_endian) !== 0
                error("Checksum error in $(image.filename).")
            end
        end

        creationDate = unix2datetime(waveHeader.creationDate + 2082844800)  # Mac HFS+ timestamp (seconds since 1904)
        modDate = unix2datetime(waveHeader.modDate + 2082844800)  # Mac HFS+ timestamp (seconds since 1904)

        for header in (binHeader, waveHeader)
            for f in fieldnames(typeof(header))
                f in header_skip_fields && continue
                val = getfield(header, f)
                if typeof(val) <: SArray{<:Tuple,UInt8,1}
                        val_str = char_to_string(val)
                elseif typeof(val) <: SArray{<:Tuple,Int32,1}
                    val_str = join(string.(val), ", ")
                else
                    val_str = string(val)
                end
                length(val_str) > IBW_MAXLENGTH_HEADER && (val_str = val_str[1:IBW_MAXLENGTH_HEADER] * " ...")
                image.header[string(f)] = val_str
            end
        end

        # set file pointer to the beginning of the wData field
        skip(f, -sizeof(waveHeader.wData))
        if !header_only
            output_info > 0 && println("Reading body of $(image.filename)")
            wave_data = load_numeric_wave_data(f, binHeader, waveHeader, switch_endian)
            if :nDim ∉ fieldnames(typeof(waveHeader))
                error("No image dimensions specified. This might not be an image file.")
            end
            nDim = [n for n in waveHeader.nDim if n > 0]
            length(wave_data) != prod(nDim) && error("Unexpected number of data points: $(length(wave_data)) (read from file) != $(join(nDim, 'x')) (dimensions specified).")
            image.data = reshape(wave_data, nDim...)
        else
            skip_numeric_wave_data(f, waveHeader)
        end

        if version == 2 || version == 3
            skip(f, 16)
        elseif version == 5
            skip(f, binHeader.formulaSize)
        end
        noteArr = Array{UInt8}(undef, binHeader.noteSize)
        read!(f, noteArr)
        notes = split(String(noteArr), x -> x in ('\r', '\n'))
        for note in notes
            if occursin(":", note)
                key, val = split(note, ":")
                key = get_unique_key(image.header, strip(key))
                image.header[key] = strip(val)
            end
        end
        if version == 5
            dataUnits = ""
            if binHeader.dataEUnitsSize > 0
                dataUnitsArr = Array{UInt8}(undef, binHeader.dataEUnitsSize)
                read!(f, dataUnitsArr)
                dataUnits = char_to_string(dataUnitsArr)
            end
            image.header["dimEUnitsSize"] = join(dataUnits, ", ")
            dimEUnits = map(binHeader.dimEUnitsSize) do s
                s == 0 && return ""
                arr = Array{UInt8}(undef, s)
                read!(f, arr)
                char_to_string(arr)
            end 
            image.header["dimEUnits"] = join(dimEUnits, ", ")

            dimLabelsRaw = map(binHeader.dimLabelsSize) do s
                s == 0 && return ""
                arr = Array{UInt8}(undef, s)
                read!(f, arr)
                String(arr)
            end
            dimLabels, dimUnits = get_dim_labels_units(dimLabelsRaw)
        end

        # number of channels should match the number of dimensions
        nDim = [n for n in waveHeader.nDim if n > 0]
        length(nDim) != 3 && error("Unexpected number of dimensions: $(length(nDim)) (read from file) != 3 (expected).")
        nDim[3] < 1 && error("Channels specified in file dimension should be greater than 0 (found $(nDim[3])).")
        i = 1
        while length(dimLabels) < nDim[3]
            push!(dimLabels, "Channel $i")
            push!(dimUnits, "V")
        end

        # todo: backward and forward channels should exist for each channel
    end

    # parse some of the extracted data
    if "FastScanSize" in keys(image.header) && "SlowScanSize" in keys(image.header)
        image.scansize = parse.(Float64, [image.header["FastScanSize"], image.header["SlowScanSize"]]) * 1e9
    elseif "ScanSize" in keys(image.header)
        s = parse(Float, image.header["ScanSize"]) * 1e9
        image.scansize = [s, s]
    end
    image.scansize_unit = "nm"

    if "XOffset" in keys(image.header) && "YOffset" in keys(image.header)
        image.center = parse.(Float64, [image.header["XOffset"], image.header["YOffset"]]) * 1e9
    end
    if "ScanAngle" in keys(image.header)
        image.angle = parse(Float64, image.header["ScanAngle"])
    end
    if "ScanPoints" in keys(image.header) && "ScanLines" in keys(image.header)
        image.pixelsize = parse.(Int, [image.header["ScanPoints"], image.header["ScanLines"]]) 
    end
    if "BottomLine" in keys(image.header) && "TopLine" in keys(image.header)
        image.scan_direction = parse(Int, image.header["BottomLine"]) < parse(Int, image.header["TopLine"]) ? up : down
    end
    if "BiasVoltage" in keys(image.header)
        image.bias = parse(Float64, image.header["BiasVoltage"])
    end

    # todo: z_feedback

    if "Seconds" in keys(image.header)
        image.start_time = unix2datetime(parse(Float64, image.header["Seconds"]) + 2082844800)  # Mac HFS+ timestamp (seconds since 1904)
    else 
        image.start_time = modDate
    end

    # todo: acquisition_time, can maybe be calculated from scan size and scan speed

    image.channel_names = dimLabels
    image.channel_units = dimUnits

    return image
end


"""Calculates checksum of `binHeader` and `waveHeader`."""
function checksum(binHeader::BinHeaderAll, waveHeader::WaveHeaderAll, oldcksum::Int, checksum_diff::Int, from_endian::Symbol, switch_endian::Bool)
    for data in (binHeader, waveHeader)
        buf = IOBuffer()
        pack(buf, data, from_endian)
        seekstart(buf)
        data_int16 = read(buf, SVector{packed_sizeof(typeof(data)) ÷ 2,Int16})

        numbytes = sizeof(data_int16)
        data == waveHeader && (numbytes += checksum_diff) # Version 5 checksum does not include the wData field.
        numbytes >>= 1  # 2 bytes to a short -- ignore trailing odd byte.
        if switch_endian
            oldcksum += sum(bswap, view(data_int16, 1:numbytes))
        else
            oldcksum += sum(view(data_int16, 1:numbytes))
        end
    end
    return oldcksum & 0xffff
end


"""loads the numeric wave data from the file `f`."""
function load_numeric_wave_data(f::IO, binHeader::BinHeaderAll, waveHeader::WaveHeaderAll, switch_endian::Bool)
    T = TYPE_TABLE[waveHeader.type]
    isnothing(T) && error("Wave type $(waveHeader.type) is not supported.")

    npnts = waveHeader.npnts  # total number of elements in all dimensions.
    npnts == 0 && return Array{T}(undef, 0)
    wave_data = Array{T}(undef, npnts)

    WaveHeader = typeof(waveHeader)
    offset = sum(sizeof,WaveHeader.types[1:end-1])
    waveDataSize = binHeader.wfmSize - offset  # number of data bytes stored in the file
    binHeader.version != 5 && (waveDataSize -= 16)
    if (waveDataSize < sizeof(wave_data))  # some sort of dependency formula thing
        @warn("Unexpected data size: $waveDataSize < $(sizeof(wave_data)). Returning empty array.")
        fill!(wave_data, zero(T))
        return wave_data
    end

    read!(f, wave_data)
    switch_endian && (wave_data = bswap.(wave_data))
    return wave_data
end


"""skips the numeric wave data in the file `f`."""
function skip_numeric_wave_data(f::IO, waveHeader::WaveHeaderAll)
    T = TYPE_TABLE[waveHeader.type]
    isnothing(T) && error("Wave type $(waveHeader.type) is not supported.")

    npnts = waveHeader.npnts  # total number of elements in all dimensions.
    skip(f, npnts * sizeof(T))
    
    return nothing
end


"""converts Array of UInt8 into a string, stopping at the first null character."""
function char_to_string(c::Union{SArray{<:Tuple,UInt8,1},Vector{UInt8}})
    pos = findfirst(isequal(0), c)
    if isnothing(pos)
        return String(c)
    else
        return String(c[1:pos-1])
    end
end


"""gets a unique key for a dictionary, by appending a number to the key if it already exists."""
function get_unique_key(dict::AbstractDict, key::Union{String,SubString})
    if haskey(dict, key)
        i_key = 1
        while haskey(dict, key * " $(i_key)")
            i_key += 1
        end
        key *= " $(i_key)"
    end
    return key
end


"""gets the dimension labels and units from the raw strings."""
function get_dim_labels_units(rawLabels::SArray{<:Tuple,String})
    dimLabels = Vector{String}(undef, 0)
    dimUnits = Vector{String}(undef, 0)
    for s in rawLabels
        for start in 1:32:length(s)
            stop = min(start+31, length(s))    # labels are saved in chunks of 32 bytes
            label = s[start:stop]
            pos0 = findfirst(isequal('\0'), label)
            if !isnothing(pos0)
                pos0 == 1 && continue
                label = label[1:pos0-1]
            end
            endswith(label, "Trace") && (label = label[1:end-5])
            endswith(label, "Retrace") && (label = label[1:end-7] * " bwd")
            push!(dimLabels, label)
        end
    end
    return dimLabels, get_dim_units(dimLabels)
end


"""
Gets the dimension units from the dimension labels.
This is how it is done in Gwyddion:
https://github.com/christian-sahlmann/gwyddion/blob/master/modules/file/igorfile.c
in the function `channel_title_to_units`.
"""
function get_dim_units(dimLabels::Vector{String})
    dimUnits = map(dimLabels) do l
        startswith(l, "DAC") && (l = l[4:end])
        startswith(l, "Height") && return "m"
        startswith(l, "ZSensor") && return "m"
        startswith(l, "Deflection") && return "m"
        startswith(l, "Amplitude") && return "m"
        startswith(l, "Phase") && return "deg"
        startswith(l, "Current") && return "A"
        startswith(l, "Frequency") && return "Hz"
        startswith(l, "Capacitance") && return "F"
        startswith(l, "Potential") && return "V"
        startswith(l, "Count") && return ""
        startswith(l, "QFactor") && return ""
        # Everything else is in Volts.
        return "V"
    end 
    return dimUnits
end