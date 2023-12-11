using SpmImages
using Test

@testset "File loading sxm" begin
    ima = load_image("Image_445.sxm")
    @test ima.filename == "Image_445.sxm"
    @test ima.channel_names == ["Z", "LI X", "LI Y", "Bias", "Current", "Phase", "Amplitude", "Frequency Shift", "Excitation"]
    @test ima.channel_units == ["m", "V", "V", "V", "A", "deg", "m", "Hz", "V"]
    @test ima.scansize == [2.0, 4.0]
    @test ima.scansize_unit == "nm"
    @test ima.pixelsize == [96, 192]
    @test ima.scan_direction == SpmImages.up
    @test ima.bias ≈ 5.000E-2

    @test ima.z_feedback_setpoint ≈ 1.0e-11
    @test ima.z_feedback_setpoint_unit == "A"
    @test ima.z ≈ -1.43049e-8

    @test ima.acquisition_time == 4877.1
    @test ima.header["Oscillation Control>Reference Phase (deg)"] == "65.16E+0"
    @test ima.data[12,42,3] ≈ -3.6583905f0
    @test ima.data[9,18,10] ≈ 3.0707468f-11
    @test ima.data[69,69,6] ≈ -0.7235545f0

    @test get_channel(ima, "Frequency shift").name == "Frequency Shift"
    @test get_channel(ima, "exCitation bwd").name == "Excitation"
    @test get_channel(ima, "exCitation bwd").direction == SpmImages.bwd
    @test get_channel(ima, "current").data[30, 55] ≈ 7.2952433f-10
    @test get_channel(ima, "current bwd").data[30, 55] ≈ 6.729436f-10
    @test get_channel(ima, "current").data[12, 90] ≈ 2.5479948f-11
    @test get_channel(ima, "current bwd").data[12, 90] ≈ 2.3916345f-11
    @test get_channel(ima, "Frequency Shift").data[78, 69] ≈ -2.944235f0
    @test get_channel(ima, "frequency shift [bwd] ").data[78, 49] ≈ -0.581533f0
    @test get_channel(ima, "frequency shift b ").data[78, 49] ≈ -0.581533f0
    @test get_channel(ima, "frequency shift backwards").data[78, 49] ≈ -0.581533f0
end

@testset "File loading nc" begin

    fnames = ["2023_12_11_007-M-Xp-Topo.nc", "2023_12_11_007-Xm-MixmIn-0mITunnel.nc", "2023_12_11_007-Xm-Topo.nc"]
    names, units, files_fwd, files_bwd = SpmImages.get_channel_names_units_netCDF(fnames)
    @test names == ["Topo", "0mITunnel"]
    @test units == ["V", "V"] 
    @test files_fwd == Dict("Topo" => fnames[1])
    @test files_bwd == Dict("0mITunnel" => fnames[2], "Topo" => fnames[3])

    fnames = ["2023_12_11_007-M-Xp-Topo.nc", "2023_12_11_007-Xm-0mITunnel.nc", "2023_12_11_007-Xm-Topo.nc"]
    names, units, files_fwd, files_bwd = SpmImages.get_channel_names_units_netCDF(fnames)
    @test names == ["Topo", "0mITunnel"]
    @test units == ["V", "V"] 
    @test files_fwd == Dict("Topo" => fnames[1])
    @test files_bwd == Dict("0mITunnel" => fnames[2], "Topo" => fnames[3])

    ima = load_image("nc/chigwell009-M-Xp-Topo.nc")
    @test ima.filename == "nc/chigwell009-M-Xp-Topo.nc"
    @test ima.channel_names == ["Topo"]
    @test ima.channel_units == ["V"]
    @test ima.scansize == [10000.0, 10000.0]
    @test ima.scansize_unit == "nm"
    @test ima.pixelsize == [64, 64]
    @test ima.scan_direction == SpmImages.down
    @test ima.bias ≈ -5.0

    fnames = filter(x -> occursin("chigwell", x), readdir("nc", join=true))
    ima = load_image(fnames)

    @test length(ima.filename) == length(fnames)
    @test ima.filename[1] == fnames[1]
    @test ima.filename[6] == fnames[6]
    @test ima.channel_names == ["Topo", "ADC1", "ADC2", "ADC4", "ADC6", "ADC7"]
    @test ima.channel_units == ["V", "V", "V", "V", "V", "V"]
    @test ima.scansize == [10000.0, 10000.0]
    @test ima.scansize_unit == "nm"
    @test ima.pixelsize == [64, 64]
    @test ima.scan_direction == SpmImages.down
    @test ima.bias ≈ -5.0

    @test ima.acquisition_time == 382.0
    @test ima.header["View Offset Z [Ang]"] == "0.0"
    @test ima.header["SRanger: Scan Event Xp Bias [Volt]"] == "-5.135874684612657e305"
    @test ima.header["t_end"] == "1638465030"

    @test ima.header["InstrumentName"] == "LT-AFM"
    @test ima.header["InstrumentType"] == "AFM"
    @test ima.header["HardwareCtrlType"] == "SRangerMK2:SPM"
    @test ima.header["Version"] == "3.48.0"
    @test ima.header["HardwareConnectionDev"] == "/dev/sranger_mk2_0"
    @test ima.header["username"] == "Nobody"
    @test ima.header["# Pixels in X, contains X-Pos Lookup"] == "Float32[-50000.0, -48412.7, -46825.4, -45238.094, -43650.793, -42063.492, -40476.19, -38888.89, -37301.586, -35714.285, -34126.984, -32539.682, -30952.38, -29365.08, -27777.777, -26190.477, -24603.174, -23015.873, -21428.572, -19841.27] ..."

    @test ima.data[12,42,3] ≈ 0.34779367f0
    @test ima.data[9,18,10] ≈ 0.051706143f0
    @test ima.data[60,60,6] ≈ 1.1762803f0

    @test get_channel(ima, "Topo").name == "Topo"
    @test get_channel(ima, "Adc1 bwd").name == "ADC1"
    @test get_channel(ima, "AdC1 backward  ").direction == SpmImages.bwd
    @test get_channel(ima, "topo").data[30, 55] ≈ -1104.3302f0
    @test get_channel(ima, "topo bwd").data[30, 55] ≈ -1082.6204f0
    @test get_channel(ima, "topo").data[12, 9] ≈ -4565.998f0
    @test get_channel(ima, "topo bwd").data[12, 9] ≈ -4533.6323f0
    @test get_channel(ima, "ADC2").data[8, 60] ≈ 1.1837914f0
    @test get_channel(ima, "adc2 [bwd] ").data[8, 49] ≈ 1.177694f0
    @test get_channel(ima, "adc4").data[49, 1] ≈ 1.8622661f0
    @test get_channel(ima, "adC4 backwards").data[49, 1] ≈ 1.8630849f0

    ima = load_image("nc/2023_12_11_009-M-Xp-Topo.nc")
    @test ima.angle ≈ 40.0
    @test get_channel(ima, "Topo").name == "Topo"
    @test ima.channel_names == ["Topo"]
    @test ima.channel_units == ["V"]
    @test get_channel(ima, "Topo").data[22, 50] ≈ -1255.0563f0
    @test isnan(get_channel(ima, "Topo bwd").data[22, 50])
end

@testset "Coordinate conversion" begin
    ima = load_image("Image_445.sxm")
    @test pixels_to_nm(ima, [1, 1]) ≈ [0., 0.]
    @test pixels_to_nm(ima, [96, 192]) ≈ [2.0, 4.0]
    @test nm_to_pixels(ima, [0., 0.]) ≈ [1., 1.]
    @test nm_to_pixels(ima, [2.0, 4.0]) ≈ [96, 192]
end

@testset "Background correction" begin
    d = [1 2 3; 2 3 4; 3 4 5]
    @test all(correct_background(d, plane_linear_fit) .< 2e-15)
    @test correct_background(d, line_average, false) == [-1 0 1; -1 0 1; -1 0 1]
    @test correct_background(d, vline_average, false) == [-1 -1 -1; 0 0 0; 1 1 1]
    @test correct_background(d, line_linear_fit, false) == [0 0 0; 0 0 0; 0 0 0]
    @test correct_background(d, vline_linear_fit, false) == [0 0 0; 0 0 0; 0 0 0]
    @test correct_background(d, line_average, true) == [0 1 2; 0 1 2; 0 1 2]
    @test correct_background(d, vline_average, true) == [0 0 0; 1 1 1; 2 2 2]
    @test correct_background(d, subtract_minimum, true) == [0 1 2; 1 2 3; 2 3 4]

    d = [1.0 2.0 3.0; 4.0 3.0 2.0; 3.0 4.0 5.0]
    @test correct_background(d, line_linear_fit) == [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    @test all(correct_background(d, vline_linear_fit, false) - [-0.666667 0.0 0.666667; 1.33333 0.0 -1.33333; -0.666667 0.0 0.666667] .< 1e-5)
end

@testset "Line profiles" begin
    # todo: check line profiles, also for non-square images and images with anisotropic pixel-densities
end

@testset "Drift" begin
    im1 = load_image("01514.sxm")
    im2 = load_image("01514.sxm")
    im2.angle = 12.1
    err = nothing
    try
        calc_drift_xy(im1, im2)
    catch err
    end
    @test err isa Exception
    @test contains(sprint(showerror, err), "same rotation")

    im1 = load_image("01514.sxm")
    im2 = load_image("01514.sxm")

end
