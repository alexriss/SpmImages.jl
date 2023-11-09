using SpmImages
using Test

@testset "File loading sxm" begin
    ima = load_image("Image_445.sxm")
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
