using SpmImages
using Test

@testset "File loading" begin
    ima = load_image("Image_445.sxm")
    @test ima.channel_names == ["Z", "LI X", "LI Y", "Bias", "Current", "Phase", "Amplitude", "Frequency Shift", "Excitation"]
    @test ima.channel_units == ["m", "V", "V", "V", "A", "deg", "m", "Hz", "V"]
    @test ima.scansize == [2.0, 4.0]
    @test ima.scansize_unit == "nm"
    @test ima.pixelsize == [96, 192]
    @test ima.scan_direction == SpmImages.up
    @test ima.acquisition_time == 4877.1
    @test ima.header["Oscillation Control>Reference Phase (deg)"] == "65.16E+0"
    @test ima.data[12,42,3] == -3.6583905f0
    @test ima.data[9,18,10] == 3.0707468f-11
    @test ima.data[69,69,6] == -0.7235545f0

    @test get_channel(ima, "Frequency shift").name == "Frequency Shift"
    @test get_channel(ima, "exCitation bwd").name == "Excitation"
    @test get_channel(ima, "exCitation bwd").direction == SpmImages.bwd
    @test get_channel(ima, "current bwd").data[30, 55] == 6.729436f-10
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
