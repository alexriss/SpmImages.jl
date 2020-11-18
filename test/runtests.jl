using SpmImages
using Test

@testset "SpmImages.jl" begin
    ima = load_image("Image_445.sxm")
    @test ima.channel_names == ["Z", "LI_X", "LI_Y", "Bias", "Current", "Phase", "Amplitude", "Frequency_Shift", "Excitation"]
    @test ima.channel_units == ["m", "V", "V", "V", "A", "deg", "m", "Hz", "V"]
    @test ima.scansize == [2.0, 4.0]
    @test ima.scansize_unit == "nm"
    @test ima.pixelsize == [96, 192]
    @test ima.scan_direction == SpmImages.up
    @test ima.acquisition_time == 4877.1
    @test ima.header["oscillation_control>reference_phase_(deg)"] == "65.16E+0"
    @test ima.data[12,42,3] == -3.6583905f0
    @test ima.data[9,18,10] == 3.0707468f-11
    @test ima.data[69,69,6] == -0.7235545f0
end
