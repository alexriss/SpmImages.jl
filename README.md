<img width="100" height="100" src="logo.svg?raw=true" />

# SpmImages.jl

[![Build Status](https://github.com/alexriss/SpmImages.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alexriss/SpmImages.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/alexriss/SpmImages.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/alexriss/SpmImages.jl)
<a href="https://twitter.com/00alexx"><img src="https://img.shields.io/twitter/follow/00alexx?style=social" alt="Twitter"></a>

A julia library to read and display SPM (scanning probe microscopy) images. Currently, only [Nanonis](https://www.specs-group.com/nanonis/products/) files are supported.

The commands are somewhat similar to the ones in the respective python version [imag*ex*](https://github.com/alexriss/imagex).

## Installation

Install from the Julia package registry:

```julia
pkg> add SpmImages
```

## Getting started

```julia
using SpmImages
using Plots  # needs to be manually imported for plotting functions

# load an SPM image
ima = load_image("Image_445.sxm");

# Plot a specific image channel directly from the image
plot_channel(ima, "frequency shift")

# Get a channel as 2D data and directly plot it
c = get_channel(ima, "current");
plot_data(c.data, color=:lapaz, title="current", legend=:none)

# Plot using the Images library
using Images
d = get_channel(ima, "current bwd", origin="upper").data
d = d .- minimum(d)
d = d ./ maximum(d)
Gray.(d)
```

Code examples can be found in the [example notebook](demo/example.ipynb).

## Get in touch

Please post issues, suggestions, and pull requests on github. <a href="https://twitter.com/00alexx">Follow me on twitter</a> for updates and more information about this project: 
<a href="https://twitter.com/00alexx"><img src="https://img.shields.io/twitter/follow/00alexx?style=social" alt="Twitter"></a>

## Related projects

- [SpmImageTycoon.jl](https://github.com/alexriss/SpmImageTycoon.jl): App to organize SPM images and spectra.
- [SpmSpectroscopy.jl](https://github.com/alexriss/SpmSpectroscopy.jl): Julia library to read and analyze SPM spectra.
- [SpmGrids.jl](https://github.com/alexriss/SpmGrids.jl): Julia library to read and analyze SPM grid spectroscopy.
- [imag*ex*](https://github.com/alexriss/imagex): Python scripts to analyze scanning probe images.
- [grid*ex*](https://github.com/alexriss/gridex): Python scripts to analyze 3D grid data.
