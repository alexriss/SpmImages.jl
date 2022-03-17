# SpmImages.jl

[![Build Status](https://github.com/alexriss/SpmImages.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alexriss/SpmImages.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/alexriss/SpmImages.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/alexriss/SpmImages.jl)

A julia library to read and display SPM (scanning probe microscopy) images. Currently, only [Nanonis](https://www.specs-group.com/nanonis/products/) files are supported.

The commands are somewhat similar to the ones in the respective python version [imag*ex*](https://github.com/alexriss/imagex).

## Installation

Install from the Julia package registry:

```julia
pkg> add SpmImages
```

## Getting started

Code examples can be found in the [example notebook](demo/example.ipynb).


## Related projects

- [SpmImageTycoon.jl](https://github.com/alexriss/SpmImageTycoon.jl): App to organize SPM images and spectra.
- [SpmSpectroscopy.jl](https://github.com/alexriss/SpmSpectroscopy.jl): Julia library to read and analyze SPM spectra.
- [SpmGrids.jl](https://github.com/alexriss/SpmGrids.jl): Julia library to read and analyze SPM grid spectroscopy.
- [imag*ex*](https://github.com/alexriss/imagex): Python scripts to analyze scanning probe images.
- [grid*ex*](https://github.com/alexriss/gridex): Python scripts to analyze 3D grid data.
