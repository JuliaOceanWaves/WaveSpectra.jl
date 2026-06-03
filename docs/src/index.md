# WaveSpectra.jl

## Overview

WaveSpectra.jl provides a unit-aware tool for working with ocean wave spectra.
This package enables input, conversion, creation, and characterization
for cartesian and polar spectra, omnidirectional spectra, and some parametric spectra.

## Package Guide

The [Package Guide](@ref Quickstart) is the main documentation for the package and includes
usage details and examples of _WaveSpectra.jl's_ capabilities. Including:

- [Directional wave spectra](@ref dir_spectra)
    - In cartesian coordinates
    - In polar coordinates
- [Omnidirectional wave spectra](@ref omnidir_spectra)
- [Functions for characterization](@ref spectra_funcs)
- [Parametric wave spectra](@ref param_spectra)
- and [other auxilary functions](@ref other_funcs)

## [Supported Spectral-Variables](@id spectral_var_cube)

As part of the unit-aware capabilities, the tool only accepts the following 
spectral-variables types: temporal/spatial, frequency/period, and linear/angular 
combinations. Represented below:

![dispersion_cube_sample](./assets/Commutative_diagram_of_harmonic_wave_properties.svg)

## Reference

See the [API reference](@ref API) for all exported functions.
