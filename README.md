# FlyPs.jl
fly dipoles through electric fields

## Install

Clone,

    $ git clone "https://github.com/PositroniumSpectroscopy/FlyPs.jl"

And import using the path to the repo, e.g.,

    julia> push!(LOAD_PATH, "/home/adam/Git/FlyDi.jl")
    julia> using FlyDi
    
Prerequisite packages:
HDF5
Statistics
DataFrames
Plotting package is PyPlot

To install these you neeed to type into Julia terminal:
    julia> using Pkg
    julia> Pkg.add("PACKAGE NAME")
