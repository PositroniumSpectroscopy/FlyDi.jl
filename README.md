# FlyPs.jl
fly dipoles through electric fields

## Install

Clone,

    $ git clone "https://github.com/PositroniumSpectroscopy/FlyPs.jl"

And import using the path to the repo, e.g.,

    julia> push!(LOAD_PATH, "/home/adam/Git/FlyPs.jl")
    julia> using FlyPs
    
Prerequisite packages:

    HDF5
    Statistics
    DataFrames
    ProgressMeter
    Dates
    SpecialFunctions
    Random, Distributions
    (Plotting package is PyPlot)

To install these you neeed to type into Julia terminal:

    julia> using Pkg
    julia> Pkg.add("PACKAGE NAME")
    
    

#Version History
0.1 (29/11/19) - First basic version but with many bugs still to fix:

            -t0 of laser selected distribution is skewed to left for some reason
            -Laser selected vx bounds correct but shape not identical to python version
            -Laser selected y not selected enough
            -fly_last/traj_last not tested
            -fly_Para not tested
            -importing of h5 files still not clear
            -Saving as feather not working for array of dataframes that flydi ouputs
