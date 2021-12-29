#!/bin/bash


julia Simulations/simulation_main.jl Epanechnikov 500 4 false
julia Simulations/simulation_main.jl Epanechnikov 1000 4 false
julia Simulations/simulation_main.jl Epanechnikov 1500 4 false
julia Simulations/simulation_main.jl Epanechnikov 5000 4 false

julia Simulations/simulation_main.jl Biweight 500 4 false
julia Simulations/simulation_main.jl Biweight 1000 4 false
julia Simulations/simulation_main.jl Biweight 1500 4 false
julia Simulations/simulation_main.jl Biweight 5000 4 false

julia Simulations/simulation_main.jl Triweight 500 4 false
julia Simulations/simulation_main.jl Triweight 1000 4 false
julia Simulations/simulation_main.jl Triweight 1500 4 false
julia Simulations/simulation_main.jl Triweight 5000 4 false

julia Simulations/simulation_main.jl Epanechnikov 500 4 true
julia Simulations/simulation_main.jl Epanechnikov 1000 4 true
julia Simulations/simulation_main.jl Epanechnikov 1500 4 true
julia Simulations/simulation_main.jl Epanechnikov 5000 4 true

julia Simulations/simulation_main.jl Biweight 500 4 true
julia Simulations/simulation_main.jl Biweight 1000 4 true
julia Simulations/simulation_main.jl Biweight 1500 4 true
julia Simulations/simulation_main.jl Biweight 5000 4 true

julia Simulations/simulation_main.jl Triweight 500 4 true
julia Simulations/simulation_main.jl Triweight 1000 4 true
julia Simulations/simulation_main.jl Triweight 1500 4 true
julia Simulations/simulation_main.jl Triweight 5000 4 true


julia Simulations/simulation_main.jl Epanechnikov 500 5 true
julia Simulations/simulation_main.jl Epanechnikov 1000 5 true
julia Simulations/simulation_main.jl Epanechnikov 1500 5 true
julia Simulations/simulation_main.jl Epanechnikov 5000 5 true

julia Simulations/simulation_main.jl Biweight 500 5 true
julia Simulations/simulation_main.jl Biweight 1000 5 true
julia Simulations/simulation_main.jl Biweight 1500 5 true
julia Simulations/simulation_main.jl Biweight 5000 5 true

julia Simulations/simulation_main.jl Triweight 500 5 true
julia Simulations/simulation_main.jl Triweight 1000 5 true
julia Simulations/simulation_main.jl Triweight 1500 5 true
julia Simulations/simulation_main.jl Triweight 5000 5 true
