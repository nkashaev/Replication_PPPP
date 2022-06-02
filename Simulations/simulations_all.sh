#!/bin/bash


julia Simulations/simulation_main.jl Epanechnikov 500 4
julia Simulations/simulation_main.jl Epanechnikov 1000 4
julia Simulations/simulation_main.jl Epanechnikov 1500 4

