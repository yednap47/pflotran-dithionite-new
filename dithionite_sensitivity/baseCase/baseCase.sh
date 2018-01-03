#!/bin/bash

mpirun -np 8 ../../src/pflotran/pflotran -pflotranin 1d-allReactions-10m-uniformVelocity.in
mpirun -np 8 ../../src/pflotran/pflotran -pflotranin 1d-allReactions-10m-uniformVelocity-old.in

julia plot.jl
julia debug_oldvnew.jl
