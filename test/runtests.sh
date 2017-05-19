#!/bin/sh

julia runtests_serial.jl
mpirun -np 2 julia runtests_parallel.jl
