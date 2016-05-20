# runtest.jl
include("../src/b-splines.jl")

# Testing functions based on Ex2.3 from the NURBS book. considering only one
# value.
U = [0,0,0,1,2,3,4,4,5,5,5]
u = [2.5]
order = 3
degree = order - 1
n_control_pts = 7 # Figure 2.6 in the NURBS book
P = 0:5/6:5
C = [0.0]




evalCurve(u, U, order, P, C)
# Analytical test
C_test = 0.125*(P[2] + 6*P[3] + P[4])

println(" C = $C")
println("C_test = $C_test")
