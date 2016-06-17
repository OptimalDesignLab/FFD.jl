# runtest.jl
using ArrayViews
using FactCheck

# Source file includes
include("../src/mapping.jl")

#Tests
include("./test_bsplines.jl")
include("./test_linearFFD.jl")

#=
# Testing functions based on Ex2.3 from the NURBS book. considering only one
# value.
# U = [0,0,0,0, 0.1,0.2,0.5,1.0,1.0,1.0,1.0]
U = [0,0,0,0,2/5,3/5,3/5,1,1,1,1]
println("U = $U")

order = 4
degree = order - 1
n_control_pts = length(U) - order
P = 0:1/(n_control_pts-1):1

println("n_control_pts = $n_control_pts")
println("P = $P\nlength(P) = ", length(P))
# Analytical test
# C_0 = ((order)/U[order+1]) *  (P[2] - P[1])
C_1 = order*(P[n_control_pts] - P[n_control_pts-1])/(1-U[length(U)-order])

u = 0.3

jth_deriv = 1
bvalue = derivValue(U, order, u, P, n_control_pts, jth_deriv)

C = zeros(Float64,length(u))

# evalCurve(u, U, order, P, C)
@printf("P = [")
for i = 1:length(P)
  @printf("%f, ", P[i])
end
println("]\nC = $C")




println("bvalue = $bvalue")
println("C_1 = $C_1")


include("../src/knot.jl")
order = 3
nctl = 7
X = 0:0.1:1
U = zeros(order+nctl)
calcKnot(order, nctl, X, U)
println("U = $U")
=#
