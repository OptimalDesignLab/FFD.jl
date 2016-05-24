# runtest.jl
include("../src/b-splines.jl")

# Testing functions based on Ex2.3 from the NURBS book. considering only one
# value.
# U = [0,0,0,0, 0.1,0.2,0.5,1.0,1.0,1.0,1.0]
U = [0,0,0,1,2,3,4,4,5,5,5]
println("U = $U")
# u = [2.5]
order = 3
degree = order - 1
#n_control_pts = 7 # Figure 2.6 in the NURBS book
# P = [1,2,3,4,5]
n_control_pts = length(U) - order
println("n_control_pts = $n_control_pts")
P = 0:5/(n_control_pts-1):5
println("P = $P\nlength(P) = ", length(P))
# Analytical test
C_test = 0.125*(P[3] + 6*P[4] + P[5])


u = 0:5
# println("val = ", n_control_pts + order - length(U))

C = zeros(Float64,length(u))

evalCurve(u, U, order, P, C)
@printf("P = [")
for i = 1:length(P)
  @printf("%f, ", P[i])
end
println("]\nC = $C")
# println("C analytical = $C_test")
