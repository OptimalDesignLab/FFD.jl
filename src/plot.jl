# plotting.jl
include("./b-splines.jl")

U = [0,0,0,0, 0.1,0.2,0.5,1.0,1.0,1.0,1.0]
u = 0:0.1:1
order = 4
degree = order - 1
nctl = length(U) - order
active_bases = zeros(Float64, order)
concerned_index = 3 # Index of the basis function to be plotted
N = zeros(Float64, length(u))

for i = 1:length(u)
  span = findSpan(u[i], nctl, U, order)
  basisFunctions(U, order, u[i], span, active_bases)
  for j = 1:length(active_bases)
    index = span - order + j
    if index == concerned_index
      N[i] = active_bases[j]
    end # end if index == concerned_index
  end # end for j = 1:length(active_bases)
end  # end for i = 1:length(u)

println("N = $N")
