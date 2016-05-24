# plotting.jl
using Gadfly
using Cairo
using Fontconfig

include("./b-splines.jl")

function evaluateBasis(U, u, order, nctl, concerned_index, N)

  active_bases = zeros(Float64, order)

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

  return nothing
end


U = [0,0,0,0, 0.1,0.2,0.5,1.0,1.0,1.0,1.0]
u = 0:0.01:1
order = 4
degree = order - 1
nctl = length(U) - order
nbases = length(U) - order # Total number of basis functions for a knot vector
concerned_index = 3 # Index of the basis function to be plotted
N = zeros(Float64, length(u), nbases)
for i = 1:size(N,2)
  evaluateBasis(U, u, order, nctl, concerned_index, N[:,i])
end

# println("N = $N")

myplot = plot(
  layer(x=u, y=N[1], Geom.line, Theme(default_color=colorant"green")),
  layer(x=u, y=N[2], Geom.line, Theme(default_color=colorant"blue")),
  layer(x=u, y=N[3], Geom.line, Theme(default_color=colorant"red")),
  layer(x=u, y=N[4], Geom.line, Theme(default_color=colorant"black")),
  layer(x=u, y=N[5], Geom.line, Theme(default_color=colorant"brown")),
  layer(x=u, y=N[6], Geom.line, Theme(default_color=colorant"pink")),
  layer(x=u, y=N[7], Geom.line, Theme(default_color=colorant"yellow")),
  Guide.XLabel("u"),
  Guide.YLabel("Basis Function Values"),
  Guide.yticks(ticks=collect(0:0.2:1)),
)

draw(PDF("Bases.pdf", 10inch, 10inch), myplot)
