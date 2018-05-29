using FFD
using Base.Test

function not_isapprox(args...; kwargs...)
  return !isapprox(args...; kwargs...)
end


# write your own tests here
@test 1 == 1

include("test_new.jl")
