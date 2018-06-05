#  Note:  evalVolume is tested under test_linearFFD.jl since it is needed
#         everytime a control point is perturbed to compute the new location of
#         the geometric point


@testset "--- Checking B-spline Formulation ---" begin

  U = [0,0,0,0,2/5,3/5,3/5,1,1,1,1]
  order = 4
  p = order - 1
  u = 0.42
  nctl = length(U) - order

  @test ( nctl )== 7

  @testset "Checking knot span index evaluation" begin
    span = FFD.findSpan(u, U, order, nctl)
    @test ( span )== 5
  end

  @testset "Checking basis function evaluations" begin
    N = zeros(order)
    span = FFD.findSpan(u, U, order, nctl)
    FFD.basisFunctions(U, order, u, span, N)
    @test isapprox( N[1], 0.081) atol= 1e-15
    @test isapprox( N[2], 0.405) atol= 1e-15
    @test isapprox( N[3], 0.5136666666666666) atol= 1e-15
    @test isapprox( N[4], 0.00033333333333333153) atol= 1e-15
  end

  @testset "Checking curve point evaluation" begin
    P = 0:1/(nctl-1):1
    C = [0.0]
    FFD.evalCurve([u], U, order, P, C)
    @test isapprox( C[1], 0.40555555555555556) atol= 1e-16
  end

end  # End facts("--- Check Basis Function Calculations ---")

@testset "--- Checking B-spline Derivatives ---" begin

  @testset "Checking Basis Function Derivatives" begin

    U = [0,0,0,0,2/5,3/5,3/5,1,1,1,1]
    order = 4
    p = order - 1
    u = 0.42
    nctl = length(U) - order
    N = zeros(AbstractFloat, order)
    Nderiv = zeros(AbstractFloat, order)
    span = FFD.findSpan(u, U, order, nctl)
    FFD.derivBasisFunctions(u, U, order, span, N, Nderiv)

    @test isapprox( N[1], 0.081) atol= 1e-15
    @test isapprox( N[2], 0.405) atol= 1e-15
    @test isapprox( N[3], 0.5136666666666666) atol= 1e-15
    @test isapprox( N[4], 0.00033333333333333153) atol= 1e-15
    @test isapprox( Nderiv[1], -1.3500000000000003) atol= 1e-15
    @test isapprox( Nderiv[2], -2.25) atol= 1e-15
    @test isapprox( Nderiv[3], 3.55) atol= 1e-15
    @test isapprox( Nderiv[4], 0.04999999999999982) atol= 1e-15

  end  # End context("Checking Basis Function Derivatives")

end  # End facts("--- Checking B-spline Derivatives ---")
