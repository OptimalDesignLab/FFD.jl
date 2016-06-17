#  Note:  evalVolume is tested under test_linearFFD.jl since it is needed
#         everytime a control point is perturbed to compute the new location of
#         the geometric point


facts("--- Checking B-spline Formulation ---") do

  U = [0,0,0,0,2/5,3/5,3/5,1,1,1,1]
  order = 4
  p = order - 1
  u = 0.42
  nctl = length(U) - order

  @fact nctl --> 7

  context("Checking knot span index evaluation") do
    span = findSpan(u, U, order, nctl)
    @fact span --> 5
  end

  context("Checking basis function evaluations") do
    N = zeros(order)
    span = findSpan(u, U, order, nctl)
    basisFunctions(U, order, u, span, N)
    @fact N[1] --> roughly(0.081, atol = 1e-15)
    @fact N[2] --> roughly(0.405, atol = 1e-15)
    @fact N[3] --> roughly(0.5136666666666666, atol = 1e-15)
    @fact N[4] --> roughly(0.00033333333333333153, atol = 1e-15)
  end

  context("Checking curve point evaluation") do
    P = 0:1/(nctl-1):1
    C = [0.0]
    evalCurve([u], U, order, P, C)
    @fact C[1] --> roughly(0.40555555555555556, atol = 1e-16)
  end

end  # End facts("--- Check Basis Function Calculations ---")
