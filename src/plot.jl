# plotting.jl
# plot the B-spline basis function
include("./b-splines2.jl")

#=
calls  bsplvb, interv
      integer i,j,k,left,leftmk,mflag,n,npoint
      real dx,t(10),values(7),x,xl

      data n,k /7,3/, t /3*0.,2*1.,3.,4.,3*6./
      data values /7*0./, npoint /31/

      xl = t(k)
      dx = (t(n+1)-t(k))/float(npoint-1)

      print 600,(i,i=1,5)
  600 format('1  x',8x,5('b',i1,'(x)',7x))
      do 10 i=1,npoint
         x = xl + float(i-1)*dx
         call interv ( t, n+1, x, left, mflag )
         leftmk = left - k
         call bsplvb ( t, k, 1, x, left, values(leftmk+1) )
         print 610, x, (values(j),j=3,7)
  610    format(f7.3,5f12.7)
         do 10 j=1,k
   10       values(leftmk+j) = 0.
                                        stop
      end

=#

dimension = 7
order = 3
U = [0.0, 0.0, 0.0, 1.0, 1.0, 3.0, 4.0, 6.0, 6.0, 6.0] # start and termination of knot vectors
                                                       # have a multiplicity of "order"

values = zeros(Float64, 7)
npoint = 31
# Set leftmost evaluation pont x_l and spacing dx
x_l = U[order]
dx = (U[dimension + 1] - U[dimension])/(npoint - 1)
nctl = dimension + 1

for i = 1:npoint
  x = x_l + (i-1)*dx
  span = findSpan(x, nctl, U, order-1) # findSpan(u[i], nctl, U, p)
  left_mk = span - order
  # basisFunctions(span, x, order-1, U, values) # i, u, p, U, N
  basisFunctions(U, order, 1, x, span, values)  #T, jhigh, index, x, left, N
  println("x = $(round(x,2)), N = $(values[3:7])")
  fill!(values, 0.0) # zero out for the next loop
end  # End for i = 1:npoint
