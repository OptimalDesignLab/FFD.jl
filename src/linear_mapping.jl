@doc """
### Linear mapping

Create a linear mapping between physical and parametric coordinate systems.
*Physical coordinate* system refers to the domian in which the actual geometry/
mesh exists. Its axes are denoted by (X,Y,Z).
*Parametric coordinate* systems refers to the local coordinate system for the
bounding box comprising of the control points. Its axes are denoted by (s,t,u)

The objective of linear mapping is to generate (s,t,u) coordinates of points in
the (x,y,z) space. Presently, Only a bounding box which is a cube varying
between [0,1] in each of the s,t,u direction.

"""->
function linearMap(map, box)

  return nothing
end
