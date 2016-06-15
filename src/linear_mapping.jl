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
function void()
  return nothing
end

@doc """
### linearMap

Take in one coordinate in (x,y,z) space and convert it to parametric (s,t,u)
space

**Arguments**

*  `map` : Object of mapping type
*  `box` : BoundingBox object
*  `X`   : Point coordinate in (x,y,z) space
*  `pX`   : corresponding coordinates in (s,t,u) space

"""->

function linearMap(map, box, X, pX)

  # The assumption of the mapping for the bounding box presently is that
  # coordinate transformation from physical x,y,z coordinate transformation to
  # parametric s,t,u coordinate involved only translation and scaling. There is
  # no rotation.

  # Get the x,y,z coordinates for the origin of the s,t,u system
  origin = box.origin
  S = box.unitVector[:,1]
  T = box.unitVector[:,2]
  U = box.unitVector[:,3]
  # lengthS = box.lengths[1]
  # lengthT = box.lengths[2]
  # lengthU = box.lengths[3]

  XmX0 = X - origin

  # calculate s
  # Division by lengths[i] to ensure s,t,u lie between [0,1]
  TcrossU = cross(T,U)
  s = dot(TcrossU,XmX0)/(dot(TcrossU,S)*box.lengths[1])

  # Calculate t
  ScrossU = cross(S,U)
  t = dot(ScrossU,XmX0)/(dot(ScrossU,T)*box.lengths[2])

  # calculate u
  ScrossT = cross(S,T)
  u = dot(ScrossT,XmX0)/(dot(ScrossT,U)*box.lengths[3])

  pX[:] = [s,t,u]

  return nothing
end

function calcParametricLinear()

  retrun nothing
end  # End function calcParametricLinear
