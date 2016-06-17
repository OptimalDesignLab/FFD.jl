# control_point.jl
@doc """
### controlPoint

Genertates a uniformly distributed control points in the physical space.
Intended only for initializing of Mapping object as of this commit.

**Inputs**

*  `map`  : Object of Mapping type
*  `box`  : Object of BoundingBox type

"""->

function controlPoint(map, box)

  # Remember B-spline control points also lie on the interior of the volume
  nctl = view(map.nctl, :)
  origin = view(box.origin, :)
  S = view(box.unitVector, :, 1)*box.lengths[1]
  T = view(box.unitVector, :, 2)*box.lengths[2]
  U = view(box.unitVector, :, 3)*box.lengths[3]

  for k = 0:map.nctl[3]-1
    for j = 0:map.nctl[2]-1
      for i = 0:map.nctl[1]-1
        map.cp_xyz[i+1,j+1,k+1,:] = box.origin + (i/(nctl[1]-1))*S +
                                    (j/(nctl[2]-1))*T + (k/(nctl[3]-1))*U
      end
    end
  end

  return nothing
end  # End function controlPoint
