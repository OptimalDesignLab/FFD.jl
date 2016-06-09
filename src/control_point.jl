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

  # Get control point locations for the edges
  # In X-direction
  nctl = map.nctl[1]
  low = box.box_coord[1,1]
  high = box.box_coord[2,1]
  step = (high-low) / (nctl-1)
  for i = 1:map.nctl[1]
    map.cp_xyz[i,:,:,1] = low
    low += step
  end

  # In the Y-direction
  if box.ndim >= 2
    nctl = map.nctl[2]
    low = box.box_coord[1,2]
    high = box.box_coord[2,2]
    step = (high-low) / (nctl-1)
    for i = 1:map.nctl[2]
      map.cp_xyz[:,i,:,1] = low
      low += step
    end
  end

  # In the Z-direction
  if box.ndim == 3
    nctl = map.nctl[3]
    low = box.box_coord[1,3]
    high = box.box_coord[2,3]
    step = (high-low) / (nctl-1)
    for i = 1:map.nctl[3]
      map.cp_xyz[:,:,i,1] = low
      low += step
    end
  end

  return nothing
end
