# control_point.jl
@doc """
### controlPoint

Genertates a uniformly distributed control points in the physical space.
Intended only for initializing of Mapping object as of this commit.

**Inputs**

*  `map`  : Object of Mapping type
*  `box`  : Object of BoundingBox type

"""->

function controlPoint(map::AbstractMappingType, box::AbstractBoundingBox)

  # Remember B-spline control points also lie on the interior of the volume
  nctl =sview(map.nctl, :)
  origin =sview(box.origin, :)
  S =sview(box.unitVector, :, 1)*box.lengths[1]
  T =sview(box.unitVector, :, 2)*box.lengths[2]
  U =sview(box.unitVector, :, 3)*box.lengths[3]

  for k = 0:map.nctl[3]-1
    for j = 0:map.nctl[2]-1
      for i = 0:map.nctl[1]-1
        map._cp_xyz[:,i+1,j+1,k+1] = box.origin + (i/(nctl[1]-1))*S +
                                    (j/(nctl[2]-1))*T + (k/(nctl[3]-1))*U
      end
    end
  end

  return nothing
end  # End function controlPoint

@doc """
### writeControlPointsVTS

Writes a `*.vts` file to be viewed in Paraview. Necessary for visualizing the
control points for testing

**Input**

* `map` : Object of AbstractMappingType

"""->

function writeControlPointsVTS(map::AbstractMappingType,
                               fname::AbstractString="control_points")

  vtsfile = vtk_grid(fname, real(map.cp_xyz))
  outfiles = vtk_save(vtsfile)

  return nothing
end
