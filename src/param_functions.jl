# Functions for parametric representations
@doc """
### para3D

Calculates the parametric locations for a 3D block Taken from GAETAN's code

**Inputs**

*  `X` : Coordinates, size(l,m,n,ndim)
*  `S` : The u, v, w parametric positions

**Outputs**

*  `u` : Averaged u parameters size(l)
*  `v` : Averaged v parameters size(m)
*  `w` : Averaged w parameters size(n)

REFERENCE: Gaetan parameterizations.f90

"""->

function para3D(X, S)

  # Zeros 3 low end faces (or edges if one plane is specified)
  S[1,:,:,1] = 0.0
  S[:,1,:,2] = 0.0
  S[:,:,1,3] = 0.0

  # Zeros out before working
  u = 0.0
  v = 0.0
  w = 0.0

  # Set up the low-edge lines since they are missed by the following loops over
  # most of the low end faces
  for i = 2:l
    S[i,1,1,1] = S[i-1,1,1,1] + del_i(i,1,1)
  end

  for j = 2:m
    S[1,j,1,2] = S[1,j-1,1,2] + del_j(1,j,1)
  end

  for k = 2:n
    S[1,1,k,3] = S[1,1,k-1,3] + del_k(1,1,k)
  end

  # Set up the rest of the low end face lines beacuse theu are missed by the
  # main loop over most of the volume
  for k = 2:n
    for j = 2:m
      S[1,j,k,2] = S[1,j-1,k,2] + del_j(1,j,k)
      S[1,j,k,3] = S[1,j,k-1,3] + del_k(1,j,k)
    end
  end

  for i = 2:l
    for k = 2:n
      S[i,1,k,1] = S[i-1,1,k,1] + del_i(i,1,k)
      S[i,1,k,3] = S[i,1,k-1,3] + del_k(i,1,k)
    end
  end

  for i = 2:l
    for j = 2:m
      S[i,j,1,1] = S[i-1,j,1,1] + del_i(i,j,1)
      S[i,j,1,2] = S[i,j-1,1,2] + del_j(i,j,1)
    end
  end

  # Traverse the block once for all lines except those within the low-end faces
  for i = 2:l
    for j = 2:m
      for k = 2:l
        S[i,j,k,1] = S[i-1,j,k,1] + del_i(i,j,k)
        S[i,j,k,2] = S[i,j-1,k,2] + del_j(i,j,k)
        S[i,j,k,3] = S[i,j,k-1,3] + del_k(i,j,k)
      end
    end
  end

  # Normalizing requires another pass through the volume.
  # Handle lines of zero length first by inserting uniform distributions. Then
  # the standard normalizations can be applied safely everywhere.

  for j = 1:m
    for k = 1:n
      if S[l,j,k,1] == 0.0 # zero length lines in `i` direction
        for i = 2:l
          S[i,j,k,1] = i - 1
        end
      end
    end
  end

  for i = 1:l
    for k = 1:n
      if S[i,m,k,2] == 0.0 # zero length lines in `j` direction
        for j = 2:m
          S[i,j,k,2] = j-1
        end
      end
    end
  end

  for i = 1:l
    for j = 1:m
      if S[i,j,l,3] == 0.0 # zero length lines in `k` direction
        for k = 2:l
          S[i,j,k,3] = k-1
        end
      end
    end
  end

  # Normalize
  for i = 1:l
    for j = 1:m
      for k = 1:n
        S[i,j,k,1] = S[i,j,k,1] / S[n,j,k,1]
        S[i,j,k,2] = S[i,j,k,2] / S[i,m,k,2]
        S[i,j,k,3] = S[i,j,k,3] / S[i,j,n,3]
      end
    end
  end

  # Finally, precise 1s for the three high end faces
  for j = 1:m
    for k = 1:n
      S[l,j,k,1] = 1.0
    end
  end

  for i = 1:l
    for k = 1:n
      S[i,m,k,2] = 1.0
    end
  end

  for i = 1:l
    for j = 1:m
      S[i,j,l,3] = 1.0
    end
  end

  # Get an average u, v, w
  # average u
  for j = 1:m
    for k = 1:n
      u += S[:,j,k,1]
    end
  end
  u = u/(m*n)

  # Average v
  for i = 1:l
    for k = 1:n
      v += S[i,:,k,2]
    end
  end
  v = v/(l*n)

  # Average w
  for i = 1:l
    for j = 1:m
      w += S[i,j,:,3]
    end
  end
  w = w/(l*m)

  return u, v, w
end

@doc """
### del_i, del_j, del_k

Functions to compute the distance between 2 points along a particular parametric
direction.

**Inputs**

*  `X` : Coordinates
*  `i`, `j`, `k` : Indices in the 3 directions

** Outputs**

*  `deli`, `delj`, `delk` : distances along the 3 directions
"""->

function del_i(X, i, j, k)

  deli = sqrt( (X[i,j,k,1] - X[i-1,j,k,1])^2 + (X[i,j,k,2] - X[i-1,j,k,2])^2 +
               (X[i,j,k,3] - X[i-1,j,k,3])^2 )

  return deli
end

function del_j(X, i, j, k)

  delj = sqrt( (X[i,j,k,1] - X[i,j-1,k,1])^2 + (X[i,j,k,2] - X[i,j-1,k,2])^2 +
               (X[i,j,k,3] - X[i,j-1,k,3])^2 )

  return delj
end

function del_k(X, i, j, k)

  delk = sqrt( (X[i,j,k,1] - X[i,j,k-1,1])^2 + (X[i,j,k,2] - X[i,j,k-1,2])^2 +
               (X[i,j,k,3] - X[i,j,k-1,3])^2 )

  return delk
end
