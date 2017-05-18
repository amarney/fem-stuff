# Create a composite type (aka a c-struct) for the mesh:
type Mesh
  x::Array{Float64,1}
  eleConn::Array{Int64,2}
  iUnknown::Array{Int64,1}
  iControl::Array{Int64,1}
end

function mesh1d(domain, nEle)
## -----------------------------------------------------------------------------
# mesh1d: makes an evenly-spaced 1d mesh
#
#   Written by Angelo Marney in the summer of 2017.
#
#   Usage: mesh = mesh1d(domain, nEle)
#
#   Inputs: domain: a two-element vector containing the endpoints domain
#           nEle: the number of elements of the mesh
#
#   Outputs: mesh: a data structure containing four fields
#             mesh.x: a vector containing the node coordinates
#             mesh.eleConn: an element connectivity matrix
#             mesh.iUnknown: a vector containing unknown node indices
#             mesh.iControl: a vector containing control nodes indices
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if length(domain) != 2
    error("mesh1d, domain must have exactly two endpoints")
  end
  if nEle < 1 || length(nEle) != 1
    error("mesh1d, number of elements must be a positive integer")
  end

  # Make the element connectivity matrix
  eleConn = Array(Int64, nEle, 2)
  for k = 1:nEle
    eleConn[k,:] = [k, k + 1]
  end

  # Make the coordinate vector
  x = Array(Float64, nEle + 1, 1)
  x = linspace(domain[1], domain[2], nEle + 1)

  # Make the index vectors
  iUnknown = Array(Int64,1)
  iUnknown = [2:nEle;]

  iControl = Array(Int64,1)
  iControl = [1; nEle + 1]

  # Return the mesh
  return Mesh(x, eleConn, iUnknown, iControl)
end
