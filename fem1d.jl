# Create a composite type for the finite element space:
type Fem
  x::Array{Float64,1}
  feConn::Array{Int64,2}
  iUnknown::Array{Int64,1}
  iControl::Array{Int64,1}
  degree::Int64
end

function fem1d(mesh::Mesh, degree::Int64)
## -----------------------------------------------------------------------------
# fem1d: makes a 1d finite element space data structure
#
#   Written by Angelo Marney in the summer of 2017.
#
#   Usage: fem = fem1d(mesh, degree)
#
#   Inputs: mesh: a data structure for the mesh (need x and eleConn)
#           degree: the degree (aka order) of the finite element space
#
#   Outputs: fem: a data structure containing five fields
#             fem.x: a vector containing the FE space coordinates
#             fem.feConn: an FE point connectivity matrix
#             fem.iUnknown: a vector containing unknown node indices
#             fem.iControl: a vector containing control nodes indices
#             fem.degree: the degree of the finite element space
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if degree < 1 || length(degree) != 1
    error("fem1d, finite element space must be of degree 1 or greater")
  end
  if !isdefined(mesh, :x) || !isdefined(mesh, :eleConn)
    error("fem1d, mesh.x and mesh.eleConn must exist")
  end

  # Create the finite element point connectivity matrix
  feConn = Array(Int64, size(mesh.eleConn, 1), degree + 1)
  feConn[1, :] = [1:(degree + 1);].'
  for k = 2:size(mesh.eleConn, 1)
    feConn[k, :] = feConn[k - 1, :] + degree
  end

  # Create the coordinate vector
  numPoint = feConn[end,end]
  x = Array(Float64, (degree + 1) + (size(mesh.eleConn, 1) - 1)*degree, 1)
  x = linspace(mesh.x[1], mesh.x[end], numPoint)

  # Create the index arrays
  iUnknown = [2:numPoint-1;]
  iControl = [1; numPoint]

  # Return the fem data structure
  return Fem(x, feConn, iUnknown, iControl, degree)
end
