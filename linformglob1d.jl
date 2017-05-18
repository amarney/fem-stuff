function linformglob1d(fName, mesh, fem, derivOrder = 0, rule = fem.degree + 1)
## -----------------------------------------------------------------------------
# linformglob1d: computes global "load vector" from \int { fName*test }
#
#   Written by Angelo Marney in the summer of 2017.
#
#   Usage: vecGlob = linformglob1d(fName, mesh, fem, derivOrder, rule)
#
#   Inputs: fName: the function in the linear form
#           mesh: a data structure for the mesh (need x and eleConn)
#           fem: a data structure for the finite element space (need feConn)
#           derivOrder: the derivative order (default 0)
#           rule: the quadrature rule to use (default fem.degree + 1)
#
#   Outputs: vecGlob: the global "load vector" associated with a linear form
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if !isdefined(mesh, :x) || !isdefined(mesh, :eleConn)
    error("linformglob1d, mesh.x and mesh.eleConn must exist")
  end
  if !isdefined(fem, :feConn)
    error("linformglob1d, fem.feConn must exist")
  end
  if derivOrder < 0
    error("linformglob1d, derivative order can't be negative")
  end
  if rule < 0
    error("linformglob1d, rule must be a positive integer")
  end

  # Loop over each element, and accumulate the local load vectors
  degree = fem.degree
  vecGlob = zeros(length(fem.x), 1)
  for k = 1:size(mesh.eleConn, 1)
    element = mesh.x[mesh.eleConn[k, :]]
    vecLoc = linformloc1d(fName, element, degree, derivOrder, rule)
    vecGlob[fem.feConn[k, :]] = vecGlob[fem.feConn[k, :]] + vecLoc
  end

  # Return the local "load vector"
  return vecGlob
end
