function biformglob1d(kerName, mesh, ufem, vfem, derivOrders = [0; 0],
  rule = max(ufem.degree, vfem.degree) + 1)
## -----------------------------------------------------------------------------
# biformglob1d: computes global "stiffness matrix" from \int{kerName*phi*test}
#
#   Written by Angelo Marney in the summer of 2017.
#
#   Usage: matGlob = biformglob1d(kerName, mesh, ufem, vfem, derivOrders, rule)
#
#   Inputs: kerName: the "coefficient" function in the bilinear form
#           mesh: a data structure for the mesh (need x and eleConn)
#           ufem: a data structure for the u-fe space (need feConn)
#           vfem: a data structure for the v-fe space (need feConn)
#           derivOrders: a two-element array of derivative orders (default 0)
#           rule: the quadrature rule to use (default degree + 1)
#
#   Outputs: matGlob: global "stiffness matrix" associated with a bilinear form
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if (!isdefined(mesh, :x) || !isdefined(mesh, :eleConn))
    error("biformglob1d, mesh.x and mesh.eleConn must exist")
  end
  if ((ufem.degree < 1) || (vfem.degree < 1))
    error("biformglob1d, finite element spaces can't have degree less
    than one")
  end
  if ((derivOrders[1] < 0) || (derivOrders[2] < 0))
    error("biformglob1d, derivative orders can't be negative")
  end
  if rule < 0
    error("biformglob1d, rule must be a positive integer")
  end

  # Set the degrees array
  degrees = [ufem.degree, vfem.degree]

  # Loop over each element, and accumulate the local stiffness mattrices
  matGlob = spzeros(length(ufem.x), length(vfem.x))
  for k = 1:size(mesh.eleConn, 1)
    element = mesh.x[mesh.eleConn[k, :]]
    matLoc = biformloc1d(kerName, element, degrees, derivOrders, rule)
    matGlob[ufem.feConn[k, :], vfem.feConn[k, :]] = matGlob[ufem.feConn[k, :],
    vfem.feConn[k, :]] + matLoc
  end
  # Return the global "stiffness matrix"
  return matGlob
end
