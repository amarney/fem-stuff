function errorglob1d(uName, uGlobal, mesh, fem, derivOrder = 0,
  rule = fem.degree + 1)
## -----------------------------------------------------------------------------
# errorglob1d: computes the global error measured by the L^2 norm
#
#   Written by Angelo Marney in the summer of 2017.
#
#   Usage: errorGlob = errorglob1d(uName, uGlobal, mesh, fem, derivOrder,
#                                 rule)
#
#   Inputs: uName: the true solution function (or its derivative)
#           uGlobal: an array for global coefficients (for entire fe space)
#           mesh: a data structure for the mesh (need x and eleConn)
#           fem: a data structure for the finite element space (need feConn)
#           derivOrder: the derivative order (default 0)
#           rule: the quadrature rule to use (default fem.degree + 1)
#
#   Outputs: errorGlob: the global error measured by the L^2 norm
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if !isdefined(mesh, :x) || !isdefined(mesh, :eleConn)
    error("errorglob1d, mesh.x and mesh.eleConn must exist")
  end
  if !isdefined(fem, :feConn)
    error("errorglob1d, fem.feConn must exist")
  end
  if derivOrder < 0
    error("errorglob1d, derivative order can't be negative")
  end
  if length(uGlobal) != length(fem.x)
        error("errorglob1d, number of global coefficients must equal the total
        number of shape functions")
  end
  if rule < 0
    error("errorglob1d, rule must be a positive integer")
  end

  # Accumulate all the local errors
  errorGlob = 0
  for k = 1:size(mesh.eleConn,1) # standard way to loop over elements
    element = mesh.x[mesh.eleConn[k, :]]
    uLocal = uGlobal[fem.feConn[k, :]]
    errorLoc = errorloc1d(uName, uLocal, element, degree, derivOrder, rule)
    errorGlob = errorGlob + errorLoc
  end
  return sqrt(errorGlob)
end
