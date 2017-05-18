function boundaryclassify1d(fem, domain, indexBC)
## -----------------------------------------------------------------------------
# boundaryclassify1d: creates data structures to aid in solving BVPs
#
#   Written by Angelo Marney in the summer of 2017.
#
#   Usage: dof, nt = boundaryclassify1d(fem, domain, indexBC)
#
#   Inputs: fem: a data structure for the finite element space (need feConn)
#           domain: a two-element vector containing the endpoints domain
#           indexBC: a two-element vector classifying boundary conditions
#
#   Outputs: dof: array of fem.x with essential boundary points skipped
#            nt: array of fem.x with all points classified as -1, 0, or -200
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if length(domain) != 2
    error("boundaryclassify1d, domain must have exactly two endpoints")
  end
  if (!isdefined(fem, :feConn) || !isdefined(fem, :x))
    error("boundaryclassify1d, fem.x and fem.feConn must exist")
  end
  if (indexBC[1] != 0 && indexBC[1] != -1) ||
      (indexBC[2] != 0  && indexBC[2] != -1)
    error("boundaryclassify1d, indexBC must contain 0s or -1s")
  end

  # indexBC[1] = index for left endpoint of domain
  # indexBC[2] = index for right endpoint of domain
  # if index = 0, the endpoint is essential/dirichlet
  # if index = -1, the endpoint is natural/neumann
  # otherwise it is an interior point and has index -200

  nt = Array(Int64, size(fem.x))
  tolerance = eps() # use Julia's tolerance
  for k = 1:length(fem.x)
    if abs(fem.x[k] - domain[1]) <= tolerance
      nt[k] = indexBC[1]
    elseif abs(fem.x[k] - domain[2]) <= tolerance
      nt[k] = indexBC[2]
    else
      nt[k] = -200
    end
  end

  #Make dofu
  dof = find(nt .< 0)

  #return dofu and nt
  return dof, nt
end
