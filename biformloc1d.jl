function biformloc1d(kerName, element, degrees, derivOrders = [0; 0],
  rule = maximum(degrees + 1))
## -----------------------------------------------------------------------------
# bilinformloc1d: computes local "stiffness matrix" from \int{kerName*LocShape^2}
#
#   Written by Angelo Marney in the summer of 2017.
#
#   Usage: matLoc = biformloc1d(kerName, element, degrees, derivOrders, rule)
#
#   Inputs: kerName: the "coefficient" function in the bilinear form
#           element: an array of coordinates for the element containing x
#           degrees: a two-element array of degrees for a fe-space pair
#           derivOrders: a two-element array of derivative orders (default 0)
#           rule: the quadrature rule to use (default  maximum(degrees + 1))
#
#   Outputs: matLoc: local "stiffness matrix" associated with a bilinear form
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if length(element) != 2
    error("biformloc1d, element must have exactly two endpoints")
  end
  if ((degrees[1] < 1) || (degrees[2] < 1))
    error("biformloc1d, finite element spaces can't have degree less
    than one")
  end
  if ((derivOrders[1] < 0) || (derivOrders[2] < 0))
    error("biformloc1d, derivative orders can't be negative")
  end
  if rule < 0
    error("biformloc1d, rule must be a positive integer")
  end

  # Compute the quadrature points and weights
  qPoint, qWeight = quadloc1d(element, rule)

  # Integrate over each local shape function for the given element
  matLoc = Array(Float64, degrees[1] + 1, degrees[2] + 1);
  for iShape1 = 1:(degrees[1] + 1)
    for iShape2 = 1:(degrees[2] + 1)
      integrand = kerName(qPoint).*shapeloc1d(qPoint, element, degrees[1],
      iShape1, derivOrders[1]).*shapeloc1d(qPoint, element, degrees[2], iShape2, 
      derivOrders[2])
      matLoc[iShape1, iShape2] = sum(qWeight.*integrand)
    end
  end

  # Return the local "load vector"
  return matLoc
end
