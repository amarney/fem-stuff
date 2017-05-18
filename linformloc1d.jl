function linformloc1d(fName, element, degree, derivOrder = 0,
  rule = degree + 1)
## -----------------------------------------------------------------------------
# linformloc1d: computes local "load vector" from \int { fName*LocShape }
#
#   Written by Angelo Marney in the summer of 2017.
#
#   Usage: vecLoc = linformloc1d(fName, element, degree, derivOrder, rule)
#
#
#   Inputs: fName: the function in the linear form
#           element: an array of coordinates for the element containing x
#           degree: the degree (aka order) of the finite element space
#           derivOrder: the derivative order (default 0)
#           rule: the quadrature rule to use (default degree + 1)
#
#   Outputs: vecLoc: the local "load vector" associated with a linear form
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if length(element) != 2
    error("linformloc1d, element must have exactly two endpoints")
  end
  if degree < 1
    error("linformloc1d, finite element space can't have degree less
    than one")
  end
  if derivOrder < 0
    error("linformloc1d, derivative order can't be negative")
  end
  if rule < 0
    error("linformloc1d, rule must be a positive integer")
  end

  # Compute the quadrature points and weights
  qPoint, qWeight = quadloc1d(element, rule)

  # Integrate over each local shape function for the given element
  vecLoc = Array(Float64, degree + 1, 1);
  for iShape = 1:(degree + 1)
    integrand = fName(qPoint).*shapeloc1d(qPoint, element, degree, iShape,
    derivOrder)
    vecLoc[iShape] = sum(qWeight.*integrand)
  end

  # Return the local "load vector"
  return vecLoc
end
