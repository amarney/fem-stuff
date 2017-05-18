function errorloc1d(uName, uLocal, element, degree, derivOrder = 0,
  rule = degree + 1)
## -----------------------------------------------------------------------------
# errorloc1d: computes the local error measured by the L^2 norm
#
#   Written by Angelo Marney in the summer of 2017.
#
#   Usage: errorLoc = errorloc1d(uName, uLocal, element, degree,
#                             derivOrder, rule)
#
#   Inputs: uName: the true solution function (or its derivative)
#           uLocal: an array for local coefficients for element
#           element: an array of coordinates for the element containing x
#           degree: the degree (aka order) of the finite element space
#           derivOrder: the derivative order (default 0)
#           rule: the quadrature rule to use (default degree + 1)
#
#   Outputs: errorLoc: the error on the element measured by the L^2 norm
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if length(element) != 2
    error("errorloc1d, element must have exactly two endpoints")
  end
  if degree < 1
    error("errorloc1d, finite element space can't have degree less than one")
  end
  if derivOrder < 0
    error("errorloc1d, derivative order can't be negative")
  end
  if length(uLocal) != (degree + 1)
    error("errorloc1d, number of local coefficients must equal the number of
    shape functions on a single element")
  end
  if rule < 0
    error("errorloc1d, rule must be a positive integer")
  end

  # Compute the quadrature points and weights
  qPoint, qWeight = quadloc1d(element, rule)

  # Evaluate the integrand at the quadrature points
  integrand = (uName(qPoint) - evalfe1d(qPoint, uLocal, element, degree,
  derivOrder)).^2

  # Use the quadrature formula and return the result
  return sum(qWeight.*integrand)
end
