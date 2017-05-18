function evalfe1d(x, uLocal, element, degree, derivOrder = 0)
## -----------------------------------------------------------------------------
# evalfe1d: evaluates 1d finite element function (derivatives) locally
#
#   Written by Angelo Marney in the summer of 2017.
#
#   Usage: u = evalfe1d(x, uLocal, element, degree, derivOrder)
#
#   Inputs: x: an array of coordinates to evaluate at
#           uLocal: array for local coefficients for element
#           element: an array of coordinates for the element containing x
#           degree: the degree (aka order) of the finite element space
#           derivOrder: the derivative order (default 0)
#
#   Outputs: u: finite element shape function (derivative) values at x
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if length(element) != 2
    error("evalfe1d, element must have exactly two endpoints")
  end
  if degree < 1
    error("evalfe1d, element space can't have degree less than one")
  end
  if derivOrder < 0
    error("evalfe1d, derivative order can't be negative")
  end
  if length(uLocal) != (degree + 1)
    error("evalfe1d, number of local coefficients must equal the number of shape
    functions")
  end

  # Take a linear combination of the shape functions
  u = zeros(size(x));
  for iShape = 1:(degree + 1)
    u = u + uLocal[iShape]*shapeloc1d(x, element, degree, iShape, derivOrder)
  end

  # Return the finite element shape function values at x
  return u
end
