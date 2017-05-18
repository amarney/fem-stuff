function shapeloc1d(x, element, degree, iShape, derivOrder = 0)
## -----------------------------------------------------------------------------
# shapeloc1d: evaluates 1d shape function (derivatives) locally
#
#   Written by Angelo Marney in the summer of 2017.
#
#   Usage: f = shapeloc1d(x, element, degree, derivOrder, iShape)
#
#   Inputs: x: an array of coordinates to evaluate at
#           element: an array of coordinates for the element containing x
#           degree: the degree (aka order) of the finite element space
#           iShape: index for specifying local shape function
#           derivOrder: the derivative order (default 0)
#
#   Outputs: f: shape function (derivative) values at x
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if length(element) != 2
    error("shapeloc1d, element must have exactly two endpoints")
  end
  if degree < 1
    error("shapeloc1d, finite element space can't have degree less than one")
  end
  if iShape < 1
    error("shapeloc1d, shape index cannot be less than one")
  end
  if derivOrder < 0
    error("shapeloc1d, derivative order can't be negative")
  end
  
  # Do an affine transformation to the reference element [-1, 1]
  t = (2/(element[2] - element[1]))*(x - element[1]) - 1

  # Evaluate 1d finite element shape function (derivatives) on reference element
  f = shaperef1d(t, degree, iShape, derivOrder)

  # Transform back to the local element, and return the shape function values
  return f*(2/(element[2] - element[1]))^derivOrder

end
