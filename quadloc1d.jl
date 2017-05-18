function quadloc1d(element, rule = 5)
## -----------------------------------------------------------------------------
# quadloc1d: returns quadrature points and weights in element
#
#   Written by Angelo Marney in the summer of 2017.
#
#   Usage: qPoint, qWeight = quadloc1d(rule)
#
#   Inputs: rule: the quadrature rule to use (default 5)
#
#   Outputs: qPoint: quadrature points located in [-1, 1]
#            qWeight: quadrature weights corresponding to qPoint
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if length(element) != 2
    error("quadloc1d, element must have exactly two endpoints")
  end
  if rule < 0
    error("quadloc1d, rule must be a positive integer")
  end

  # Compute qPoint and qWeight on [-1, 1]
  qPoint, qWeight = quadref1d(rule)

  # Map to given element
  qPoint = .5*(element[2] - element[1])*qPoint + .5*(element[1] + element[2])
  qWeight = .5*(element[2] - element[1])*qWeight

  # return quadrature points and weights on element
  return qPoint, qWeight
end
