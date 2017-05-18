function shaperef1d(t, degree, iShape, derivOrder)
## -----------------------------------------------------------------------------
# shaperef1d: evaluates 1d fe shape function (derivatives) on [-1, 1]
#
#   Written by Angelo Marney in the summer of 2017.
#
#   Usage: f = shaperef1d(t, degree, iShape, derivOrder)
#
#   Inputs: t: an array of coordinates (between -1 and 1) to evaluate at
#           degree: the degree (aka order) of the finite element space
#           iShape: index for specifying local shape function
#           derivOrder: the derivative order (default 0)
#
#   Outputs: f: finite element shape function (derivative) values at t
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if degree < 1
    error("shaperef1d, finite element space can't have degree less than one")
  end
#  if (t[1] != -1 || t[end] != 1)
#    error("shaperef1d, given t interval does not go from -1 to 1")
#  end    # points just need to lie in reference interval
  if iShape < 1
    error("shaperef1d, shape index can't be less than one")
  end
  if derivOrder < 0
    error("shaperef1d, derivative order can't be negative")
  end

# Evaluate the function at t using the correct formula
if degree == 1
  if iShape == 1
    if derivOrder == 0
      f = (1 - t)/2
    elseif derivOrder == 1
      f = (-1/2)*ones(size(t))
    elseif derivOrder > 1
      f = zeros(size(t))
    end
  elseif iShape == 2
    if derivOrder == 0
      f = (1 + t)/2
    elseif derivOrder == 1
      f = (1/2)*ones(size(t));
    elseif derivOrder > 1
      f = zeros(size(t));
    end
  else
    error("shaperef1d, for a first degree finite element space, iShape cannot exceed
    two")
  end
elseif degree == 2
  if iShape == 1
    if derivOrder == 0
      f = -(1/2)*(1 - t).*t
    elseif derivOrder == 1
      f = -1/2 + t
    elseif derivOrder == 2
      f = ones(size(t))
    elseif derivOrder > 2
      f = zeros(size(t))
    end
  elseif iShape == 2
    if derivOrder == 0
      f = 1 - t.^2
    elseif derivOrder == 1
      f = -2*t
    elseif derivOrder == 2
      f = -2*ones(size(t))
    elseif derivOrder > 2
      f = zeros(size(t))
    end
  elseif iShape == 3
    if derivOrder == 0
      f = (1/2)*t.*(1 + t)
    elseif derivOrder == 1
      f = 1/2 + t
    elseif derivOrder == 2
      f = ones(size(t))
    elseif derivOrder > 2
      f = zeros(size(t))
    end
  else
    error("shaperef1d, for a second degree finite element space, iShape cannot exceed
    three")
  end
elseif degree == 3
  if iShape == 1
    if derivOrder == 0
      f = (1/16)*(-1 + t + 9*t.^2 - 9*t.^3)
    elseif derivOrder == 1
      f = (1/16)*(1 + 18*t - 27*t.^2)
    elseif derivOrder == 2
      f = (9/8)*(1 - 3.*t)
    elseif derivOrder == 3
      f = -(27/8)*ones(size(t))
    elseif derivOrder > 3
      f = zeros(size(t))
    end
  elseif iShape == 2
    if derivOrder == 0
      f = (9/16)*(1 - 3*t - t.^2 + 3*t.^3)
    elseif derivOrder == 1
      f = (9/16)*(-3 - 2*t + 9*t.^2)
    elseif derivOrder == 2
      f = (9/8)*(-1 + 9.*t)
    elseif derivOrder == 3
      f = (81/8)*ones(size(t))
    elseif derivOrder > 3
      f = zeros(size(t))
    end
  elseif iShape == 3
    if derivOrder == 0
      f = -(9/16)*(-1 - 3*t + t.^2 + 3*t.^3);
    elseif derivOrder == 1
      f = -(9/16)*(-3 + 2*t + 9*t.^2);
    elseif derivOrder == 2
      f = -(9/8)*(1 + 9*t);
    elseif derivOrder == 3
      f = -(81/8)*ones(size(t));
    elseif derivOrder > 3
      f = zeros(size(t));
    end
  elseif iShape == 4
    if derivOrder == 0
      f = (1/16)*(-1 - t + 9*t.^2 + 9*t.^3);
    elseif derivOrder == 1
      f = (1/16)*(-1 + 18*t + 27*t.^2);
    elseif derivOrder == 2
      f = (9/8)*(1 + 3*t);
    elseif derivOrder == 3
      f = (27/8)*ones(size(t));
    elseif derivOrder > 3
      f = zeros(size(t));
    else
      error("shaperef1d, for a third degree finite element space, iShape cannot exceed
      four.")
    end
  end
else
    error("shaperef1d, finite element spaces of degree higher than 3
    are not yet supported")
end

  # return the function value
  return f
end
