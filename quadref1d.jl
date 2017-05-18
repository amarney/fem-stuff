function quadref1d(rule = 5)
## -----------------------------------------------------------------------------
# quadref1d: returns quadrature points and weights in [-1, 1]
#
#   Adapted from Jeff Borggaard's code in summer of 2017.
#
#   Usage: qPoint, qWeight = quadref1d(rule)
#
#   Inputs: rule: the quadrature rule to use (default 5)
#
#   Outputs: qPoint: quadrature points located in [-1, 1]
#            qWeight: quadrature weights corresponding to qPoint
## -----------------------------------------------------------------------------

  # Check if inputs make sense
  if rule < 0
    error("quadref1d, rule must be a positive integer")
  end

  # Declare variables
  qPoint = Array(Float64,rule);
  qWeight = Array(Float64,rule);

  # Compute quadrature points and weights on [-1, 1]
  if rule == 1
    qPoint[1] = 0.;
    qWeight[1] = 2.;
  elseif rule == 2     # up to order 3 polynomials exact
    qPoint[1] =-1.0 / sqrt(3.0);
    qPoint[2] =-qPoint[1];
    qWeight[1] = 1.0;
    qWeight[2] = 1.0;
  elseif rule == 3     # up to order 5 polynomials exact
    qPoint[1] =-sqrt(3.0/5.0);
    qPoint[2] = 0.0;
    qPoint[3] =-qPoint[1];
    qWeight[1] = 5.0 / 9.0;
    qWeight[2] = 8.0 / 9.0;
    qWeight[3] = qWeight[1];
  elseif rule == 4    # up to order 7 polynomials exact
    qPoint[1] =-sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0);
    qPoint[2] =-sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0);
    qPoint[3] =-qPoint[2];
    qPoint[4] =-qPoint[1];
    qWeight[1] = 0.5 - 1.0 / ( 6.0 * sqrt(6.0/5.0) );
    qWeight[2] = 0.5 + 1.0 / ( 6.0 * sqrt(6.0/5.0) );
    qWeight[3] = qWeight[2];
    qWeight[4] = qWeight[1];
  elseif rule == 5    # up to order 9 polynomials exact
#    qPoint[1] = -0.90617984593866399280
#    qPoint[2] = -0.53846931010568309104
#    qPoint[4] = -qPoint[2]
#    qPoint[5] = -qPoint[1]
#    qWeight[1] =
#    qWeight[3] =
#    qWeight[4] =
#    qWeight[5] =
    qPoint[1] =-sqrt(5.0+4.0*sqrt(5.0/14.0)) / 3.0;
    qPoint[2] =-sqrt(5.0-4.0*sqrt(5.0/14.0)) / 3.0;
    qPoint[3] = 0.0;
    qPoint[4] =-qPoint[2];
    qPoint[5] =-qPoint[1];
    qWeight[1] = 161.0/450.0-13.0/(180.*sqrt(5.0/14.0));
    qWeight[2] = 161.0/450.0+13.0/(180.*sqrt(5.0/14.0));
    qWeight[3] = 128.0/225.0;
    qWeight[4] = qWeight[2];
    qWeight[5] = qWeight[1];
  elseif rule == 6    # up to order 11 polynomials exact
    qPoint[1] = -0.2386191860831969;
    qPoint[2] = -0.6612093864662645;
    qPoint[3] = -0.9324695142031521;
    qPoint[4] = - qPoint[1];
    qPoint[5] = - qPoint[2];
    qPoint[6] = - qPoint[3];
    qWeight[1] = 0.4679139345726910;
    qWeight[2] = 0.3607615730481386;
    qWeight[3] = 0.1713244923791704;
    qWeight[4] = qWeight[1];
    qWeight[5] = qWeight[2];
    qWeight[6] = qWeight[3];
  elseif (rule == 7)    # up to order 13 polynomials exact
    qPoint[1] = -0.9491079123427585;
    qPoint[2] = -0.7415311855993945;
    qPoint[3] = -0.4058451513773972;
    qPoint[4] =  0.0000000000000000;
    qPoint[5] = - qPoint[3];
    qPoint[6] = - qPoint[2];
    qPoint[7] = - qPoint[1];
    qWeight[1] = 0.1294849661688697;
    qWeight[2] = 0.2797053914892766;
    qWeight[3] = 0.3818300505051189;
    qWeight[4] = 0.4179591836734694;
    qWeight[5] = qWeight[3];
    qWeight[6] = qWeight[2];
    qWeight[7] = qWeight[1];
  elseif (rule == 8)    # up to order 15 polynomials exact
    qPoint[1] = -0.9602898564975363;
    qPoint[2] = -0.7966664774136267;
    qPoint[3] = -0.5255324099163290;
    qPoint[4] = -0.1834346424956498;
    qPoint[5] = - qPoint[4];
    qPoint[6] = - qPoint[3];
    qPoint[7] = - qPoint[2];
    qPoint[8] = - qPoint[1];
    qWeight[1] = 0.1012285362903763;
    qWeight[2] = 0.2223810344533745;
    qWeight[3] = 0.3137066458778873;
    qWeight[4] = 0.3626837833783620;
    qWeight[5] = qWeight[4];
    qWeight[6] = qWeight[3];
    qWeight[7] = qWeight[2];
    qWeight[8] = qWeight[1];
  elseif (rule == 9)     #  up to order 17 polynomials exact
    qPoint[1] = -0.9681602395076261;
    qPoint[2] = -0.8360311073266358;
    qPoint[3] = -0.6133714327005904;
    qPoint[4] = -0.3242534234038089;
    qPoint[5] =  0.0000000000000000;
    qPoint[6] = - qPoint[4];
    qPoint[7] = - qPoint[3];
    qPoint[8] = - qPoint[2];
    qPoint[9] = - qPoint[1];
    qWeight[1] = 0.0812743883615744;
    qWeight[2] = 0.1806481606948574;
    qWeight[3] = 0.2606106964029354;
    qWeight[4] = 0.3123470770400029;
    qWeight[5] = 0.3302393550012598;
    qWeight[6] = qWeight[4];
    qWeight[7] = qWeight[3];
    qWeight[8] = qWeight[2];
    qWeight[9] = qWeight[1];
  elseif (rule == 10)    #  up to order 19 polynomials exact
    qPoint[1] = -0.9739065285171717;
    qPoint[2] = -0.8650633666889845;
    qPoint[3] = -0.6794095682990244;
    qPoint[4] = -0.4333953941292472;
    qPoint[5] = -0.1488743389816312;
    qPoint[6] = - qPoint[5];
    qPoint[7] = - qPoint[4];
    qPoint[8] = - qPoint[3];
    qPoint[9] = - qPoint[2];
    qPoint[10] = - qPoint[1];
    qWeight[1] = 0.0666713443086881;
    qWeight[2] = 0.1494513491505806;
    qWeight[3] = 0.2190863625159820;
    qWeight[4] = 0.2692667193099963;
    qWeight[5] = 0.2955242247147529;
    qWeight[6] = qWeight[5];
    qWeight[7] = qWeight[4];
    qWeight[8] = qWeight[3];
    qWeight[9] = qWeight[2];
    qWeight[10] = qWeight[1];
  elseif (rule == 11)     #  up to order 21 polynomials exact
    qPoint[1] = -0.9782286581460570;
    qPoint[2] = -0.8870625997680953;
    qPoint[3] = -0.7301520055740494;
    qPoint[4] = -0.5190961292068118;
    qPoint[5] = -0.2695431559523450;
    qPoint[6] =  0.0000000000000000;
    qPoint[7] = - qPoint[5];
    qPoint[8] = - qPoint[4];
    qPoint[9] = - qPoint[3];
    qPoint[10] = - qPoint[2];
    qPoint[11] = - qPoint[1];
    qWeight[1] = 0.0556685671161737;
    qWeight[2] = 0.1255803694649046;
    qWeight[3] = 0.1862902109277343;
    qWeight[4] = 0.2331937645919905;
    qWeight[5] = 0.2628045445102467;
    qWeight[6] = 0.2729250867779006;
    qWeight[7] = qWeight[5];
    qWeight[8] = qWeight[4];
    qWeight[9] = qWeight[3];
    qWeight[10] = qWeight[2];
    qWeight[11] = qWeight[1];
  else
    error("quadref1d, rule not supported")
  end

  # Return quadrature points and weights on [-1, 1]
  return qPoint, qWeight
end
