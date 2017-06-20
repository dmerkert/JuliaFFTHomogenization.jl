
#    (d1,d2,d3) = DepolarizationFactors(l1,l2,l3) calculates the depolarization
# factors for a ellipsoid using the lengths of the semi-axes
#
# Currently only the following special cases are supported:
#   prolate spheroid:    l1 >= l2 == l3
#   oblate spheroid:     l1 <= l2 == l3
#   elliptical cylinder: l3 == inf
#
# INPUT:
#   l1: length of semi-axis in 1-direction
#   l2: length of semi-axis in 2-direction
#   l3: length of semi-axis in 3-direction
#
# OUTPUT
#   d1: depolarization factor in 1-direction
#   d2: depolarization factor in 2-direction
#   d3: depolarization factor in 3-direction


function DepolarizationFactors{R <: Real}(l1::R,l2::R,l3::R)
  @argcheck l1 > 0
  @argcheck l2 > 0
  @argcheck l3 > 0
  @argcheck (l1 >= l2 == l3) || (l1 <= l2 == l3) || (l3 == Inf)

  if l1 == l2 && l2 == l3
    d1 = 1.0/3.0
    d2 = d1
    d3 = d1
    return (d1,d2,d3)
  elseif l1 > l2 && l2 == l3  #prolate spheroid
    epsilon = sqrt(1.0 - (l2/l1)^2)
    d1 = (1.0-epsilon^2)/epsilon^2 *
    (
     1.0/(2.0epsilon)*log((1.0+epsilon)/(1.0-epsilon))-1.0
    )
    d2 = 0.5 - 0.5d1
    d3 = d2
    return (d1,d2,d3)
  elseif l1 < l2 && l2 == l3  #oblate spheroid
    epsilon = sqrt(1.0 - (l1/l2)^2)
    d1 = 1.0/epsilon^2 * (1.0 - sqrt(1.0-epsilon^2)/epsilon*asin(epsilon))
    d2 = 0.5 - 0.5d1
    d3 = d2
    return (d1,d2,d3)
  elseif l3 == Inf  #elliptical cylinder
    d1 = l2/(l1+l2)
    d2 = l1/(l1+l2)
    d3 = 0.0
    return (d1,d2,d3)
  end

end

function CartesianToEllipsoidal{R <: Real}(x::Array{R,1},c1::R,c2::R,c3::R)
  @argcheck length(x) == 3
  @argcheck c1 > 0
  @argcheck c2 > 0
  @argcheck c3 > 0

  x1 = x[1]
  x2 = x[2]
  x3 = x[3]

  if (c3 == Inf)
    #polynomial coefficients of polynomial a*p^2+b*p+c=0
    a = -1.0
    b = x1.^2 + x2.^2 - c1^2 - c2^2
    c = x1.^2.*c2^2 + x2.^2.*c1^2 - c1^2*c2^2

    y = real(roots(Poly([c,b,a])))
    p = max(y...)
    @assert p+c1^2 >= -1e-12

    #account for numerical errors
    (p < -c1^2) && (p = -c1^2)
  else
    #polynomial coefficients of polynomial a*p^3+b*p^2+c*p+d=0
    a = -1
    b = x1.^2+x2.^2+x3.^2-c1^2-c2^2-c3^2
    c = x1.^2*c2^2 + x1.^2*c3^2 +
    x2.^2*c1^2 + x2.^2*c3^2 +
    x3.^2*c1^2 + x3.^2*c2^2 -
    c1^2*c2^2 - c1^2*c3^2 - c2^2*c3^2
    d = x1.^2*c2^2*c3^2 +
    x2.^2*c1^2*c3^2 +
    x3.^2*c1^2*c2^2 -
    c1^2*c2^2*c3^2


    y = real(roots(Poly([d,c,b,a])))
    p = max(y...)

    @assert all(abs(
                    x1.^2./(c1^2+p) +
                    x2.^2./(c2^2+p) +
                    x3.^2./(c3^2+p) -
                    1.0
                   ) < 1e-12)

  end
end

