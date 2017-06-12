function RotationMatrixNormal{R <: Real}(
                                         fromNormal :: Array{R,1},
                                         toNormal :: Array{R,1}
                                        )

  @argcheck length(fromNormal) == 3
  @argcheck length(toNormal) == 3

  fromNormal = fromNormal/norm(fromNormal)
  toNormal = toNormal/norm(toNormal)

  isapprox(fromNormal,toNormal,rtol=0.0,atol=1e-12) && return eye(3)

  v  = cross(fromNormal,toNormal)
  s  = norm(v,2)
  c  = dot(fromNormal,toNormal)
  vX = [0 -v[3] v[2];v[3] 0 -v[1];-v[2] v[1] 0]

  eye(3) + vX + vX^2*(1.0-c)/(s^2)
end

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

      # reduction to depressed cubic by Tschirnhaus transformation to polynomial
      # t^3 + p*t + q = 0
      # with x = t - b/(3a)

      #     p = (3.0.*a.*c-b.^2)./(3.0.*a.^2);
      #     q = (2.0.*b.^3-9.0.*a.*b.*c+27.0.*a.^2.*d)./(27.0.*a.^3);
      #     
      #     %formula by Viete for the special case of pure real roots
        #     t0 = ...
        #         2.0.*sqrt(-p./3.0).*cos(1.0./3.0.*acos(3.0./2.0.*q./p.*sqrt(-3.0./p)) ...
        #                                 %         -2.0*pi*0.0/3.0);
        #     
        #     t1 = ...
        #         2.0.*sqrt(-p./3.0).*cos(1.0./3.0.*acos(3.0./2.0.*q./p.*sqrt(-3.0./p)) ...
        #                                 %         -2.0*pi*1.0/3.0);
        #     
        #     t2 = ...
        #         2.0.*sqrt(-p./3.0).*cos(1.0./3.0.*acos(3.0./2.0.*q./p.*sqrt(-3.0./p)) ...
        #                                 %         -2.0*pi*2.0/3.0);
        #     
        #     %transform back
        #     y1 = t0 - b./(3.0.*a);
        #     y2 = t1 - b./(3.0.*a);
        #     y3 = t2 - b./(3.0.*a);
        #     
        #     %"radius" is the largest solution
        #     p = max([y1;y2;y3]);

        @assert all(abs(
                        x1.^2./(c1^2+p) +
                        x2.^2./(c2^2+p) +
                        x3.^2./(c3^2+p) -
                        1.0
                       ) < 1e-12)

      end
    end

