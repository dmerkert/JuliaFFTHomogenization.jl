export RotationMatrixNormal

function RotationMatrixNormal(
                              fromNormal :: Array{R,1},
                              toNormal :: Array{R,1}
                             ) where {R <: Real}

  @argcheck length(fromNormal) == 3
  @argcheck length(toNormal) == 3

  fromNormal = fromNormal/norm(fromNormal)
  toNormal = toNormal/norm(toNormal)

  isapprox(fromNormal,toNormal,rtol=0.0,atol=1e-12) && return eye(3)

  v  = cross(fromNormal,toNormal)
  s  = norm(v,2)
  c  = dot(fromNormal,toNormal)
  vX = [0 -v[3] v[2];v[3] 0 -v[1];-v[2] v[1] 0]

  (eye(3) + vX + vX^2*(1.0-c)/(s^2)) :: Array{R,2}
end
