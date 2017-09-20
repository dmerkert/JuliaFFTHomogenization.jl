export toVoigtSSymFourthOrder!,
       toInverseVoigtSSymFourthOrder!,
       toVoigtSSymIdentity,
       toMandelFromVoigt!,
       toVoigtFromMandel!

function toVoigtSSymFourthOrder!(T :: Array{R,2},f :: Function) where {R}
  T[1,1] = f(1,1,1,1) :: R
  T[1,2] = f(1,1,2,2) :: R
  T[1,3] = f(1,1,3,3) :: R
  T[1,4] = f(1,1,2,3) :: R
  T[1,5] = f(1,1,1,3) :: R
  T[1,6] = f(1,1,1,2) :: R

  T[2,1] = T[1,2]
  T[2,2] = f(2,2,2,2) :: R
  T[2,3] = f(2,2,3,3) :: R
  T[2,4] = f(2,2,2,3) :: R
  T[2,5] = f(2,2,1,3) :: R
  T[2,6] = f(2,2,1,2) :: R

  T[3,1] = T[1,3]
  T[3,2] = T[2,3]
  T[3,3] = f(3,3,3,3) :: R
  T[3,4] = f(3,3,2,3) :: R
  T[3,5] = f(3,3,1,3) :: R
  T[3,6] = f(3,3,1,2) :: R

  T[4,1] = T[1,4]
  T[4,2] = T[2,4]
  T[4,3] = T[3,4]
  T[4,4] = f(2,3,2,3) :: R
  T[4,5] = f(2,3,1,3) :: R
  T[4,6] = f(2,3,1,2) :: R

  T[5,1] = T[1,5]
  T[5,2] = T[2,5]
  T[5,3] = T[3,5]
  T[5,4] = T[4,5]
  T[5,5] = f(1,3,1,3) :: R
  T[5,6] = f(1,3,1,2) :: R

  T[6,1] = T[1,6]
  T[6,2] = T[2,6]
  T[6,3] = T[3,6]
  T[6,4] = T[4,6]
  T[6,5] = T[5,6]
  T[6,6] = f(1,2,1,2) :: R

  T
end

function toInverseVoigtSSymFourthOrder!(T :: Array{R,2}, f :: Function) where {R}
  toVoigtSSymFourthOrder!(T,f)
  T[4:6,1:6] .*= 2.0
  T[1:6,4:6] .*= 2.0

  T
end

function toMandelFromVoigt!(T::Array{R,2}) where {R}
  T[:,3:6] .*= sqrt(2.0)
  T[3:6,:] .*= sqrt(2.0)
  T
end

function toVoigtFromMandel!(T::Array{R,2}) where {R}
  T[:,3:6] ./= sqrt(2.0)
  T[3:6,:] ./= sqrt(2.0)
  T
end
