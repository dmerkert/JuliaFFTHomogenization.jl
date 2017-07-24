export laminateStiffness

function laminateStiffness(
                            C1 :: T1,
                            C2 :: T2,
                            amountC1 :: R,
                            amountC2 :: R,
                            n :: Array{R,1}
                           ) where {
                                    T1 <: StiffnessTensor,
                                    T2 <: StiffnessTensor,
                                    R <: AbstractFloat
                                   }

  amountC1 /= (amountC1+amountC2)
  amountC2 /= (amountC1+amountC2)

  n /= norm(n)

  P = zeros(6,6)

  P[1,1] = _evaluateProjection(n,1,1,1,1)
  P[1,2] = _evaluateProjection(n,1,1,2,2)
  P[1,3] = _evaluateProjection(n,1,1,3,3)
  P[1,4] = 2.0_evaluateProjection(n,1,1,2,3)
  P[1,5] = 2.0_evaluateProjection(n,1,1,1,3)
  P[1,6] = 2.0_evaluateProjection(n,1,1,1,2)

  P[2,1] = _evaluateProjection(n,2,2,1,1)
  P[2,2] = _evaluateProjection(n,2,2,2,2)
  P[2,3] = _evaluateProjection(n,2,2,3,3)
  P[2,4] = 2.0_evaluateProjection(n,2,2,2,3)
  P[2,5] = 2.0_evaluateProjection(n,2,2,1,3)
  P[2,6] = 2.0_evaluateProjection(n,2,2,1,2)

  P[3,1] = _evaluateProjection(n,3,3,1,1)
  P[3,2] = _evaluateProjection(n,3,3,2,2)
  P[3,3] = _evaluateProjection(n,3,3,3,3)
  P[3,4] = 2.0_evaluateProjection(n,3,3,2,3)
  P[3,5] = 2.0_evaluateProjection(n,3,3,1,3)
  P[3,6] = 2.0_evaluateProjection(n,3,3,1,2)

  P[4,1] = 2.0_evaluateProjection(n,2,3,1,1)
  P[4,2] = 2.0_evaluateProjection(n,2,3,2,2)
  P[4,3] = 2.0_evaluateProjection(n,2,3,3,3)
  P[4,4] = 4.0_evaluateProjection(n,2,3,2,3)
  P[4,5] = 4.0_evaluateProjection(n,2,3,1,3)
  P[4,6] = 4.0_evaluateProjection(n,2,3,1,2)

  P[5,1] = 2.0_evaluateProjection(n,1,3,1,1)
  P[5,2] = 2.0_evaluateProjection(n,1,3,2,2)
  P[5,3] = 2.0_evaluateProjection(n,1,3,3,3)
  P[5,4] = 4.0_evaluateProjection(n,1,3,2,3)
  P[5,5] = 4.0_evaluateProjection(n,1,3,1,3)
  P[5,6] = 4.0_evaluateProjection(n,1,3,1,2)

  P[6,1] = 2.0_evaluateProjection(n,1,2,1,1)
  P[6,2] = 2.0_evaluateProjection(n,1,2,2,2)
  P[6,3] = 2.0_evaluateProjection(n,1,2,3,3)
  P[6,4] = 4.0_evaluateProjection(n,1,2,2,3)
  P[6,5] = 4.0_evaluateProjection(n,1,2,1,3)
  P[6,6] = 4.0_evaluateProjection(n,1,2,1,2)

  Id = eye(6)
  Id[4,4] = 0.5
  Id[5,5] = 0.5
  Id[6,6] = 0.5

  eigC1 = eig(C1)
  eigC2 = eig(C2)

  λ = max(max(eigC1...),max(eigC2...)) + 1.0

  CEffective =
  amountC1*
    inv(
        P + λ*inv(convert(AnisotropicStiffnessTensor,C1).C - λ*Id)
       ) +
  amountC2*
    inv(
        P + λ*inv(convert(AnisotropicStiffnessTensor,C2).C - λ*Id)
    )

  CEffective = inv(CEffective)
  CEffective -= P
  CEffective ./= λ
  CEffective = inv(CEffective)
  CEffective += λ*Id
  AnisotropicStiffnessTensor(CEffective)
end

function _evaluateProjection(
                              n :: Array{R,1},
                              i :: I,
                              j :: I,
                              l :: I,
                              m :: I
                             ) where {
                                      R <: AbstractFloat,
                                      I <: Integer
                                     }
  0.5(
      (j==l)*n[i]*n[m] + 
      (j==m)*n[i]*n[l] + 
      (i==l)*n[j]*n[m] + 
      (i==m)*n[j]*n[l]
     ) - n[i]*n[j]*n[l]*n[m]
end

