abstract GreenOperator

type Gamma0 <: GreenOperator end

function mult!{C <: Complex, R <: Real}(
                             strain :: Strain{C},
                             stress :: Stress{C},
                             referenceStiffness :: IsotropicStiffnessTensor{R},
                             FourierIndex :: Array{R,1}
                            )

  mu = referenceStiffness.mu.val
  lambda = referenceStiffness.lambda.val
  SSum = sumabs2(FourierIndex)

  factor1 = 1.0/(4.0mu*SSum)
  factor2 = (lambda+mu)/(mu*(lambda+2.0mu)*SSum^2)

  strain.val = zeros(stress.val);

  for i in [(1,1,1) (2,2,2) (3,3,3) (4,2,3) (5,1,3) (6,1,2)]
    strain.val[i[1]] =  factor2 * 
    (
     FourierIndex[i[2]] * FourierIndex[i[3]] * 
     (
      FourierIndex[1]^2 * stress.val[1] +
      FourierIndex[2]^2 * stress.val[2] +
      FourierIndex[3]^2 * stress.val[3] +
      FourierIndex[2]*FourierIndex[3] * stress.val[4] +
      FourierIndex[1]*FourierIndex[3] * stress.val[5] +
      FourierIndex[1]*FourierIndex[1] * stress.val[6]
     )
    )
  end

  for i = 1:3
    strain.val[i] -= factor1 * 4.0FourierIndex[i]^2 * stress.val[i]
  end

  strain.val[4] -= factor1 *
  (
   FourierIndex[2]*FourierIndex[2]*stress.val[2] +
   FourierIndex[3]*FourierIndex[3]*stress.val[3] +
   2.0FourierIndex[2]*FourierIndex[3]*stress.val[4]
  )

  strain.val[5] -= factor1 *
  (
   FourierIndex[1]*FourierIndex[1]*stress.val[1] +
   FourierIndex[3]*FourierIndex[3]*stress.val[3] +
   2.0FourierIndex[1]*FourierIndex[3]*stress.val[5]
  )

  strain.val[6] -= factor1 *
  (
   FourierIndex[1]*FourierIndex[1]*stress.val[1] +
   FourierIndex[2]*FourierIndex[2]*stress.val[2] +
   2.0FourierIndex[1]*FourierIndex[2]*stress.val[6]
  )
end
