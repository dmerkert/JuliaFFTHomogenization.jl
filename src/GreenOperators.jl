abstract GreenOperator

type Gamma0 <: GreenOperator end

function applyGammaHat!{C <: Complex, R <: Real}(
                             strain :: Strain{C},
                             gamma :: Gamma0,
                             stress :: Stress{C},
                             referenceStiffness :: IsotropicStiffnessTensor{R},
                             FourierIndex :: Array{R,1}
                            )

  if norm(FourierIndex) == 0.0
    strain.val = zeros(strain.val)
    return
  end

  mu = referenceStiffness.mu.val
  lambda = referenceStiffness.lambda.val

  SSum = sumabs2(FourierIndex)
  factor1 = 1.0/(4.0mu*SSum)
  factor2 = (lambda+mu)/((lambda+2.0mu)*mu*SSum^2)

  strain.val = zeros(strain.val)

  gamma = zeros((6,6))

  gamma[1,1] = _evalGamma0Hat(factor1,factor2,FourierIndex,1,1,1,1)
  gamma[1,2] = _evalGamma0Hat(factor1,factor2,FourierIndex,1,1,2,2)
  gamma[1,3] = _evalGamma0Hat(factor1,factor2,FourierIndex,1,1,3,3)
  gamma[1,4] = _evalGamma0Hat(factor1,factor2,FourierIndex,1,1,2,3)
  gamma[1,5] = _evalGamma0Hat(factor1,factor2,FourierIndex,1,1,1,3)
  gamma[1,6] = _evalGamma0Hat(factor1,factor2,FourierIndex,1,1,1,2)

  gamma[2,1] = gamma[1,2]
  gamma[2,2] = _evalGamma0Hat(factor1,factor2,FourierIndex,2,2,2,2)
  gamma[2,3] = _evalGamma0Hat(factor1,factor2,FourierIndex,2,2,3,3)
  gamma[2,4] = _evalGamma0Hat(factor1,factor2,FourierIndex,2,2,2,3)
  gamma[2,5] = _evalGamma0Hat(factor1,factor2,FourierIndex,2,2,1,3)
  gamma[2,6] = _evalGamma0Hat(factor1,factor2,FourierIndex,2,2,1,2)

  gamma[3,1] = gamma[1,3]
  gamma[3,2] = gamma[2,3]
  gamma[3,3] = _evalGamma0Hat(factor1,factor2,FourierIndex,3,3,3,3)
  gamma[3,4] = _evalGamma0Hat(factor1,factor2,FourierIndex,3,3,2,3)
  gamma[3,5] = _evalGamma0Hat(factor1,factor2,FourierIndex,3,3,1,3)
  gamma[3,6] = _evalGamma0Hat(factor1,factor2,FourierIndex,3,3,1,2)

  gamma[4,1] = gamma[1,4]
  gamma[4,2] = gamma[2,4]
  gamma[4,3] = gamma[3,4]
  gamma[4,4] = _evalGamma0Hat(factor1,factor2,FourierIndex,2,3,2,3)
  gamma[4,5] = _evalGamma0Hat(factor1,factor2,FourierIndex,2,3,1,3)
  gamma[4,6] = _evalGamma0Hat(factor1,factor2,FourierIndex,2,3,1,2)

  gamma[5,1] = gamma[1,5]
  gamma[5,2] = gamma[2,5]
  gamma[5,3] = gamma[3,5]
  gamma[5,4] = gamma[4,5]
  gamma[5,5] = _evalGamma0Hat(factor1,factor2,FourierIndex,1,3,1,3)
  gamma[5,6] = _evalGamma0Hat(factor1,factor2,FourierIndex,1,3,1,2)

  gamma[6,1] = gamma[1,6]
  gamma[6,2] = gamma[2,6]
  gamma[6,3] = gamma[3,6]
  gamma[6,4] = gamma[4,6]
  gamma[6,5] = gamma[5,6]
  gamma[6,6] = _evalGamma0Hat(factor1,factor2,FourierIndex,1,2,1,2)

  #Corrections for Voigt notation
  gamma[4:6,1:6] = 2.0gamma[4:6,1:6]
  gamma[1:6,4:6] = 2.0gamma[1:6,4:6]

  strain.val = gamma*stress.val
end

@inline function _evalGamma0Hat(factor1,factor2,FourierIndex, i,j,k,h)
  -(
    factor1*(
             (k==i)*FourierIndex[h]*FourierIndex[j] +
             (h==i)*FourierIndex[k]*FourierIndex[j] +
             (k==j)*FourierIndex[h]*FourierIndex[i] +
             (h==j)*FourierIndex[k]*FourierIndex[i]
            ) - factor2*(
                         FourierIndex[i]*
                         FourierIndex[j]*
                         FourierIndex[k]*
                         FourierIndex[h]
                        )
   )
end
