export Gamma0,
       mult!

immutable Gamma0 <: GreenOperator end

function mult!(
               strain :: Strain{C},
               _gamma :: Gamma0,
               stress :: Stress{C},
               referenceStiffness :: IsotropicStiffnessTensor,
               FourierIndex :: Array{IR,1}
              ) where {
                       C <: Complex,
                       IR <: Union{Integer, AbstractFloat}
                      }

  if norm(FourierIndex) == 0
    strain.val = zeros(strain.val)
  else

    mu = referenceStiffness.mu
    lambda = referenceStiffness.lambda

    factor1 = 1.0/(4.0mu*sum(abs2,FourierIndex)) :: Float64
    factor2 = (lambda+mu)/((lambda+2.0mu)*mu*sum(abs2,FourierIndex)^2) ::
    Float64

    gamma = zeros((6,6))

    gamma[1,1] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,1,1,1,1)
    gamma[1,2] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,1,1,2,2)
    gamma[1,3] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,1,1,3,3)
    gamma[1,4] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,1,1,2,3)
    gamma[1,5] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,1,1,1,3)
    gamma[1,6] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,1,1,1,2)

    gamma[2,1] = gamma[1,2]
    gamma[2,2] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,2,2,2,2)
    gamma[2,3] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,2,2,3,3)
    gamma[2,4] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,2,2,2,3)
    gamma[2,5] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,2,2,1,3)
    gamma[2,6] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,2,2,1,2)

    gamma[3,1] = gamma[1,3]
    gamma[3,2] = gamma[2,3]
    gamma[3,3] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,3,3,3,3)
    gamma[3,4] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,3,3,2,3)
    gamma[3,5] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,3,3,1,3)
    gamma[3,6] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,3,3,1,2)

    gamma[4,1] = gamma[1,4]
    gamma[4,2] = gamma[2,4]
    gamma[4,3] = gamma[3,4]
    gamma[4,4] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,2,3,2,3)
    gamma[4,5] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,2,3,1,3)
    gamma[4,6] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,2,3,1,2)

    gamma[5,1] = gamma[1,5]
    gamma[5,2] = gamma[2,5]
    gamma[5,3] = gamma[3,5]
    gamma[5,4] = gamma[4,5]
    gamma[5,5] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,1,3,1,3)
    gamma[5,6] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,1,3,1,2)

    gamma[6,1] = gamma[1,6]
    gamma[6,2] = gamma[2,6]
    gamma[6,3] = gamma[3,6]
    gamma[6,4] = gamma[4,6]
    gamma[6,5] = gamma[5,6]
    gamma[6,6] = _evalGamma0HatElasticity(factor1,factor2,FourierIndex,1,2,1,2)

    #Corrections for Voigt notation
    gamma[4:6,1:6] = 2.0gamma[4:6,1:6]
    gamma[1:6,4:6] = 2.0gamma[1:6,4:6]

    strain.val = gamma*stress.val
  end
  strain
end

@inline function _evalGamma0HatElasticity(
                                          factor1 :: R,
                                          factor2 :: R,
                                          FourierIndex :: Array{IR,1},
                                          i :: I,
                                          j :: I,
                                          k :: I,
                                          h :: I
                                         ) :: R where {
                                                  R <: AbstractFloat,
                                                  I <: Integer,
                                                  IR <: Union{Integer, AbstractFloat}
                                                 }
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
   ) :: R
end
