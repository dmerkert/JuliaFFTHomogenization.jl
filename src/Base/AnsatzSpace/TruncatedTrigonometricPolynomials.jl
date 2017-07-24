export TruncatedTrigonometricPolynomials

type TruncatedTrigonometricPolynomials <: AnsatzSpace end

function coefficients2Values!(
                              pointField :: SolutionTensorField{R,G,N},
                              space :: TruncatedTrigonometricPolynomials,
                              L :: Lattice
                             ) where {
                                      R <: AbstractFloat,
                                      G,
                                      N
                                     }
  #nothing to do here
end

function values2Coefficients!(
                              pointField :: SolutionTensorField{R,G,N},
                              space :: TruncatedTrigonometricPolynomials,
                              L :: Lattice
                             ) where {
                                      R <: AbstractFloat,
                                      G,
                                      N
                                     }
  #nothing to do here
end
