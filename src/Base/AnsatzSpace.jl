export AnsatzSpace,
transform!,
transformInverse!,
setAveragingFrequency!,
coefficients2Values!,
values2Coefficients!

abstract type AnsatzSpace end

function transform!(
                    frequencyField :: SolutionTensorField{C,F,N},
                    pointField :: SolutionTensorField{R,G,N},
                    space :: Space,
                    L :: Lattice
                   ) where {
                            C <: Complex,
                            R <: Real,
                            F,
                            G,
                            N,
                            Space <: AnsatzSpace
                           }

  frequencyField.val =
  patternfft(pointField.val,L,[2:(L.rank+1)...])
  frequencyField
end

function transformInverse!(
                           pointField :: SolutionTensorField{R,G,N},
                           frequencyField :: SolutionTensorField{C,F,N},
                           space :: Space,
                           L :: Lattice
                          ) where {
                                   C <: Complex,
                                   R <: Real,
                                   F,
                                   G,
                                   N,
                                   Space <: AnsatzSpace
                                  }

  pointField.val = real(patternifft(frequencyField.val,L,[2:(L.rank+1)...]))
  pointField
end

function setAveragingFrequency!(frequencyField :: SolutionTensorField{C,F,N},
                                tensor :: F,
                                space :: Space,
                                L :: Lattice
                               ) where {
                                        C <: Complex,
                                        F,
                                        N,
                                        Space <: AnsatzSpace
                                       }

  setFourierCoefficient!(frequencyField.val,L,tensor.val*L.m,ones(Int, L.rank),[2:(L.rank+1)...])
  frequencyField
end

function coefficients2Values!(
                              pointField :: SolutionTensorField{R,G,N},
                              space :: Space,
                              L :: Lattice
                             ) where {
                                      R <: AbstractFloat,
                                      G,
                                      N,
                                      Space <: AnsatzSpace
                                     }
  error("Ansatz space not implemented.")
end

function values2Coefficients!(
                              pointField :: SolutionTensorField{R,G,N},
                              space :: Space,
                              L :: Lattice
                             ) where {
                                      R <: AbstractFloat,
                                      G,
                                      N,
                                      Space <: AnsatzSpace
                                     }
  error("Ansatz space not implemented.")
end

function mult!(c :: SolutionTensorField{R,F1,M},
               G :: GreenOperator,
               b :: SolutionTensorField{R,F2,M},
               referenceCoefficient :: C,
               L :: Lattice,
               space :: Space
              ) where {
                       R,
                       F1,
                       F2,
                       M,
                       C<:CoefficientTensor,
                       Space <: AnsatzSpace
                      }
  error("mult! not implemented.")
end
