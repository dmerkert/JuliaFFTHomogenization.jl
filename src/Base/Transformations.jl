export FFTTransformation,
       transform!,
       transformInverse!,
       setAveragingFrequency!

immutable FFTTransformation <: Transformation end

function transform!(
                    frequencyField :: SolutionTensorField{C,F,N},
                    pointField :: SolutionTensorField{R,G,N},
                    transformation :: FFTTransformation,
                    L :: Lattice
                   ) where {
                            C <: Complex,
                            R <: Real,
                            F,
                            G,
                            N
                           }

  frequencyField.val =
  FFT!(frequencyField.val,pointField.val,L,[2:(L.rank+1)...])
  frequencyField
end

function transformInverse!(
                           pointField :: SolutionTensorField{R,G,N},
                           frequencyField :: SolutionTensorField{C,F,N},
                           transformation :: FFTTransformation,
                           L :: Lattice
                          ) where {
                                   C <: Complex,
                                   R <: Real,
                                   F,
                                   G,
                                   N
                                  }

  pointField.val = IFFT!(pointField.val,frequencyField.val,L,[2:(L.rank+1)...])
  pointField
end

function setAveragingFrequency!(
                                frequencyField ::
                                SolutionTensorField{C,F,N},
                                tensor :: F,
                                transformation ::
                                FFTTransformation,
                                L :: Lattice
                               ) where {
                                        C <: Complex,
                                        F,
                                        N
                                       }

  setFourierCoefficient!(frequencyField.val,L,tensor.val*L.m,[2:(L.rank+1)...])
  frequencyField
end
