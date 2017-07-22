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
  patternfft(pointField.val,L,[2:(L.rank+1)...])
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

  pointField.val = real(patternifft(frequencyField.val,L,[2:(L.rank+1)...]))
  pointField
end

function setAveragingFrequency!{C <: Complex, F}(frequencyField ::
                                                 SolutionTensorField{C,F},
                                                 tensor :: F,
                                                 transformation ::
                                                 FFTTransformation,
                                                 L :: Lattice
                                                )
# TODO: Update this function call to not read Zeorth anymore? Also: How large is the zeroth-vectro here?
  setFourierCoefficient!(frequencyField.val,L,tensor.val*L.m,ones(Int, L.rank),[2:(L.rank+1)...])
  frequencyField
end
