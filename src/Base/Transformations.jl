export FFTTransformation,
       transform!,
       transformInverse!,
       setAveragingFrequency!

immutable FFTTransformation <: Transformation end

function transform!{C <: Complex,
                    R <: Real,
                    F,
                    G}(frequencyField :: SolutionTensorField{C,F},
                       pointField :: SolutionTensorField{R,G},
                       transformation :: FFTTransformation,
                       L :: Lattice
                      )

  frequencyField.val =
  patternfft(pointField.val,L,LastDimensionsFFT)
  frequencyField
end

function transformInverse!{C <: Complex,
                           R <: Real,
                           F,
                           G}(pointField :: SolutionTensorField{R,G},
                              frequencyField :: SolutionTensorField{C,F},
                              transformation :: FFTTransformation,
                              L :: Lattice
                             )

  pointField.val = patternifft(frequencyField.val,L,LastDimensionsFFT)
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
  setZerothFourierCoefficient!(frequencyField.val,L,tensor.val*L.m,LastDimensionsFFT)
  frequencyField
end
