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
  FFT!(frequencyField.val,pointField.val,L,LastDimensionsFFT)
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

  pointField.val = IFFT!(pointField.val,frequencyField.val,L,LastDimensionsFFT)
  pointField
end

function setAveragingFrequency!{C <: Complex, F}(frequencyField ::
                                                 SolutionTensorField{C,F},
                                                 tensor :: F,
                                                 transformation ::
                                                 FFTTransformation,
                                                 L :: Lattice
                                                )
  setZerothFourierCoefficient!(frequencyField.val,L,tensor.val*L.m,LastDimensionsFFT)
  frequencyField
end
