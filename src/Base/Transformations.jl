immutable FFT <: Transformation end

function transform!{C <: Complex,
                    R <: Real,
                    F,
                    G}(frequencyField :: SolutionTensorField{C,F},
                       pointField :: SolutionTensorField{R,G},
                       transformation :: FFT,
                       L :: Lattice
                      )

  frequencyField.val = FFT!(frequencyField.val,pointField.val,L)
  frequencyField
end

function transformInverse!{C <: Complex,
                           R <: Real,
                           F,
                           G}(pointField :: SolutionTensorField{R,G},
                              frequencyField :: SolutionTensorField{C,F},
                              transformation :: FFT,
                              L :: Lattice
                             )

  pointField.val = IFFT!(pointField.val,frequencyField.val,L)
  pointField
end

function setAveragingFrequency!{C <: Complex, F}(frequencyField ::
                                                 SolutionTensorField{C,F},
                                                 tensor :: F,
                                                 transformation :: FFT,
                                                 L :: Lattice
                                                )
  setZerothFourierCoefficient!(frequencyField.val,L,tensor.val)
  frequencyField
end
