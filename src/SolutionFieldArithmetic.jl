function _mult!(
                fluxField::FluxSolutionTensorField,
                gradientField::GradientSolutionTensorField,
                coefficientField::CoefficientTensorField,
                gradient::GradientSolutionTensor,
                flux::FluxSolutionTensor,
                Rpre
               )

  fluxField.val = gradientField.val

  for Ipre in Rpre
    get!(gradient,gradientField,Ipre)
    mult!(flux,get(coefficientField,Ipre),gradient)
    set!(flux,fluxField,Ipre)
  end
  fluxField
end

function _multShifted!(
                fluxField::FluxSolutionTensorField,
                gradientField::GradientSolutionTensorField,
                coefficientField::CoefficientTensorField,
                coefficient::CoefficientTensor,
                gradient::GradientSolutionTensor,
                flux::FluxSolutionTensor,
                Rpre
               )

  fluxField.val = gradientField.val

  for Ipre in Rpre
    get!(gradient,gradientField,Ipre)
    mult!(flux,get(coefficientField,Ipre)-coefficient,gradient)
    set!(flux,fluxField,Ipre)
  end
  fluxField
end


function _avg!(sum::SolutionTensor,field::SolutionTensorField,Iter)
  for I in Iter
    add!(sum,field,I)
  end
  sum = sum/prod(last(Iter).I)
end

function _FFT!(
               fluxFrequencyField :: FluxSolutionTensorField,
               fluxField :: FluxSolutionTensorField,
               dims
              )
  fluxFrequencyField.val = fft(fluxField.val,dims)
  fluxFrequencyField
end

function _IFFT!(
                gradientField :: GradientSolutionTensorField,
                gradientFrequencyField :: GradientSolutionTensorField,
                dims
               )
  gradientField.val = ifft(gradientFrequencyField.val,dims)
  gradientField
end

function _setZerothFourierCoefficient!(
                                       gradientFrequencyField ::
                                       GradientSolutionTensorField,
                                       gradient :: GradientSolutionTensor,
                                       IRest
                                      )
  
  set!(gradientFrequencyField,gradient,ones(1:lenght(size(gradientFrequencyField))))
  gradientFrequencyField
end

function +{T <: SolutionTensor}(A :: T,B :: T)
  T(A.val+B.val)
end

function -{T <: SolutionTensor}(A :: T,B :: T)
  T(A.val-B.val)
end

function +{T <: SolutionTensorField}(A :: T,B :: T)
  T(A.val+B.val)
end

function -{T <: SolutionTensorField}(A :: T,B :: T)
  T(A.val-B.val)
end

function norm{T <: SolutionTensorField}(A :: T)
  sqrt(sum(abs(A.val[:]).^2))
end
