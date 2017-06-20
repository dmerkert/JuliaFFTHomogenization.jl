function mult!{R}(
                  stressField::StressField{R},
                  strainField::StrainField{R},
                  stiffnessField::StiffnessTensorField
              )

  @argcheck size(strainField) == size(stiffnessField)
  stressField.val = strainField.val
  strain = Strain(R)
  stress = Stress(R)

  Rpre = CartesianRange(size(strainField))
  _mult!(stressField,strainField,stiffnessField,strain,stress,Rpre)
end

function multShifted!{R}(
                  stressField::StressField{R},
                  strainField::StrainField{R},
                  stiffnessField::StiffnessTensorField,
                  stiffnessTensor::StiffnessTensor
              )

  size(strainField) == size(stiffnessField) || error("Sizes not compatible")
  stressField.val = strainField.val
  strain = Strain(R)
  stress = Stress(R)

  Rpre = CartesianRange(size(strainField))
  _multShifted!(stressField,
                strainField,
                stiffnessField,
                stiffnessTensor,
                strain,
                stress,
                Rpre)
end




function avg{R}(strainField::StrainField{R})
  sum = Strain(R)
  Iter = CartesianRange(size(strainField))
  _avg!(sum,strainField,Iter)
  sum
end

function avg{R}(stressField::StressField{R})
  sum = Stress(R)
  Iter = CartesianRange(size(stressField))
  _avg!(sum,stressField,Iter)
  sum
end

function FFT!{C <: Complex, R <: Real}(
                                       stressFrequencyField :: StressField{C},
                                       stressField :: StressField{R}
                                      )
  _FFT!(stressFrequencyField,stressField,1:(length(size(stressField.val))-1))
  stressFrequencyField
end

function IFFT!{C <: Complex, R <: Real}(
                                       strainField :: StrainField{R},
                                       strainFrequencyField :: StrainField{C}
                                      )
  _FFT!(strainField,strainFrequencyField,1:(length(size(strainField.val))-1))
  strainField
end

function applyGammaHat{C <: Complex}(
                             stressFrequencyField :: StressField{C},
                             gamma :: GreenOperator,
                             referenceStiffness :: IsotropicStiffnessTensor,
                             lattice :: Lattice
                            )

  strainFrequencyField = StrainField(stressFrequencyField.val)
  strain = Strain(Complex128)
  stress = Stress(Complex128)
  for index in getFrequencyIterator(lattice)
    FourierIndex = getFrequencyPoint(lattice,index)
    stress = stressFrequencyField[index]
    applyGammaHat!(strain,gamma,stress,referenceStiffness,FourierIndex)
    strainFrequencyField[index] = strain
  end
  strainFrequencyField
end

function computePolarization{R <: Real}(strainField :: StrainField{R},
                                        stiffnessField :: StiffnessTensorField,
                                        referenceStiffness :: StiffnessTensor
                                       )
  stressField = StressField(strainField.val)
  multShifted!(stressField,
               strainField,
               stiffnessField,
               referenceStiffness
              )
  stressField
end

function transform{R <: Real}(stressField::StressField{R})
  stressFrequencyField = StressField(Complex128,size(stressField))
  FFT!(stressFrequencyField,stressField)
  stressFrequencyField
end

function inverseTransform{C <: Complex}(strainFrequencyField :: StrainField{C})
  strainField = StrainField(Float64,size(strainFrequencyField))
  IFFT!(strainField,strainFrequencyField)
  strainField
end

function setAverage!{C <: Complex}(strainFrequencyField :: StrainField{C},
                                  macroscopicStrain :: Strain{C}
                                 )

  index = size(strainFrequencyField)
  strainFrequencyField.val[index,:] = macroscopicStrain.val
  strainFrequencyField
end
