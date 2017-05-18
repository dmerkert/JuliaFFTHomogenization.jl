function mult!{R}(stress::Stress{R},stiffness::IsotropicStiffnessTensor{R},strain::Strain{R})
  stress.val[1] = (2.0stiffness.mu.val+stiffness.lambda.val)*strain.val[1]+
  stiffness.lambda.val*(strain.val[2]+strain.val[3])
  stress.val[2] = (2.0stiffness.mu.val+stiffness.lambda.val)*strain.val[2]+
  stiffness.lambda.val*(strain.val[1]+strain.val[3])
  stress.val[3] = (2.0stiffness.mu.val+stiffness.lambda.val)*strain.val[3]+
  stiffness.lambda.val*(strain.val[1]+strain.val[2])
  stress.val[4] = stiffness.mu.val*strain.val[4]
  stress.val[5] = stiffness.mu.val*strain.val[5]
  stress.val[6] = stiffness.mu.val*strain.val[6]
  stress
end

function mult!{R}(stress::Stress{R},stiffness::TransversalIsotropicZ{R},strain::Strain{R})
  gamma = 1.0 / (1.0 - stiffness.nu_p.val^2 -
                  2.0stiffness.nu_pt.val*stiffness.nu_tp.val -
                  2.0stiffness.nu_p.val*stiffness.nu_pt.val*stiffness.nu_tp.val);
  c11 = stiffness.E_p.val*(1.0-stiffness.nu_pt.val*stiffness.nu_tp.val)*gamma;
  c33 = stiffness.E_t.val*(1.0-stiffness.nu_p.val^2)*gamma;
  c12 = stiffness.E_p.val*(stiffness.nu_p.val+stiffness.nu_pt.val*stiffness.nu_tp.val)*gamma;
  c13 = stiffness.E_p.val*(stiffness.nu_tp.val+stiffness.nu_p.val*stiffness.nu_tp.val)*gamma;
  c44 = stiffness.mu_t.val;

  stress.val[1] = c11*strain.val[1]+c12*strain.val[2]+c13*strain.val[3]
  stress.val[2] = c12*strain.val[1]+c11*strain.val[2]+c13*strain.val[3]
  stress.val[3] = c13*strain.val[1]+c13*strain.val[2]+c33*strain.val[3]
  stress.val[4] = c44/2.0*strain.val[4]
  stress.val[5] = c44/2.0*strain.val[5]
  stress.val[6] = (c11-c12)/4.0*strain.val[6]
  stress
end

function mult!{R}(
                  stressField::StressField{R},
                  strainField::StrainField{R},
                  stiffnessField::StiffnessTensorField
              )

  size(strainField) == size(stiffnessField) || error("Sizes not compatible")
  stressField.val = strainField.val
  strain = Strain(R)
  stress = Stress(R)

  Rpre = CartesianRange(size(strainField))
  _mult!(stressField,strainField,stiffnessField,strain,stress,Rpre)
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
