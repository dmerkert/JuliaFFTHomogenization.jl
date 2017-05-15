  abstract CoefficientTensor
  abstract StiffnessTensor <: CoefficientTensor

  abstract CoefficientVector
  abstract ElasticityVector <: CoefficientVector

  abstract CoefficientTensorField
  abstract ElasticityVectorField

  type Strain{R<:AbstractFloat} <: ElasticityVector
    val::Array{R,1}

    function Strain{R}(val::Array{R,1})
      length(val) == 6 || error("Strain vector has wrong size")
      new(val)
    end
  end
  Strain{R}(val::Array{R,1}) = Strain{R}(val)

  function Strain()
    Strain(Array(Float64,6))
  end
  function Strain{R<:AbstractFloat}(val::Array{R,2})
    size(val) == (3,3) || error("Wrong size")
    issymmetric(val) || error("Only symmetric strains supported")

    Strain([val[[1,5,9]];2*val[[8,7,4]]])
  end

  type Stress{R<:AbstractFloat} <: ElasticityVector
    val::Array{R,1}

    function Stress{R}(val::Array{R,1})
      length(val) == 6 || error("Stress vector has wrong size")
      new(val)
    end

  end
  Stress{R}(val::Array{R,1}) = Stress{R}(val)
  function Stress()
    Stress(Array(Float64,6))
  end
  function Stress{R<:AbstractFloat}(val::Array{R,2})
    size(val) == (3,3) || error("Wrong size")
    issymmetric(val) || error("Only symmetric stresses supported")

    Stress([val[[1,5,9]];val[[8,7,4]]])
  end

  type Displacement{R<:AbstractFloat} <: ElasticityVector
    val::Array{R,1}

    function Displacement{R}(val::Array{R,1})
      length(val) == 3 || error("Displacement vector has wrong size")
      new(val)
    end
  end
  Displacement{R}(val::Array{R,1}) = Displacement{R}(val)
  function Displacement()
    Displacement(Array(Float64,3))
  end

  immutable IsotropicStiffnessTensor{R<:AbstractFloat} <: StiffnessTensor
    lambda::R
    mu::R

    function IsotropicStiffnessTensor{R}(lambda::R,mu::R)
      lambda >= 0 || warn("Lambda is negative")
      mu >= 0 || error("Mu must be positive")

      new(lambda,mu)
    end
  end
  IsotropicStiffnessTensor{R}(lambda::R,mu::R) =
    IsotropicStiffnessTensor{R}(lambda,mu)

    immutable TransversalIsotropicZ{R<:AbstractFloat} <: StiffnessTensor
      E_p::R
      E_t::R
      nu_p::R
      nu_pt::R
      nu_tp::R
      mu_t::R
    end
    TransversalIsotropicZ{R}(E_p::R, E_t::R, nu_p::R, nu_pt::R, nu_tp::R,
                             mu_t::R) = TransversalIsotropicZ{R}( E_p, E_t,
                                                                 nu_p, nu_pt,
                                                                 nu_tp, mu_t)


  type StiffnessTensorField{R <: StiffnessTensor} <: CoefficientTensorField
    val::Array{R}
  end
  StiffnessTensorField(dims) = 
    StiffnessTensorField{StiffnessTensor}(Array{StiffnessTensor}(dims))
  StiffnessTensorField(R,dims) = 
    StiffnessTensorField{R}(Array{R}(dims))

  type StrainStressField{R <: AbstractFloat} <: ElasticityVectorField
    val::Array{R}
  end
  StrainStressField(R,dims) = StrainStressField{R}(tuple(dims...,6))


  function mult!{R}(stress::Stress{R},stiffness::IsotropicStiffnessTensor{R},strain::Strain{R})
    stress.val[1] = (2.0*stiffness.mu+stiffness.lambda)*strain.val[1]+
      stiffness.lambda*(strain.val[2]+strain.val[3])
    stress.val[2] = (2.0*stiffness.mu+stiffness.lambda)*strain.val[2]+
      stiffness.lambda*(strain.val[1]+strain.val[3])
    stress.val[3] = (2.0*stiffness.mu+stiffness.lambda)*strain.val[3]+
      stiffness.lambda*(strain.val[1]+strain.val[2])
    stress.val[4] = stiffness.mu*strain.val[4]
    stress.val[5] = stiffness.mu*strain.val[5]
    stress.val[6] = stiffness.mu*strain.val[6]
    stress
  end

  function mult!{R}(stress::Stress{R},stiffness::TransversalIsotropicZ{R},strain::Strain{R})
    gamma = 1.0 ./ (1.0 - stiffness.nu_p.^2 -
                    2.0.*stiffness.nu_pt.*stiffness.nu_tp -
                    2.0.*stiffness.nu_p.*stiffness.nu_pt.*stiffness.nu_tp);
    c11 = stiffness.E_p.*(1.0-stiffness.nu_pt.*stiffness.nu_tp).*gamma;
    c33 = stiffness.E_t.*(1.0-stiffness.nu_p.^2).*gamma;
    c12 = stiffness.E_p.*(stiffness.nu_p+stiffness.nu_pt.*stiffness.nu_tp).*gamma;
    c13 = stiffness.E_p.*(stiffness.nu_tp+stiffness.nu_p.*stiffness.nu_tp).*gamma;
    c44 = stiffness.mu_t;

    stress.val[1] = c11*strain.val[1]+c12*strain.val[2]+c13*strain.val[3]
    stress.val[1] = c12*strain.val[1]+c11*strain.val[2]+c13*strain.val[3]
    stress.val[1] = c13*strain.val[1]+c13*strain.val[2]+c33*strain.val[3]
    stress.val[4] = c44./2.0*strain.val[4]
    stress.val[5] = c44./2.0*strain.val[5]
    stress.val[6] = (c11-c12)./4.0*strain.val[6]
    stress
  end

  function mult!(strainStressField::StrainStressField,
                      stiffnessField::StiffnessTensorField)
    size(strainStressField.val)[1:end-1] == size(stiffnessField.val) || error("Sizes not compatible")

    Rpre = CartesianRange(size(strainStressField.val)[1:end-1])
    _mult!(strainStressField,stiffnessField,Rpre)
  end

  function _mult!(strainStressField::StrainStressField,
                      stiffnessField::StiffnessTensorField,Rpre)
    strain = Strain()
    stress = Stress()

    for Ipre in Rpre
      strain.val = strainStressField.val[Ipre,:]
      mult!(stress,stiffnessField.val[Ipre],strain)
      strainStressField.val[Ipre,:] = stress.val
    end
    strainStressField
  end

  function avg{R}(strainStressField::StrainStressField{R})
    sum = Array{R}(6)
    Iter = CartesianRange(size(strainStressField.val)[1:end-1])
    _avg!(sum,strainStressField,Iter)
    sum
  end

  function _avg!(sum,strainStressField,Iter)
    for I in Iter
      sum += strainStressField.val[Iter,:]
    end
    sum = sum/prod(last(Iter).I)
  end
