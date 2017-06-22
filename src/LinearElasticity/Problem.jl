type LinearElasticityProblem <: MacroscopicGradientProblem
  lattice :: Lattice
  stiffness :: CoefficientTensorField
  macroscopicStrain :: Strain
  strain :: Nullable{StrainField{Float64}}
  averageStress :: Nullable{Stress}

  function LinearElasticityProblem{S <: StiffnessTensor}(
                                                         lattice, stiffness ::
                                                         CoefficientTensorField{S},
                                                         macroscopicStrain
                                                        )
    @argcheck size(stiffness) == lattice.size
    @argcheck norm(macroscopicStrain.val) != 0.0

    new(lattice,
        stiffness,
        macroscopicStrain,
        Nullable{StrainField{Float64}}(),
        Nullable{Stress}()
       )
  end
end

function initializeProblem!(P :: LinearElasticityProblem)
  P.strain = StrainField(Float64,P.lattice.size,0.0)
  P.averageStress = Stress(Float64)
  P
end

function unpackProblem(P :: LinearElasticityProblem)
  stress = StressField(Float64,get(P.strain).val)
  strainFourier = StrainField(Complex128,size(get(P.strain)))
  stressFourier = StressField(Complex128,strainFourier.val)

  (get(P.strain),
   stress,
   strainFourier,
   stressFourier,
   P.stiffness,
   P.macroscopicStrain)
end

type LinearElasticityHomogenizationProblem <: EffectiveTensorProblem
  lattice :: Lattice
  stiffness :: CoefficientTensorField
  effectiveStiffness :: Nullable{AnisotropicStiffnessTensor}

  function LinearElasticityHomogenizationProblem{S <: StiffnessTensor}(lattice,
                                                                       stiffness
                                                                       ::
                                                                       CoefficientTensorField{S}
                                                                      )
    @argcheck size(stiffness) == lattice.size

    new(lattice,
        stiffness,
        Nullable{AnisotropicStiffnessTensor}()
       )
  end
end

function initializeProblem!(P :: LinearElasticityHomogenizationProblem)
  P.effectiveStiffness = AnisotropicStiffnessTensor(zeros(Float64,6,6))
  P
end

length(P :: LinearElasticityHomogenizationProblem) = 6

function Base.getindex(P :: LinearElasticityHomogenizationProblem, i :: Int)
  @argcheck 1 <= i <= 6
  vector = zeros(Float64,6)
  vector[i] = 1.0
  macroscopicStrain = Strain(vector)
  LinearElasticityProblem(P.lattice,
                          P.stiffness,
                          macroscopicStrain
                         )
end

function Base.setindex!(P :: LinearElasticityHomogenizationProblem,
                        S :: LinearElasticityProblem,
                        i :: Int)
  @argcheck 1 <= i <= 6
  if isnull(P.effectiveStiffness)
    initializeProblem!(P)
  end
  P.effectiveStiffness[:,i] = get(S.averageStress).val
  P
end

