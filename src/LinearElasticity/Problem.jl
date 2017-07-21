export LinearElasticityProblem,
       copy,
       postprocess!,
       printNumericalError,
       initializeProblem!,
       unpackProblem,
       LinearElasticityHomogenizationProblem


type LinearElasticityProblem{N,T,M,R} <: MacroscopicGradientProblem
 lattice :: Lattice
 stiffness :: CoefficientTensorField{T,N}
 macroscopicStrain :: Strain{R}
 strain :: Nullable{StrainField{R,M}}
 averageStress :: Nullable{Stress}

 function LinearElasticityProblem(
                                  lattice :: Lattice,
                                  stiffness :: CoefficientTensorField{T,N},
                                  macroscopicStrain :: Strain{R}
                                 ) where {
                                          T <: StiffnessTensor,
                                          N,
                                          R
                                         }
    @argcheck size(stiffness) == lattice.size
    @argcheck norm(macroscopicStrain.val) != 0.0

    new{N,T,N+1,R}(
                   lattice,
                   stiffness,
                   macroscopicStrain,
                   Nullable{StrainField{R,N+1}}(),
                   Nullable{Stress}()
                  )
  end
end

function copy(
              P :: LinearElasticityProblem{N,T,M,R}
             ) where {N,T,M,R}
  LinearElasticityProblem(P.lattice,P.stiffness,P.macroscopicStrain)
end

function postprocess!(
                      P :: LinearElasticityProblem{N,T,M,R}
                     ) where {N,T,M,R}
  @argcheck !isnull(P.strain)

  averageStress = Stress{R}()
  averageStress.val = zeros(R,6)
  stressTmp = Stress{R}()
  strainTmp = Strain{R}()
  for i in CartesianRange(size(P.stiffness))
    stressTmp = mult!(stressTmp,P.stiffness[i],get(P.strain)[i])
    averageStress += stressTmp
  end
  averageStress.val = averageStress.val/P.lattice.m
  P.averageStress = averageStress
end

function printNumericalError(
                             PAnalytic :: LinearElasticityProblem{N,T,M,R},
                             PNumeric  :: LinearElasticityProblem{N,T,M,R}
                            ) where {N,T,M,R}
  @argcheck PAnalytic.lattice.M == PNumeric.lattice.M
  @argcheck !isnull(PAnalytic.strain)
  @argcheck !isnull(PAnalytic.averageStress)
  @argcheck !isnull(PNumeric.strain)
  @argcheck !isnull(PNumeric.averageStress)

  errorAverageStress =
  norm(get(PAnalytic.averageStress)-get(PNumeric.averageStress))/
  norm(get(PAnalytic.averageStress))

  errorL2Strain =
  norm(get(PAnalytic.strain)-get(PNumeric.strain))/norm(get(PAnalytic.strain))

  print("Error in average stress     : $(errorAverageStress)\n")
  print("l2-error in strain          : $(errorL2Strain)\n")
  (errorAverageStress,errorL2Strain)
end

function initializeProblem!(
                            P :: LinearElasticityProblem{N,T,M,R}
                           ) where {N,T,M,R}
  P.strain = StrainField{R}(P.lattice.size,0.0)
  P.averageStress = Stress{R}()
  P
end

function unpackProblem(
                       P :: LinearElasticityProblem{N,T,M,R}
                      ) where {N,T,M,R}
  stress = StressField{R}(get(P.strain).val)
  strainFourier = StrainField{Complex128}(size(get(P.strain)))
  stressFourier = StressField{Complex128}(strainFourier.val)

  (get(P.strain),
   stress,
   strainFourier,
   stressFourier,
   P.stiffness,
   P.macroscopicStrain)
end

type LinearElasticityHomogenizationProblem{T,N} <: EffectiveTensorProblem
  lattice :: Lattice
  stiffness :: CoefficientTensorField{T,N}
  effectiveStiffness :: Nullable{AnisotropicStiffnessTensor}

  function LinearElasticityHomogenizationProblem(
                                                 lattice :: Lattice,
                                                 stiffness ::
                                                 CoefficientTensorField{T,N}
                                                ) where {T,N}
    @argcheck size(stiffness) == lattice.size

    new{T,N}(lattice,
        stiffness,
        Nullable{AnisotropicStiffnessTensor}()
       )
  end
end

type LinearElasticityHomogenizationProblemWithFields{N,T,M} <: EffectiveTensorProblem
  lattice :: Lattice
  stiffness :: CoefficientTensorField{T,N}
  effectiveStiffness :: Nullable{AnisotropicStiffnessTensor}
  strain :: Array{Nullable{StrainField{Float64,M}},1}


  function LinearElasticityHomogenizationProblemWithFields(
                                                           lattice :: Lattice,
                                                           stiffness
                                                           ::
                                                           CoefficientTensorField{T,N}
                                                          ) where {T,N}
    @argcheck size(stiffness) == lattice.size

    new{N,T,N+1}(lattice,
        stiffness,
        Nullable{AnisotropicStiffnessTensor}(),
        Array{Nullable{StrainField{Float64,N+1}}}(6)
       )
  end
end


function copy(
              P :: LinearElasticityHomogenizationProblem{T,N}
             ) where {T,N}
  LinearElasticityHomogenizationProblem(P.lattice,P.stiffness)
end

function copy(
              P :: LinearElasticityHomogenizationProblemWithFields{T,N}
             ) where {T,N}
  LinearElasticityHomogenizationProblemWithFields(P.lattice,P.stiffness)
end

function initializeProblem!(
                            P :: LinearElasticityHomogenizationProblem{T,N}
                           ) where {T,N}
  P.effectiveStiffness = AnisotropicStiffnessTensor(zeros(Float64,6,6))
  P
end

function initializeProblem!(
                            P ::
                            LinearElasticityHomogenizationProblemWithFields{T,N}
                           ) where {T,N}
  P.effectiveStiffness = AnisotropicStiffnessTensor(zeros(Float64,6,6))
  P
end

length(P :: LinearElasticityHomogenizationProblem{T,N}) where {T,N} = 6
length(P :: LinearElasticityHomogenizationProblemWithFields{T,N}) where {T,N} = 6

function Base.getindex(
                       P :: LinearElasticityHomogenizationProblem{T,N},
                       i :: Int
                      ) where {T,N}
  @argcheck 1 <= i
  @argcheck i <= 6
  vector = zeros(Float64,6)
  vector[i] = 1.0
  macroscopicStrain = Strain(vector)
  LinearElasticityProblem(P.lattice,
                          P.stiffness,
                          macroscopicStrain
                         )
end

function Base.getindex(
                       P :: LinearElasticityHomogenizationProblemWithFields{T,N}
                       , i :: Int
                      ) where {T,N}
  @argcheck 1 <= i
  @argcheck i <= 6
  vector = zeros(Float64,6)
  vector[i] = 1.0
  macroscopicStrain = Strain(vector)
  LinearElasticityProblem(P.lattice,
                          P.stiffness,
                          macroscopicStrain
                         )
end

function Base.setindex!(
                        P :: LinearElasticityHomogenizationProblem{T,N},
                        S :: LinearElasticityProblem{T,N},
                        i :: Int
                       ) where {T,N}
  @argcheck 1 <= i
  @argcheck i <= 6
  if isnull(P.effectiveStiffness)
    initializeProblem!(P)
  end
  get(P.effectiveStiffness).C[:,i] = get(S.averageStress).val
  P
end

function Base.setindex!(
                        P ::
                        LinearElasticityHomogenizationProblemWithFields{T,N},
                        S :: LinearElasticityProblem{T,N},
                        i :: Int
                       ) where {T,N}
  @argcheck 1 <= i
  @argcheck i <= 6
  if isnull(P.effectiveStiffness)
    initializeProblem!(P)
  end
  get(P.effectiveStiffness).C[:,i] = get(S.averageStress).val
  P.strain[i] = get(S.strain)
  P
end


