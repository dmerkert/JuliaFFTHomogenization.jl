module JuliaBasicScheme

import Base.*
import Base./
import Base.+
import Base.-
import Base.getindex
import Base.norm
import Base.promote_rule
import Base.convert
using ArgCheck
using MPAWL
import Base.size
using Polynomials
using IterationManagers

FFTW.set_num_threads(Sys.CPU_CORES)

export
  #Types
  CoefficientTensor,
  SolutionTensor,
  GradientSolutionTensor,
  FluxSolutionTensor,
  PrimarySolutionTensor,
  GreenOperator,
  Transformation,
  Solver,
  Problem,
  MacroscopicGradientProblem,
  EffectiveTensorProblem,
  SolutionTensorField,
  CoefficientTensorField,
  size,
  mult!,
  Pattern,
  #StiffnessParameters
  StiffnessParameter,
  BulkModulus,
  YoungsModulus,
  LamesFirstParameter,
  ShearModulus,
  PoissonsRatio,
  convert,
  #StrainStressDisplacement
  Strain,
  Stress,
  Displacement,
  fromVoigt!,
  toVoigt!,
  #StrainStressDisplacementFields
  StrainField,
  StressField,
  DisplacementField,
  init!,
  #StiffnessTensors
  StiffnessTensor,
  IsotropicStiffnessTensor,
  TransversalIsotropicZStiffnessTensor,
  DiagonalStiffnessTensor,
  AnisotropicStiffnessTensor,
  #GreenOperators
  GreenOperator,
  Gamma0,
  #Transformations
  FFT,
  #ApproximationMethod
  ApproximationMethod,
  #Problem
  LinearElasticityProblem,
  LinearElasticityHomogenizationProblem,
  #Solvers
  BasicScheme,
  solve!


  



include("Types.jl")
include("StiffnessParameters.jl")
include("StrainStressDisplacement.jl")
include("StrainStressDisplacementFields.jl")
include("StiffnessTensors.jl")
include("GreenOperators.jl")
include("Transformations.jl")
include("ApproximationMethod.jl")
include("Problem.jl")
include("Solvers.jl")
include("ReferenceTensor.jl")




#include("TypesLinearElasticity.jl")
#include("SolutionFieldArithmetic.jl")
#include("LinearElasticity.jl")
#include("Misc.jl")
#include("Elasticity_Hashin_ellipsoid.jl")
#include("Ellipsoids.jl")

end # module
