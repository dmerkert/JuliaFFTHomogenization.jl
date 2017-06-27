module JuliaFFTHomogenization

import Base.*
import Base./
import Base.+
import Base.-
import Base.getindex
import Base.norm
import Base.promote_rule
import Base.convert
import Base.length
import Base.copy
import Base.copy!
import Base.setindex!
import Base.getindex
using ArgCheck
using MPAWL
import Base.size
using Polynomials

FFTW.set_num_threads(Sys.CPU_CORES)

export
  #MPAWL
  Lattice,
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
  average,
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
  printNumericalError,
  #Solvers
  BasicScheme,
  solve!,
  #referenceTensor
  getReferenceTensor,
  #ElasticityHashinEllipsoid
  ElasticityHashinEllipsoid,
  ElasticityHashinEllipsoidAnalytic




include("Base/Types.jl")
include("Base/Transformations.jl")
include("Base/ApproximationMethod.jl")
include("Base/Solvers.jl")

include("LinearElasticity/StiffnessParameters.jl")
include("LinearElasticity/StrainStressDisplacement.jl")
include("LinearElasticity/StrainStressDisplacementFields.jl")
include("LinearElasticity/StiffnessTensors.jl")
include("LinearElasticity/GreenOperators.jl")
include("LinearElasticity/Problem.jl")
include("LinearElasticity/ReferenceTensor.jl")

include("Misc/Misc.jl")
include("Misc/Ellipsoids.jl")

include("ExampleGeometries/ElasticityHashinEllipsoid.jl")

end # module
