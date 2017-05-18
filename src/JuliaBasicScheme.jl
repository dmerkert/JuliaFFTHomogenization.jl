module JuliaBasicScheme

import Base.*
FFTW.set_num_threads(Sys.CPU_CORES)

export
  #Types Linear Elasticity
  StiffnessTensor,
  BulkModulus,
  YoungsModulus,
  LamesFirstParameter,
  ShearModulus,
  PoissonsRatio,
  Strain,
  Stress,
  Displacement,
  IsotropicStiffnessTensor,
  TransversalIsotropicZ,
  StiffnessTensorField,
  StrainField,
  StressField,
  Displacement,
  #Methods Linear Elasticity
  add!,
  set!,
  get!,
  mult!,
  avg,
  FFT!,
  IFFT!,
  Gamma0


include("Types.jl")
include("TypesLinearElasticity.jl")
include("SolutionFieldArithmetic.jl")
include("LinearElasticity.jl")
include("GreenOperators.jl")

end # module
