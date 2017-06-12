module JuliaBasicScheme

import Base.*
import Base.getindex
using ArgCheck
using MPAWL
import Base.size
using Polynomials

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
  Gamma0,
  RotationMatrixNormal,
  DepolarizationFactors,
  fromVoigt!,
  toVoigt!,
  applyGammaHat!,
  Elasticity_Hashin_ellipsoid,
  Elasticity_Hashin_ellipsoid_analytic,
  Lattice


include("Types.jl")
include("TypesLinearElasticity.jl")
include("SolutionFieldArithmetic.jl")
include("LinearElasticity.jl")
include("GreenOperators.jl")
include("Misc.jl")
include("Elasticity_Hashin_ellipsoid.jl")

end # module
