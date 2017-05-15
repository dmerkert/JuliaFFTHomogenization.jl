module JuliaBasicScheme

import Base.*

export
  #Types
  StiffnessTensor,
  IsotropicStiffnessTensor,
  TransversalIsotropicZ,
  ElasticityVector,
  Strain,
  Stress,
  Displacement,
  StiffnessTensorField,
  StrainStressField,

  #methods
  mult!,
  avg

include("Coefficients.jl")

end # module
