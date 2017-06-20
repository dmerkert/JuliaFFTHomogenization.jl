function getReferenceTensor{S <: IsotropicStiffnessTensor}(stiffness ::
                                                           CoefficientTensorField{S},
                                                           solver :: BasicScheme)
  mu = stiffness[:].mu
  @show mu

end
