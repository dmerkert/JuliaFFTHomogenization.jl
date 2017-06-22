function getReferenceTensor{S <: IsotropicStiffnessTensor}(stiffness ::
                                                           CoefficientTensorField{S},
                                                           solver :: BasicScheme)
  minMu = Inf
  maxMu = -Inf
  minLambda = Inf
  maxLambda = -Inf
  for i in CartesianRange(size(stiffness))
    minMu = min(minMu,stiffness[i].mu)
    maxMu = max(maxMu,stiffness[i].mu)
    minLambda = min(minLambda,stiffness[i].lambda)
    maxLambda = max(maxLambda,stiffness[i].lambda)
  end

  IsotropicStiffnessTensor(0.5(minLambda+maxLambda), 0.5(minMu+maxMu))
end

function getReferenceTensor{S <: StiffnessTensor}(stiffness ::
                                                  CoefficientTensorField{S},
                                                  solver :: BasicScheme)

  (minOfMinEig,maxOfMinEig,minOfMaxEig,maxOfMaxEig) = _eigBounds(stiffness)

  #@show mu = ShearModulus((minOfMinEig + maxOfMinEig)/4.0)
  mu = ShearModulus(0.5*0.5(minOfMinEig + maxOfMaxEig))
  lambda = LamesFirstParameter(0.0)

  IsotropicStiffnessTensor(lambda,mu)
end

function _eigBounds{S <: StiffnessTensor}(stiffness :: CoefficientTensorField{S})
  minOfMinEig = Inf
  maxOfMinEig = -Inf
  minOfMaxEig = Inf
  maxOfMaxEig = -Inf

  for i in CartesianRange(size(stiffness))
    eigs = eig(stiffness[i])
    minE = min(eigs...)
    maxE = max(eigs...)

    (minOfMinEig > minE) && (minOfMinEig = minE)
    (maxOfMinEig < minE) && (maxOfMinEig = minE)
    (minOfMaxEig > maxE) && (minOfMaxEig = maxE)
    (maxOfMaxEig < maxE) && (maxOfMaxEig = maxE)
  end
  (minOfMinEig,maxOfMinEig,minOfMaxEig,maxOfMaxEig)
end
