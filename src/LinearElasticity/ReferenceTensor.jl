export getReferenceTensor

function getReferenceTensor(
                            stiffness ::
                            CoefficientTensorField{T,N},
                            solver :: BasicScheme
                           ) where {T <: IsotropicStiffnessTensor,N}

  (minMu ,maxMu ,minLambda ,maxLambda) = 
  _getReferenceTensorMinMax!(
                             stiffness,
                             solver,
                             CartesianRange(size(stiffness))
                            )

  IsotropicStiffnessTensor(0.5(minLambda+maxLambda), 0.5(minMu+maxMu))
end

function _getReferenceTensorMinMax!(
                                    stiffness ::
                                    CoefficientTensorField{T,N},
                                    solver :: BasicScheme,
                                    R :: CartesianRange{CartesianIndex{N}}
                                   ) where {T <: IsotropicStiffnessTensor,N}

  minMu = typemax(Float64)
  maxMu = typemin(Float64)
  minLambda = typemax(Float64)
  maxLambda = typemin(Float64)
  for i in R
    minMu = min(minMu,stiffness[i].mu.val::Float64)
    maxMu = max(maxMu,stiffness[i].mu.val::Float64)
    minLambda = min(minLambda,stiffness[i].lambda.val::Float64)
    maxLambda = max(maxLambda,stiffness[i].lambda.val::Float64)
  end
  (minMu,maxMu,minLambda,maxLambda)
end

function getReferenceTensor(
                            stiffness :: CoefficientTensorField{S,N},
                            solver :: BasicScheme
                           ) where {S <: StiffnessTensor, N}

  (minOfMinEig,maxOfMinEig,minOfMaxEig,maxOfMaxEig) =
  _eigBounds(stiffness,
             CartesianRange(size(stiffness))
            )

  #@show mu = ShearModulus((minOfMinEig + maxOfMinEig)/4.0)
  mu = ShearModulus(0.5*0.5(minOfMinEig + maxOfMaxEig))
  lambda = LamesFirstParameter(0.0)

  IsotropicStiffnessTensor(lambda,mu)
end

function _eigBounds(
                    stiffness :: CoefficientTensorField{S,N},
                    R :: CartesianRange{CartesianIndex{N}}
                   ) where {S <: StiffnessTensor, N}

  minOfMinEig = typemax(Float64)
  maxOfMinEig = typemin(Float64)
  minOfMaxEig = typemax(Float64)
  maxOfMaxEig = typemin(Float64)

  for i in R
    eigs = eig(stiffness[i])::Array{Float64}
    minE = min(eigs...)::Float64
    maxE = max(eigs...)::Float64

    minOfMinEig = min(minOfMinEig,minE)
    maxOfMinEig = max(maxOfMinEig,minE)
    minOfMaxEig = min(minOfMaxEig,maxE)
    maxOfMaxEig = max(maxOfMaxEig,maxE)
  end
  (minOfMinEig,maxOfMinEig,minOfMaxEig,maxOfMaxEig)
end
