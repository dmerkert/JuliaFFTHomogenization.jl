function Transform!{R}(solutionTensor :: SolutionTensor, A :: Array{R,2})
  @argcheck size(A) == (3,3)
  @argcheck typeof(solutionTensor.val[1]) == R

  solutionTensorMatrix = zeros(A)

  fromVoigt!(solutionTensorMatrix,solutionTensor)
  solutionTensorMatrix = A*solutionTensorMatrix*A'
  toVoigt!(solutionTensor,solutionTensorMatrix)
end
