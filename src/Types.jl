abstract CoefficientTensor
abstract SolutionTensor{R <: Number}

abstract GradientSolutionTensor{R} <: SolutionTensor{R}
abstract FluxSolutionTensor{R} <: SolutionTensor{R}
abstract PrimarySolutionTensor{R} <: SolutionTensor{R}

abstract GreenOperator
abstract Transformation
abstract Solver
abstract Problem
abstract MacroscopicGradientProblem <: Problem
abstract EffectiveTensorProblem <: Problem

type SolutionTensorField{R <: Number, F <: SolutionTensor}
  val :: Array{R}

  function SolutionTensorField(r :: Type{R},f :: Type{F},dims)
    new(Array{R}(dims))
  end

  function SolutionTensorField(r :: Type{R},f :: Type{F},dims,val :: R)
    new(val*ones(R,dims))
  end
end
SolutionTensorField{R,F}(r :: Type{R},f :: Type{F},dims) =
SolutionTensorField{R,F}(r,f,dims)
SolutionTensorField{R,F}(r :: Type{R},f :: Type{F},dims,val :: R) =
SolutionTensorField{R,F}(r,f,dims,val)

type CoefficientTensorField{T <: CoefficientTensor}
  val :: Array{T}

  function CoefficientTensorField(t :: Type{T}, dims)
    new(Array{T}(dims))
  end
end

size(A :: CoefficientTensorField) = size(A.val)

function Base.setindex!(field::CoefficientTensorField,
                           tensor::CoefficientTensor,
                           I::Vararg{Int})
    field.val[I...] = tensor
    field
end
Base.setindex!(field :: CoefficientTensorField,
               tensor :: CoefficientTensor,
               I::CartesianIndex) = (field[I.I...] = tensor)

Base.getindex(field::CoefficientTensorField, I::Vararg{Int}) = field.val[I...]
Base.getindex(field::CoefficientTensorField, I::CartesianIndex) = field[I.I...]

for op in [:+,:-]
  @eval ($op){R, F}(A :: SolutionTensorField{R,F},
                    B :: SolutionTensorField{R,F}) = 
  SolutionTensorField(R,F,($op)(A.val,B.val))

  @eval ($op){R,S,F}(A :: SolutionTensorField{R,F},
                     B :: SolutionTensorField{S,F}) =
  SolutionTensorField(promote_type(R,S),F,($op)(A.val,B.val))

  @eval function ($op){T,S}(A :: CoefficientTensorField{T},
                            B :: CoefficientTensorField{S})
    @argcheck size(A) == size(B)
    C = CoefficientTensorField(promote_type(T,S),size(A))
    for i in CartesianRange(ones(size(A)),size(A))
      C[i] = ($op)(A[i],B[i])
    end
  end
end

function mult!(c :: SolutionTensorField,
               A :: CoefficientTensorField,
               b :: SolutionTensorField)
  @argcheck size(A) == size(b) == size(c)

  for i in CartesianRange(ones(size(A)),size(A))
    c[i] = mult!(c[i],A[i],b[i])
  end
  c
end

function mult!(c :: SolutionTensorField,
               G :: GreenOperator,
               b :: SolutionTensorField,
               referenceCoefficient :: CoefficientTensor,
               L :: Lattice
  )
  @argcheck size(b) == size(c)

  for i in getFrequencyIterator(L)
    frequency = getFrequencyPoint(M,i)
    c[i] = mult!(c[i],G,b[i],referenceCoefficient,frequency)
  end
  c
end


abstract Pattern
