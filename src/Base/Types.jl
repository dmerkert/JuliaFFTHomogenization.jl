export CoefficientTensor,
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
       Pattern,
       SolutionTensorField,
       CoefficientTensorField,
       copy!,
       copy,
       norm,
       mult!

abstract type CoefficientTensor end
abstract type SolutionTensor{R <: Number} end

abstract type GradientSolutionTensor{R} <: SolutionTensor{R} end
abstract type FluxSolutionTensor{R} <: SolutionTensor{R} end
abstract type PrimarySolutionTensor{R} <: SolutionTensor{R} end

abstract type GreenOperator end
abstract type Transformation end
abstract type Solver end
abstract type Problem end
abstract type MacroscopicGradientProblem <: Problem end
abstract type EffectiveTensorProblem <: Problem end

abstract type Pattern end

type SolutionTensorField{R <: Number, F <: SolutionTensor}
  val :: Array{R}

  SolutionTensorField{R,F}(val :: Array{R}) where {R,F} = new{R,F}(val)
end

function copy!{R,F}(dest :: SolutionTensorField{R,F},
                    source :: SolutionTensorField{R,F})
  copy!(dest.val,source.val)
  dest
end

function copy{R,F}(source :: SolutionTensorField{R,F})
  dest = SolutionTensorField{R,F}(similar(source.val))
  copy!(dest,source)
end

function norm(S :: SolutionTensorField)
  norm(S.val[:])
end

function norm(S :: SolutionTensor)
  norm(S.val[:])
end


type CoefficientTensorField{T <: CoefficientTensor}
  val :: Array{T}

  CoefficientTensorField{T}(dims :: Tuple) where {T} = new{T}(Array{T}(dims))
  CoefficientTensorField{T}(val :: Array{T}) where {T} = new{T}(val)
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
  SolutionTensorField{R,F}(($op).(A.val,B.val))

  @eval ($op){R,S,F}(A :: SolutionTensorField{R,F},
                     B :: SolutionTensorField{S,F}) =
  SolutionTensorField{promote_type(R,S),F}(($op).(A.val,B.val))

  @eval function ($op){T,S}(A :: CoefficientTensorField{T},
                            B :: CoefficientTensorField{S})
    @argcheck size(A) == size(B)
    C = CoefficientTensorField{promote_type(T,S)}(size(A))
    for i in CartesianRange(size(A))
      C[i] = ($op)(A[i],B[i])
    end
    C
  end

  @eval function ($op){T,S <: CoefficientTensor}(A :: CoefficientTensorField{T},
                                                 B :: S)
    C = CoefficientTensorField{promote_type(T,S)}(size(A))
    for i in CartesianRange(size(A))
      C[i] = ($op)(A[i],B)
    end
    C
  end
end

function mult!(c :: SolutionTensorField,
               A :: CoefficientTensorField,
               b :: SolutionTensorField)
  @argcheck size(A) == size(b)
  @argcheck size(b) == size(c)

  for i in CartesianRange(size(A))
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
    frequency = getFrequencyPoint(L,i)
    c[i] = mult!(c[i],G,b[i],referenceCoefficient,2.0pi*frequency)
  end
  c
end


