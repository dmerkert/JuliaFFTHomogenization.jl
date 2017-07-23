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

type SolutionTensorField{R <: Number, F <: SolutionTensor, N}
  val :: Array{R,N}

  function SolutionTensorField{R,F}(val :: Array{R,N}) where {R,F,N} 
    new{R,F,N}(val)
  end
end

function copy!(
               dest :: SolutionTensorField{R,F,N},
               source :: SolutionTensorField{R,F,N}
              ) where {R,F,N}
  copy!(dest.val,source.val)
  dest
end

function copy(
              source :: SolutionTensorField{R,F,N}
             ) where {R,F,N}
  dest = SolutionTensorField{R,F}(similar(source.val))
  copy!(dest,source)
end

function norm(S :: SolutionTensorField{R,F,N}) where {R,F,N}
  norm(S.val[:])
end

function norm(S :: SolutionTensor{R}) where {R}
  norm(S.val[:])
end

function init!(field :: SolutionTensorField{R,F,N}, tensor :: F) where {N,F,R}
  field.val .= tensor.val
end



type CoefficientTensorField{T <: CoefficientTensor,N}
  val :: Array{T,N}

  CoefficientTensorField{T}(dims :: NTuple{N,Int}) where {T,N} =
  new{T,N}(Array{T,N}(dims))
  CoefficientTensorField(val :: Array{T,N}) where {T,N} = new{T,N}(val)
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
  @eval ($op)(
              A :: SolutionTensorField{R,F,N},
              B :: SolutionTensorField{R,F,N}
             ) where {R,F,N} =
  SolutionTensorField{R,F}(($op).(A.val,B.val))

  @eval ($op)(
              A :: SolutionTensorField{R,F,N},
              B :: SolutionTensorField{S,F,N}
             ) where {R,S,F,N} =
  SolutionTensorField{promote_type(R,S),F}(($op).(A.val,B.val))

  @eval function ($op)(
                       A :: CoefficientTensorField{T,N},
                       B :: CoefficientTensorField{S,N}
                      ) where {T,S,N}
    @argcheck size(A) == size(B)
    C = CoefficientTensorField{promote_type(T,S)}(size(A))
    for i in CartesianRange(size(A))
      C[i] = ($op)(A[i],B[i])
    end
    C
  end

  @eval function ($op)(
                       A :: CoefficientTensorField{T,N},
                       B :: S
                      ) where {T,S <: CoefficientTensor,N}
    C = CoefficientTensorField{promote_type(T,S)}(size(A))
    for i in CartesianRange(size(A))
      C[i] = ($op)(A[i],B)
    end
    C
  end
end

function mult!(
               c :: SolutionTensorField{R,F1,M},
               A :: CoefficientTensorField{T,N},
               b :: SolutionTensorField{R,F2,M}
              ) where {R,F1,F2,T,N,M}
  @argcheck size(A) == size(b)
  @argcheck size(b) == size(c)

  for i in CartesianRange(size(A))
    c[i] = mult!(c[i],A[i],b[i])
  end
  c
end

function mult!(c :: SolutionTensorField{R,F1,M},
               G :: GreenOperator,
               b :: SolutionTensorField{R,F2,M},
               referenceCoefficient :: C,
               L :: Lattice
              ) where {R,F1,F2,M,C<:CoefficientTensor}
  @argcheck size(b) == size(c)

  tmpGamma = initMult(G)
  tmpc = c[start(getFrequencyIterator(L))] :: F1{R}
  tmpb = b[start(getFrequencyIterator(L))] :: F2{R}
  tmpFrequency = getFrequencyPoint(L,start(getFrequencyIterator(L))) :: Array{Int,1}
  tmpFrequencyFloat = zeros(Float64,size(tmpFrequency)) :: Array{Float64,1}

  _mult!(c,G,b,referenceCoefficient,L,tmpGamma,tmpc,tmpb,tmpFrequency,tmpFrequencyFloat,getFrequencyIterator(L))
end

function _mult!(c :: SolutionTensorField{R,F1,M},
               G :: GreenOperator,
               b :: SolutionTensorField{R,F2,M},
               referenceCoefficient :: C,
               L :: Lattice,
               tmpGamma :: TMPG,
               tmpc :: F3,
               tmpb :: F4,
               tmpFrequency :: Array{I,1},
               tmpFrequencyFloat :: Array{Float64,1},
               range :: CartesianRange{CartesianIndex{N}}
              ) where {R,F1,F2,M,C<:CoefficientTensor,TMPG,N,F3,F4,I}

  for i in range
    getFrequencyPoint!(L,i,tmpFrequency,tmpFrequencyFloat)
    tmpb.val .= b.val[:,i]
    mult!(tmpc,G,tmpb,referenceCoefficient,tmpFrequency,tmpGamma)
    c.val[:,i] .= tmpc.val
  end
  c
end
