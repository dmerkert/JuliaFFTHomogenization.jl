type Strain{R <: Number} <: GradientSolutionTensor{R}
  val :: Array{R,1}

  function Strain{R}(val :: Array{R,1})
    @argcheck length(val) == 6
    new(val)
  end
end
Strain{R}(val :: Array{R,1}) = Strain{R}(val)
Strain(R) =  Strain{R}(Array(R,6))
Strain() = Strain(Float64)

function Strain{R <: Number}(val :: Array{R,2})
  @argcheck size(val) == (3,3)
  @argcheck issymmetric(val)

  Strain([val[[1,5,9]];2*val[[8,7,4]]])
end


type Stress{R <: Number} <: FluxSolutionTensor{R}
  val :: Array{R,1}

  function Stress{R}(val :: Array{R,1})
    @argcheck length(val) == 6
    new(val)
  end
end
Stress{R}(val :: Array{R,1}) = Stress{R}(val)
Stress(R) = Stress(Array(R,6))
Stress() = Stress(Float64)

function Stress{R <: Number}(val :: Array{R,2})
  @argcheck size(val) == (3,3)
  @argcheck issymmetric(val)

  Stress([val[[1,5,9]];val[[8,7,4]]])
end

type Displacement{R <: Number} <: PrimarySolutionTensor{R}
  val :: Array{R,1}

  function Displacement{R}(val :: Array{R,1})
    @argcheck length(val) == 3
    new(val)
  end
end
Displacement{R}(val :: Array{R,1}) = Displacement{R}(val)
Displacement(R) = Displacement{R}(Array(R,3))
Displacement() = Displacement(Float64)

for field in [:Strain,:Stress,:Displacement]
  for op in [:+,:-]
    @eval ($op)(A :: ($field), B :: ($field)) = $(field)(($op)(A.val,B.val))
  end

  @eval *{R <: Number}(a :: R, A :: ($field)) = $(field)(a.*A.val)
  @eval /{R <: Number}(A :: ($field), a :: R) = $(field)(A.val./a)
end

function Base.getindex(strain::Strain, i::Int)
  strain.val[i]
end

function Base.getindex(strain::Strain, I::Vararg{Int, 2})
  @argcheck 1 <= I[1] <= 3
  @argcheck 1 <= I[2] <= 3
  I[1] == I[2] && return strain[I[1]]
  I[1] == 1 && I[2] == 2 && return 0.5strain[6]
  I[1] == 1 && I[2] == 3 && return 0.5strain[5]
  I[1] == 2 && I[2] == 1 && return 0.5strain[6]
  I[1] == 2 && I[2] == 3 && return 0.5strain[4]
  I[1] == 3 && I[2] == 1 && return 0.5strain[5]
  I[1] == 3 && I[2] == 2 && return 0.5strain[4]
end

function Base.setindex!{R}(strain::Strain{R}, v::R, i::Int)
  strain.val[i] = v
  strain
end

function Base.setindex!{R}(strain::Strain{R}, v::R, I::Vararg{Int, 2})
  @argcheck 1 <= I[1] <= 3
  @argcheck 1 <= I[2] <= 3
  I[1] == I[2] && (strain[I[1]] = v)
  I[1] == 1 && I[2] == 2 && (strain[6] = 2.0v)
  I[1] == 1 && I[2] == 3 && (strain[5] = 2.0v)
  I[1] == 2 && I[2] == 1 && (strain[6] = 2.0v)
  I[1] == 2 && I[2] == 3 && (strain[4] = 2.0v)
  I[1] == 3 && I[2] == 1 && (strain[5] = 2.0v)
  I[1] == 3 && I[2] == 2 && (strain[4] = 2.0v)
  strain
end

function Base.getindex(stress::Stress, i::Int)
  stress.val[i]
end

function Base.getindex(stress::Stress, I::Vararg{Int, 2})
  @argcheck 1 <= I[1] <= 3
  @argcheck 1 <= I[2] <= 3
  I[1] == I[2] && return stress[I[1]]
  I[1] == 1 && I[2] == 2 && return stress[6]
  I[1] == 1 && I[2] == 3 && return stress[5]
  I[1] == 2 && I[2] == 1 && return stress[6]
  I[1] == 2 && I[2] == 3 && return stress[4]
  I[1] == 3 && I[2] == 1 && return stress[5]
  I[1] == 3 && I[2] == 2 && return stress[4]
end

function Base.setindex!{R}(stress::Stress{R}, v::R, i::Int)
  stress.val[i] = v
  stress
end

function Base.setindex!{R}(stress::Stress{R}, v::R, I::Vararg{Int, 2})
  @argcheck 1 <= I[1] <= 3
  @argcheck 1 <= I[2] <= 3
  I[1] == I[2] && (stress[I[1]] = v)
  I[1] == 1 && I[2] == 2 && (stress[6] = v)
  I[1] == 1 && I[2] == 3 && (stress[5] = v)
  I[1] == 2 && I[2] == 1 && (stress[6] = v)
  I[1] == 2 && I[2] == 3 && (stress[4] = v)
  I[1] == 3 && I[2] == 1 && (stress[5] = v)
  I[1] == 3 && I[2] == 2 && (stress[4] = v)
  stress
end

function Base.getindex(dispacement :: Displacement, i::Int)
  Displacement.val[i]
end

function Base.setindex!{R}(displacement :: Displacement{R}, v::R, i::Int)
  displacement.val[i] = v
  displacement
end

function fromVoigt!{R <: Number}(A :: Array{R,2},strain::Strain{R})
  @argcheck size(A) == (3,3)

  A[1,1] = strain[1]
  A[2,2] = strain[2]
  A[3,3] = strain[3]
  A[2,3] = 0.5strain[4]
  A[1,3] = 0.5strain[5]
  A[1,2] = 0.5strain[6]
  A[2,1] = A[1,2]
  A[3,1] = A[1,3]
  A[3,2] = A[2,3]
  A
end

function toVoigt!{R <: Number}(strain :: Strain{R}, A :: Array{R,2})
  @argcheck size(A) == (3,3)

  strain[1] = A[1,1]
  strain[2] = A[2,2]
  strain[3] = A[3,3]
  strain[4] = 2.0A[2,3]
  strain[5] = 2.0A[1,3]
  strain[6] = 2.0A[1,2]
  strain
end

function fromVoigt!{R <: Number}(A :: Array{R,2},stress::Stress{R})
  @argcheck size(A) == (3,3)

  A[1,1] = stress[1]
  A[2,2] = stress[2]
  A[3,3] = stress[3]
  A[2,3] = stress[4]
  A[1,3] = stress[5]
  A[1,2] = stress[6]
  A[2,1] = A[1,2]
  A[3,1] = A[1,3]
  A[3,2] = A[2,3]
  A
end

function toVoigt!{R <: Number}(stress :: Stress{R}, A :: Array{R,2})
  @argcheck size(A) == (3,3)

  stress[1] = A[1,1]
  stress[2] = A[2,2]
  stress[3] = A[3,3]
  stress[4] = A[2,3]
  stress[5] = A[1,3]
  stress[6] = A[1,2]
  stress
end

