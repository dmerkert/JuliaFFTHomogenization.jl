export StrainField,
       StressField,
       DisplacementField,
       average



StrainField{R,N} = SolutionTensorField{R,Strain,N}
StressField{R,N} = SolutionTensorField{R,Stress,N}
DisplacementField{R,N} = SolutionTensorField{R,Displacement,N}

StrainField{R}(dims :: NTuple{N,Int}) where {R,N} =
SolutionTensorField{R,Strain}(Array{R}((6,dims...)))

StrainField{R}(dims :: NTuple{N,Int},val :: R) where {R,N} =
SolutionTensorField{R,Strain}(fill(val,(6,dims...)))

StressField{R}(dims :: NTuple{N,Int}) where {R,N} =
SolutionTensorField{R,Stress}(Array{R}((6,dims...)))

StressField{R}(dims :: NTuple{N,Int},val :: R) where {R,N} =
SolutionTensorField{R,Stress}(fill(val,(6,dims...)))

DisplacementField{R}(dims :: NTuple{N,Int}) where {R,N} =
SolutionTensorField{R,Displacement}(Array{R}((3,dims...)))

DisplacementField{R}(dims :: NTuple{N,Int},val :: R) where {R,N} =
SolutionTensorField{R,Displacement}(fill(val,(3,dims...)))

for f in [:StrainField,:StressField,:DisplacementField]
    @eval size(field :: ($f){R,N}) where {R,N} = size(field.val)[2:end]

    @eval ($f){R}(field :: F) where {R,F <: SolutionTensorField} =
    ($f){R}(field.val)
end

function average(field :: StrainField{R,N}) where {R,N}
    Strain(sum(field.val,2:length(size(field.val)))[:]/prod(size(field)))
end

function average(field :: StressField{R,N}) where {R,N}
    Stress(sum(field.val,2:length(size(field.val)))[:]/prod(size(field)))
end

function setindex!(field::StrainField{R,N},
                   tensor::Strain{R},
                   I::Vararg{Int}) where {R,N}
    field.val[:,I...] = tensor.val
    field
end
setindex!(field :: StrainField{R,N},
          tensor :: Strain{R},
          I::CartesianIndex) where {R,N} = (field[I.I...] = tensor)

function setindex!(field::StressField{R,N},
                           tensor::Stress{R},
                           I::Vararg{Int}) where {R,N}
    field.val[:,I...] = tensor.val
    field
end
setindex!(field :: StressField{R,N},
          tensor :: Stress{R},
          I::CartesianIndex) where {R,N} = (field[I.I...] = tensor)

function setindex!(field::DisplacementField{R,N},
                   tensor::Displacement{R},
                   I::Vararg{Int}) where {R,N}
    field.val[:,I...] = tensor.val
    field
end
setindex!(field :: DisplacementField{R,N},
          tensor :: Displacement{R},
          I::CartesianIndex) where {R,N} = (field[I.I...] = tensor)


getindex(field::StrainField{R,N}, I::Vararg{Int}) where {R,N} = Strain(field.val[:,I...])
getindex(field::StrainField{R,N}, I::CartesianIndex) where {R,N} = field[I.I...]
getindex(field::StressField{R,N}, I::Vararg{Int}) where {R,N} = Stress(field.val[:,I...])
getindex(field::StressField{R,N}, I::CartesianIndex) where {R,N} = field[I.I...]
getindex(field::DisplacementField{R,N}, I::Vararg{Int}) where {R,N} = Displacement(field.val[:,I...])
getindex(field::DisplacementField{R,N}, I::CartesianIndex) where {R,N} = field[I.I...]

