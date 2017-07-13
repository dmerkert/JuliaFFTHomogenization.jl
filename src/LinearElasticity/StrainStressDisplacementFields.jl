export StrainField,
       StressField,
       DisplacementField,
       init!,
       average



StrainField{R} = SolutionTensorField{R,Strain}
StressField{R} = SolutionTensorField{R,Stress}
DisplacementField{R} = SolutionTensorField{R,Displacement}

StrainField{R}(dims :: Tuple) where {R} =
SolutionTensorField{R,Strain}(Array{R}((6,dims...)))

StrainField{R}(dims :: Tuple,val :: R) where {R} =
SolutionTensorField{R,Strain}(fill(val,(6,dims...)))

StressField{R}(dims :: Tuple) where {R} =
SolutionTensorField{R,Stress}(Array{R}((6,dims...)))

StressField{R}(dims :: Tuple,val :: R) where {R} =
SolutionTensorField{R,Stress}(fill(val,(6,dims...)))

DisplacementField{R}(dims :: Tuple) where {R} =
SolutionTensorField{R,Displacement}(Array{R}((3,dims...)))

DisplacementField{R}(dims :: Tuple,val :: R) where {R} =
SolutionTensorField{R,Displacement}(fill(val,(3,dims...)))

for f in [:StrainField,:StressField,:DisplacementField]
    @eval size(field :: ($f)) = size(field.val)[2:end]

    @eval function init!(field :: ($f), tensor)
        field.val .= tensor.val
    end

    @eval ($f){R}(field :: F) where {R,F <: SolutionTensorField} =
    ($f){R}(field.val)
end

function average{R}(field :: StrainField{R})
    Strain(sum(field.val,2:length(size(field.val)))[:]/prod(size(field)))
end

function average{R}(field :: StressField{R})
    stress = Stress(R)
    stress.val = zeros(R,6)

    for i in CartesianRange(size(field))
        stress += field[i]
    end
    stress.val = stress.val/prod(size(field))
    stress
end

function setindex!(field::StrainField{R},
                   tensor::Strain{R},
                   I::Vararg{Int}) where {R}
    field.val[:,I...] = tensor.val
    field
end
setindex!(field :: StrainField{R},
          tensor :: Strain{R},
          I::CartesianIndex) where {R} = (field[I.I...] = tensor)

function setindex!(field::StressField{R},
                           tensor::Stress{R},
                           I::Vararg{Int}) where {R}
    field.val[:,I...] = tensor.val
    field
end
setindex!(field :: StressField{R},
          tensor :: Stress{R},
          I::CartesianIndex) where {R} = (field[I.I...] = tensor)

function setindex!(field::DisplacementField{R},
                   tensor::Displacement{R},
                   I::Vararg{Int}) where {R}
    field.val[:,I...] = tensor.val
    field
end
setindex!(field :: DisplacementField{R},
          tensor :: Displacement{R},
          I::CartesianIndex) where {R} = (field[I.I...] = tensor)


getindex(field::StrainField, I::Vararg{Int}) = Strain(field.val[:,I...])
getindex(field::StrainField, I::CartesianIndex) = field[I.I...]
getindex(field::StressField, I::Vararg{Int}) = Stress(field.val[:,I...])
getindex(field::StressField, I::CartesianIndex) = field[I.I...]
getindex(field::DisplacementField, I::Vararg{Int}) = Displacement(field.val[:,I...])
getindex(field::DisplacementField, I::CartesianIndex) = field[I.I...]

