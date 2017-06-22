typealias StrainField{R} SolutionTensorField{R,Strain}
typealias StressField{R} SolutionTensorField{R,Stress}
typealias DisplacementField{R} SolutionTensorField{R,Displacement}

StrainField{R}(r :: Type{R},dims) =
SolutionTensorField{R,Strain}(R,Strain,(dims...,6))

StrainField{R}(r :: Type{R},dims,val) =
SolutionTensorField{R,Strain}(R,Strain,(dims...,6),val)

StrainField{R}(r :: Type{R},val :: Array) =
SolutionTensorField{R,Strain}(R,Strain,val)

StrainField{R}(r :: Type{R},val :: Array) =
SolutionTensorField{R,Strain}(R,Strain,val)

StressField{R}(r :: Type{R},dims) =
SolutionTensorField{R,Stress}(R,Stress,(dims...,6))

StressField{R}(r :: Type{R},dims,val) =
SolutionTensorField{R,Stress}(R,Stress,(dims...,6),val)

StressField{R}(r :: Type{R},val :: Array) =
SolutionTensorField{R,Stress}(R,Stress,val)

DisplacementField{R}(r :: Type{R},dims) =
SolutionTensorField{R,Displacement}(R,Displacement,(dims...,3))

DisplacementField{R}(r :: Type{R},dims,val) =
SolutionTensorField{R,Displacement}(R,Displacement,(dims...,3),val)

DisplacementField{R}(r :: Type{R},val :: Array) =
SolutionTensorField{R,Displacement}(R,Displacement,val)

for f in [:StrainField,:StressField,:DisplacementField]
    @eval size(field :: ($f)) = size(field.val)[1:(end-1)]

    @eval function init!(field :: ($f), tensor)
        for i in CartesianRange(size(field))
            field[i] = tensor
        end
        field
    end

    @eval ($f){R,F <: SolutionTensorField}(r :: Type{R}, field :: F) =
    ($f){R}(r,field.val)
end

function average{R}(field :: StressField{R})
    stress = Stress(R)
    stress.val = sum(field.val,size(field.val)[1:(end-1)])
    stress
end

function Base.setindex!{R}(field::StrainField{R},
                           tensor::Strain{R},
                           I::Vararg{Int})
    field.val[I...,:] = tensor.val
    field
end
Base.setindex!(field :: StrainField,
               tensor :: Strain,
               I::CartesianIndex) = (field[I.I...] = tensor)


function Base.setindex!{R}(field::StressField{R},
                           tensor::Stress{R},
                           I::Vararg{Int})
    field.val[I...,:] = tensor.val
    field
end
Base.setindex!(field :: StressField,
               tensor :: Stress,
               I::CartesianIndex) = (field[I.I...] = tensor)

function Base.setindex!{R}(field::DisplacementField{R},
                           tensor::Displacement{R},
                           I::Vararg{Int})
    field.val[I...,:] = tensor.val
    field
end
Base.setindex!(field :: DisplacementField,
               tensor :: Displacement,
               I::CartesianIndex) = (field[I.I...] = tensor)


Base.getindex(field::StrainField, I::Vararg{Int}) = Strain(field.val[I...,:])
Base.getindex(field::StrainField, I::CartesianIndex) = field[I.I...]
Base.getindex(field::StressField, I::Vararg{Int}) = Stress(field.val[I...,:])
Base.getindex(field::StressField, I::CartesianIndex) = field[I.I...]
Base.getindex(field::DisplacementField, I::Vararg{Int}) = Displacement(field.val[I...,:])
Base.getindex(field::DisplacementField, I::CartesianIndex) = field[I.I...]

