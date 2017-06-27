StrainField{R} = SolutionTensorField{R,Strain}
StressField{R} = SolutionTensorField{R,Stress}
DisplacementField{R} = SolutionTensorField{R,Displacement}

StrainField{R}(dims :: Tuple) where {R} =
SolutionTensorField{R,Strain}(Array{R}((dims...,6)))

StrainField{R}(dims :: Tuple,val :: R) where {R} =
SolutionTensorField{R,Strain}(fill(val,(dims...,6)))

StressField{R}(dims :: Tuple) where {R} =
SolutionTensorField{R,Stress}(Array{R}((dims...,6)))

StressField{R}(dims :: Tuple,val :: R) where {R} =
SolutionTensorField{R,Stress}(fill(val,(dims...,6)))

DisplacementField{R}(dims :: Tuple) where {R} =
SolutionTensorField{R,Displacement}(Array{R}((dims...,3)))

DisplacementField{R}(dims :: Tuple,val :: R) where {R} =
SolutionTensorField{R,Displacement}(fill(val,(dims...,3)))

for f in [:StrainField,:StressField,:DisplacementField]
    @eval size(field :: ($f)) = size(field.val)[1:(end-1)]

    @eval function init!(field :: ($f), tensor)
        #for i in CartesianRange(size(field))
        #    field[i] = tensor
        #end
        for i in 1:6
            field[CartesianRange(size(field)),i] = tensor.val[i]
        end
        field
    end

    @eval ($f){R}(field :: F) where {R,F <: SolutionTensorField} =
    ($f){R}(field.val)
end

function average{R}(field :: StrainField{R})
    strain = Strain{R}()
    strain.val = zeros(R,6)

    for i in CartesianRange(size(field))
        strain += field[i]
    end
    strain.val = strain.val/prod(size(field))
    strain
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
    field.val[I...,:] = tensor.val
    field
end

function setindex!(field :: StrainField{R},
                   tensor :: Strain{R},
                   I::CartesianIndex) where {R}
    field.val[I,:] = tensor.val
    field
end


function setindex!(field::StressField{R},
                   tensor::Stress{R},
                   I::Vararg{Int}) where {R}
    field.val[I...,:] = tensor.val
    field
end

function setindex!(field :: StressField{R},
                   tensor :: Stress{R},
                   I::CartesianIndex) where {R}
    field[I,:] = tensor.val
    field
end

function setindex!(field::DisplacementField{R},
                           tensor::Displacement{R},
                           I::Vararg{Int}) where {R}
    field.val[I...,:] = tensor.val
    field
end
function setindex!(field :: DisplacementField{R},
                   tensor :: Displacement{R},
                   I::CartesianIndex) where {R}
    field[I,:] = tensor.val
    field
end


Base.getindex(field::StrainField, I::Vararg{Int}) = Strain(field.val[I...,:])
Base.getindex(field::StrainField, I::CartesianIndex) = field[I.I...]
Base.getindex(field::StressField, I::Vararg{Int}) = Stress(field.val[I...,:])
Base.getindex(field::StressField, I::CartesianIndex) = field[I.I...]
Base.getindex(field::DisplacementField, I::Vararg{Int}) = Displacement(field.val[I...,:])
Base.getindex(field::DisplacementField, I::CartesianIndex) = field[I.I...]

