export StrainField,
       StressField,
       DisplacementField,
       average,
       normDifference



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

function normDifference(field1 :: StrainField{R,N1},
                        L1 :: Lattice,
                        field2 :: StrainField{R,N2},
                        L2 :: Lattice,
                        space :: Space
                       ) where {
                                R <: Real,
                                N1,
                                N2,
                                Space <: AnsatzSpace
                               }
    #assumes that L2 is a tensor product lattice
    @argcheck norm(diagm(diag(L2.M))-L2.M) < 1e-12
    #assumes that L1 generates subpattern of L2
    @argcheck isSublattice(L2,L1)


    field1Frequency = StrainField{Complex128}(size(field1))
    field2Frequency = StrainField{Complex128}(size(field2))

    transform!(
               field1Frequency,
               field1,
               space,
               L1
              )

    transform!(
               field2Frequency,
               field2,
               space,
               L2
              )

    _normDifference(field1Frequency,
                    L1,
                    field2Frequency,
                    L2,
                    space,
                    getFrequencyIterator(L2)
                   )
end

function _normDifference(field1Frequency :: StrainField{R,N1},
                         L1 :: Lattice{I1,MF1,MF12},
                        field2Frequency :: StrainField{R,N2},
                        L2 :: Lattice{I2,MF2,MF22},
                        space :: Space,
                        range :: CartesianRange{RR}
                       ) where {
                                R <: Complex,
                                N1,
                                N2,
                                I1,MF1,MF12,
                                I2,MF2,MF22,
                                RR,
                                Space <: AnsatzSpace
                               }


    diffNorm = 0.0
    refNorm = 0.0
    L1Unit = Lattice(L1.M,target="unit")
    L2Unit = Lattice(L2.M,target="unit")
    for coord in range
        freq = getFrequencyPoint(L2,coord)
        h1 = frequencyLatticeBasisDecomp(freq,L1Unit)
        h2 = frequencyLatticeBasisDecomp(freq,L2Unit)
        ck1 = ck(freq,L1,space)/L1.m
        ck2 = ck(freq,L2,space)/L2.m
        refNorm += norm(
                        field2Frequency[(h2+1)...].val*ck2
                       )^2
        diffNorm += norm(
                        field2Frequency[(h2+1)...].val*ck2 -
                        field1Frequency[(h1+1)...].val*ck1
                        )^2
    end

    sqrt(diffNorm)/sqrt(refNorm)
end

