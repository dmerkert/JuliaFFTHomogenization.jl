export subsamplingComposites

function subsamplingComposites(
                               LSuper :: Lattice,
                               LSub :: Lattice,
                               LDecomposition :: Lattice,
                               stiffness :: CoefficientTensorField{T,N},
                               compositeType :: Type{CompositeType}
                              ) where {
                                       T,
                                       N,
                                       CompositeType <: CompositeStiffnessTensor
                              }

  @argcheck isLatticeDecomposition(LSuper,LDecomposition,LSub)
  @argcheck size(stiffness) == LSuper.size

  C = CoefficientTensorField{StiffnessTensor}(LSub.size)

  for coordSub in getSamplingIterator(LSub)
    pointSub = getSamplingPoint(LSub,coordSub)
    CList =
    [
     stiffness[
               CartesianIndex(
                              ((getPatternBasisDecomp(LSuper,pointSub)+1)...)
                             )
              ]
    ]
    volumes = [0.0]

    for coordDecomposition in getSamplingIterator(LDecomposition)
      point = pointSub + LSub.M\getSamplingPoint(LDecomposition,coordDecomposition)
      coordSuper = CartesianIndex(((getPatternBasisDecomp(LSuper,point)+1)...))

      CIndex = findfirst(stiffness[coordSuper] .== CList)
      if CIndex > 0
        volumes[CIndex] += 1.0
      else
        CList = [CList;stiffness[coordSuper]]
        volumes = [volumes;1.0]
      end
    end

    C[coordSub] = compositeType(CList,volumes)
  end

  C
end

function subsamplingComposites(
                               LSuper :: Lattice,
                               LSub :: Lattice,
                               LDecomposition :: Lattice,
                               stiffness :: CoefficientTensorField{T,N},
                               compositeType ::
                               Type{CompositeLaminateStiffnessTensor}
                              ) where {
                                       T,
                                       N
                              }

  @argcheck isLatticeDecomposition(LSuper,LDecomposition,LSub)
  @argcheck size(stiffness) == LSuper.size

  C = CoefficientTensorField{StiffnessTensor}(LSub.size)

  for coordSub in getSamplingIterator(LSub)
    pointSub = getSamplingPoint(LSub,coordSub)
    CList =
    [
     stiffness[
               CartesianIndex(
                              ((getPatternBasisDecomp(LSuper,pointSub)+1)...)
                             )
              ]
    ]
    volumes = [0.0]

    for coordDecomposition in getSamplingIterator(LDecomposition)
      point = pointSub + LSub.M\getSamplingPoint(LDecomposition,coordDecomposition)
      coordSuper = CartesianIndex(((getPatternBasisDecomp(LSuper,point)+1)...))

      CIndex = findfirst(stiffness[coordSuper] .== CList)
      if CIndex > 0
        volumes[CIndex] += 1.0
      else
        CList = [CList;stiffness[coordSuper]]
        volumes = [volumes;1.0]
      end
    end

    if length(CList) > 1

      dominantStiffness = CList[findmax(volumes)[2]]

      centerDominantStiffness = zeros(size(pointSub))
      centerSubvoxel = zeros(size(pointSub))
      nDominantStiffness = 0
      nSubvoxel = 0

      for coordDecomposition in getSamplingIterator(LDecomposition)
        point = pointSub + LSub.M\getSamplingPoint(LDecomposition,coordDecomposition)
        coordSuper = CartesianIndex(((getPatternBasisDecomp(LSuper,point)+1)...))

        centerSubvoxel += point
        nSubvoxel += 1
        if stiffness[coordSuper] == dominantStiffness
          centerDominantStiffness += point
          nDominantStiffness += 1
        end
      end
      normal =
      centerDominantStiffness/nDominantStiffness -
      centerSubvoxel/nSubvoxel

      @assert !any(isnan.(normal))

      if norm(normal) ≈ 0.0
        C[coordSub] = CompositeArithmeticMeanStiffnessTensor(CList,volumes)
      else
        C[coordSub] = compositeType(CList,volumes,normal)
      end

    else
      C[coordSub] = CList[1]
    end
  end

  C
end

function subsamplingComposites(
                               LSuper :: Lattice,
                               LSub :: Lattice,
                               LDecomposition :: Lattice,
                               problem :: LinearElasticityProblem,
                               compositeType :: Type{CompositeType}
                              ) where {CompositeType}
  LinearElasticityProblem(
                          LSub,
                          subsamplingComposites(
                                                LSuper,
                                                LSub,
                                                LDecomposition,
                                                problem.stiffness,
                                                compositeType
                                               ),
                          problem.macroscopicStrain
                         )
end

function subsamplingComposites(
                               LSuper :: Lattice,
                               LSub :: Lattice,
                               problem :: LinearElasticityProblem,
                               compositeType :: Type{CompositeType}
                              ) where {CompositeType}

  LDecomposition = Lattice(round.(Int,LSuper.M*inv(LSub.M)))

  LinearElasticityProblem(
                          LSub,
                          subsamplingComposites(
                                                LSuper,
                                                LSub,
                                                LDecomposition,
                                                problem.stiffness,
                                                compositeType
                                               ),
                          problem.macroscopicStrain
                         )
end
