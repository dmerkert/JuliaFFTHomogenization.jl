export Elasticity1DLaminate

function Elasticity1DLaminate(lattice :: Lattice)
  @argcheck lattice.d == 3
  @argcheck lattice.dimension == 1
  @argcheck isdiag(lattice.M)
  @argcheck lattice.m == 10




  C = CoefficientTensorField{StiffnessTensor}(lattice.size)
  epsilon = StrainField{Float64}(lattice.size)

  C[1] = IsotropicStiffnessTensor(
                                  YoungsModulus(29.3333333333333),
                                  PoissonsRatio(0.3333333333)
                                 )
  for i in 2:10
    C[i] = IsotropicStiffnessTensor(
                                    YoungsModulus(2.666666666666),
                                    PoissonsRatio(0.3333333333)
                                   )
  end
  problem = LinearElasticityHomogenizationProblemWithFields(lattice, C)


  if lattice.M[1,1] == 10

    CEffective = AnisotropicStiffnessTensor(
                                            [4.4 2.2 2.2 0 0 0;
                                             2.2 7.1 3.1 0 0 0;
                                             2.2 3.1 7.1 0 0 0;
                                             0   0   0   2 0 0;
                                             0   0   0   0 1.1 0;
                                             0   0   0   0 0 1.1
                                            ]
                                           )

    problem.effectiveStiffness = CEffective

    problem.strain[1] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[1])[1] = Strain([0.1,0.0,0.0,0.0,0.0,0.0])
    for j in 2:10
      get(problem.strain[1])[j] = Strain([1.1,0.0,0.0,0.0,0.0,0.0])
    end

    problem.strain[2] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[2])[1] = Strain([-0.45,1.0,0.0,0.0,0.0,0.0])
    for j in 2:10
      get(problem.strain[2])[j] = Strain([0.05,1.0,0.0,0.0,0.0,0.0])
    end

    problem.strain[3] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[3])[1] = Strain([-0.45,0.0,1.0,0.0,0.0,0.0])
    for j in 2:10
      get(problem.strain[3])[j] = Strain([0.05,0.0,1.0,0.0,0.0,0.0])
    end

    problem.strain[4] = Nullable(StrainField{Float64}((10,)))
    for j in 1:10
      get(problem.strain[4])[j] = Strain([0.0,0.0,0.0,1.0,0.0,0.0])
    end

    problem.strain[5] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[5])[1] = Strain([0.0,0.0,0.0,0.0,0.1,0.0])
    for j in 2:10
      get(problem.strain[5])[j] = Strain([0.0,0.0,0.0,0.0,1.1,0.0])
    end

    problem.strain[6] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[6])[1] = Strain([0.0,0.0,0.0,0.0,0.0,0.1])
    for j in 2:10
      get(problem.strain[6])[j] = Strain([0.0,0.0,0.0,0.0,0.0,1.1])
    end


  end

  if lattice.M[2,2] == 10

    CEffective = AnisotropicStiffnessTensor(
                                            [7.1 2.2 3.1 0 0 0;
                                             2.2 4.4 2.2 0 0 0;
                                             3.1 2.2 7.1 0 0 0;
                                             0   0   0   1.1 0 0;
                                             0   0   0   0 2 0;
                                             0   0   0   0 0 1.1
                                            ]
                                           )

    problem.effectiveStiffness = CEffective

    problem.strain[1] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[1])[1] = Strain([1.0,-0.45,0.0,0.0,0.0,0.0])
    for j in 2:10
      get(problem.strain[1])[j] = Strain([1.0,0.05,0.0,0.0,0.0,0.0])
    end

    problem.strain[2] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[2])[1] = Strain([0.0,0.1,0.0,0.0,0.0,0.0])
    for j in 2:10
      get(problem.strain[2])[j] = Strain([0.0,1.1,0.0,0.0,0.0,0.0])
    end


    problem.strain[3] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[3])[1] = Strain([0.0,-0.45,1.0,0.0,0.0,0.0])
    for j in 2:10
      get(problem.strain[3])[j] = Strain([0.0,0.05,1.0,0.0,0.0,0.0])
    end

    problem.strain[4] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[4])[1] = Strain([0.0,0.0,0.0,0.1,0.0,0.0])
    for j in 2:10
      get(problem.strain[4])[j] = Strain([0.0,0.0,0.0,1.1,0.0,0.0])
    end

    problem.strain[5] = Nullable(StrainField{Float64}((10,)))
    for j in 1:10
      get(problem.strain[5])[j] = Strain([0.0,0.0,0.0,0.0,1.0,0.0])
    end

    problem.strain[6] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[6])[1] = Strain([0.0,0.0,0.0,0.0,0.0,0.1])
    for j in 2:10
      get(problem.strain[6])[j] = Strain([0.0,0.0,0.0,0.0,0.0,1.1])
    end
  end

  if lattice.M[3,3] == 10

    CEffective = AnisotropicStiffnessTensor(
                                            [7.1 3.1 2.2 0 0 0;
                                             3.1 7.1 2.2 0 0 0;
                                             2.2 2.2 4.4 0 0 0;
                                             0   0   0   1.1 0 0;
                                             0   0   0   0 1.1 0;
                                             0   0   0   0 0 2
                                            ]
                                           )

    problem.effectiveStiffness = CEffective

    problem.strain[1] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[1])[1] = Strain([1.0,0.0,-0.45,0.0,0.0,0.0])
    for j in 2:10
      get(problem.strain[1])[j] = Strain([1.0,0.0,0.05,0.0,0.0,0.0])
    end

    problem.strain[2] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[2])[1] = Strain([0.0,1.0,-0.45,0.0,0.0,0.0])
    for j in 2:10
      get(problem.strain[2])[j] = Strain([0.0,1.0,0.05,0.0,0.0,0.0])
    end


    problem.strain[3] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[3])[1] = Strain([0.0,0.0,0.1,0.0,0.0,0.0])
    for j in 2:10
      get(problem.strain[3])[j] = Strain([0.0,0.0,1.1,0.0,0.0,0.0])
    end

    problem.strain[4] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[4])[1] = Strain([0.0,0.0,0.0,0.1,0.0,0.0])
    for j in 2:10
      get(problem.strain[4])[j] = Strain([0.0,0.0,0.0,1.1,0.0,0.0])
    end

    problem.strain[5] = Nullable(StrainField{Float64}((10,)))
    get(problem.strain[5])[1] = Strain([0.0,0.0,0.0,0.0,0.1,0.0])
    for j in 2:10
      get(problem.strain[5])[j] = Strain([0.0,0.0,0.0,0.0,1.1,0.0])
    end

    problem.strain[6] = Nullable(StrainField{Float64}((10,)))
    for j in 1:10
      get(problem.strain[6])[j] = Strain([0.0,0.0,0.0,0.0,0.0,1.0])
    end
  end


problem
end
