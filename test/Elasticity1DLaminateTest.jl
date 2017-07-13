using JuliaFFTHomogenization
using Base.Test
using MPAWL

@testset "Elasticity 1D Laminate Test x" begin
  L = Lattice(diagm([10,1,1]))

  problem = Elasticity1DLaminate(L)
  problemNumeric = copy(problem)

  approximationMethod = ApproximationMethod(FFTTransformation(),Gamma0())
  solver = BasicScheme(;printSkip=1,verbose=false,maxIter=1000,tol=1e-10)
  solve!(problemNumeric,approximationMethod,solver)

  CEff = get(problem.effectiveStiffness).C
  CEffNumeric = get(problemNumeric.effectiveStiffness).C
  @test CEff ≈ CEffNumeric

  for i in 1:6
    s =  get(problem.strain[i]).val
    sNumeric = get(problemNumeric.strain[i]).val
    @test s ≈ sNumeric
  end

end

@testset "Elasticity 1D Laminate Test y" begin
  L = Lattice(diagm([1,10,1]))

  problem = Elasticity1DLaminate(L)
  problemNumeric = copy(problem)

  approximationMethod = ApproximationMethod(FFTTransformation(),Gamma0())
  solver = BasicScheme(;printSkip=1,verbose=false,maxIter=1000,tol=1e-10)
  solve!(problemNumeric,approximationMethod,solver)

  CEff = get(problem.effectiveStiffness).C
  CEffNumeric = get(problemNumeric.effectiveStiffness).C
  @test CEff ≈ CEffNumeric

  for i in 1:6
    s =  get(problem.strain[i]).val
    sNumeric = get(problemNumeric.strain[i]).val
    @test s ≈ sNumeric
  end

end

@testset "Elasticity 1D Laminate Test z" begin
  L = Lattice(diagm([1,1,10]))

  problem = Elasticity1DLaminate(L)
  problemNumeric = copy(problem)

  approximationMethod = ApproximationMethod(FFTTransformation(),Gamma0())
  solver = BasicScheme(;printSkip=1,verbose=false,maxIter=1000,tol=1e-10)
  solve!(problemNumeric,approximationMethod,solver)

  CEff = get(problem.effectiveStiffness).C
  CEffNumeric = get(problemNumeric.effectiveStiffness).C
  @test CEff ≈ CEffNumeric

  for i in 1:6
    s =  get(problem.strain[i]).val
    sNumeric = get(problemNumeric.strain[i]).val
    @test s ≈ sNumeric
  end
end
