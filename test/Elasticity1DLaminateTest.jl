using JuliaFFTHomogenization
using Base.Test
using MPAWL

function test1DLaminate(L,
                        convergenceCriterion,
                        ansatzSpace :: AnsatzSpace
                       )

  problem = Elasticity1DLaminate(L)
  problemNumeric = copy(problem)

  approximationMethod = ApproximationMethod(ansatzSpace,Gamma0())
  solver = BasicScheme(
                       printSkip = 1,
                       verbose = true,
                       maxIter = 1000,
                       tol = 1e-10,
                       convergenceCriterion = convergenceCriterion
                      )

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

function test1DLaminateDirections(
                                  convergenceCriterion,
                                  ansatzSpace
                                 )
  L = Lattice(diagm([10,1,1]))
  test1DLaminate(L,convergenceCriterion,ansatzSpace)
  L = Lattice(diagm([1,10,1]))
  test1DLaminate(L,convergenceCriterion,ansatzSpace)
  L = Lattice(diagm([1,1,10]))
  test1DLaminate(L,convergenceCriterion,ansatzSpace)
end



@testset "Elasticity 1D Laminate Test" begin
  test1DLaminateDirections(
                           CauchyConvergenceCriterion(),
                           TruncatedTrigonometricPolynomials()
                          )

  test1DLaminateDirections(
                           NormConvergenceCriterion(),
                           TruncatedTrigonometricPolynomials()
                          )

  test1DLaminateDirections(
                           CauchyConvergenceCriterion(),
                           deLaValleePoussinMeansSpace([0.25;0.4;0.4])
                          )

end
