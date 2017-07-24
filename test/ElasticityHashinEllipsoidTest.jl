using JuliaFFTHomogenization
using Base.Test
using MPAWL

@testset "Basic scheme on MacroscopicGradientProblem" begin
  L = Lattice([64 0 0;0 64 0;0 0 1])

  n = [0.5,1,0]
  problem = ElasticityHashinEllipsoid(L;
                                      c1=0.05,
                                      c2=0.35,
                                      c3=Inf,
                                      pInner=0.0,
                                      pOuter=0.09,
                                      normal=n
                                     )
  problemNumeric = copy(problem)

  approximationMethod = ApproximationMethod(FFTTransformation(),Gamma0())
  solver =
  BasicScheme(;printSkip=1,verbose=false,maxIter=200,convergenceCriterion=CauchyConvergenceCriterion())
  solve!(problemNumeric,approximationMethod,solver)
  @test isapprox(
                 average(get(problemNumeric.strain)).val,
                 problemNumeric.macroscopicStrain.val,
                 rtol=1e-12
                )
  @test isapprox(
                 get(problem.strain).val,
                 get(problemNumeric.strain).val,
                 rtol=1e-1
                )

  @test isapprox(
                 get(problem.averageStress).val,
                 get(problemNumeric.averageStress).val,
                 rtol=1e-2
                )

  solver = BasicScheme(;printSkip=1,verbose=false,maxIter=200,convergenceCriterion=NormConvergenceCriterion())
  solve!(problemNumeric,approximationMethod,solver)
  @test isapprox(
                 average(get(problemNumeric.strain)).val,
                 problemNumeric.macroscopicStrain.val,
                 rtol=1e-12
                )
  @test isapprox(
                 get(problem.strain).val,
                 get(problemNumeric.strain).val,
                 rtol=1e-1
                )

  @test isapprox(
                 get(problem.averageStress).val,
                 get(problemNumeric.averageStress).val,
                 rtol=1e-2
                )
end
