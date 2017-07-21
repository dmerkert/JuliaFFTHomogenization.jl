using JuliaFFTHomogenization
using Base.Test

@testset "CauchyConvergenceCriterion" begin
  cauchy = CauchyConvergenceCriterion()

  s = SolutionTensorField{Float64,Strain}(rand((6,3,4,5)))
  t = SolutionTensorField{Float64,Strain}(rand((6,3,4,5)))
  z = SolutionTensorField{Float64,Strain}(zeros((6,3,4,5)))
  init!(cauchy, s)
  @test get(cauchy.oldGradient).val ≈ s.val

  set!(cauchy,t)
  @test get(cauchy.oldGradient).val ≈ t.val

  computeError(cauchy,s)
  @test computeError(cauchy,t) ≈ 0.0
  @test computeError(cauchy,z) ≈ Inf
  set!(cauchy,z)
  @test isnan(computeError(cauchy,z))
end

@testset "NormConvergenceCriterion" begin
  n = NormConvergenceCriterion()

  s = SolutionTensorField{Float64,Strain}(rand((6,3,4,5)))
  t = SolutionTensorField{Float64,Strain}(rand((6,3,4,5)))
  z = SolutionTensorField{Float64,Strain}(zeros((6,3,4,5)))
  init!(n, s)
  @test n.oldNorm ≈ norm(s)

  set!(n,t)
  @test n.oldNorm ≈ norm(t)

  computeError(n,s)
  @test computeError(n,t) ≈ 0.0
  @test computeError(n,z) ≈ Inf
  set!(n,z)
  @test isnan(computeError(n,z))
end
