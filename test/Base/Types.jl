using JuliaFFTHomogenization
using MPAWL
using Base.Test
using BenchmarkTools

@testset "Types" begin
  s = SolutionTensorField{Float64,Strain}(rand((6,3,4,5)))
  t = SolutionTensorField{Float64,Strain}(rand((6,3,4,5)))
  stress = SolutionTensorField{Float64,Stress}(rand((6,3,4,5)))

  copy!(s,t)
  @test s.val == t.val
  s = copy(t)
  @test s.val == t.val
  @test norm(s) == norm(s.val[:])
  norm(s[1,1,1])

  c1 = CoefficientTensorField{IsotropicStiffnessTensor}((3,4,5))
  c2 = CoefficientTensorField(c1.val)

  l = LamesFirstParameter(0.243443)
  mu = ShearModulus(6.34331247)
  CIso = IsotropicStiffnessTensor(l,mu)

  c1[1,1,1] = CIso
  CIso2 = c1[1,1,1]
  @test CIso == CIso2

  u = s+t
  @test u.val == s.val+t.val
  tC = SolutionTensorField{Complex128,Strain}(rand(Complex128,(6,3,4,5)))
  u = s+tC
  @test u.val == s.val+tC.val

  for i in CartesianRange((3,4,5))
    c1[i] = CIso
    c2[i] = CIso
  end

  c1+c2
  c1+CIso

  mult!(stress,c1,s)

  strain = SolutionTensorField{Complex128,Strain}(rand(Complex128,(6,3,4,5)))
  stress = SolutionTensorField{Complex128,Stress}(rand(Complex128,(6,3,4,5)))
  L = Lattice(diagm([3;4;5]))

  mult!(strain,Gamma0(),stress,CIso,L)

  init!(StrainField{Float64}((2,3,4)),Strain())
  init!(StressField{Float64}((2,3,4)),Stress())
  init!(DisplacementField{Float64}((2,3,4)),Displacement())

  @test_throws MethodError init!(StrainField{Float64}((2,3,4)),Stress())
  @test_throws MethodError init!(StressField{Float64}((2,3,4)),Strain())
  @test_throws MethodError init!(DisplacementField{Float64}((2,3,4)),Stress())
end

