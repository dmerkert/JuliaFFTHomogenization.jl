using JuliaFFTHomogenization
using Base.Test

@testset "Strain and Stress Fields" begin
  StrainField{Float64}((2,3,4))
  StressField{Float64}((2,3,4))
  DisplacementField{Float64}((2,3,4))

  strain = StrainField{Float64}((2,3,4),3.0)
  stress = StressField{Float64}((2,3,4),3.0)
  disp = DisplacementField{Float64}((2,3,4),3.0)

  @test StrainField{Float64}(stress).val == stress.val
  @test StressField{Float64}(strain).val == strain.val
  @test DisplacementField{Float64}(disp).val == disp.val

  @test size(StrainField{Float64}((2,3,4),3.0)) == (2,3,4)
  @test size(StressField{Float64}((2,3,4),3.0)) == (2,3,4)
  @test size(DisplacementField{Float64}((2,3,4),3.0)) == (2,3,4)

  @test sum(abs.(average(strain).val)) ≈ 6*3.0
  @test sum(abs.(average(stress).val)) ≈ 6*3.0

  strain[1,1,1] = Strain()
  strain[1,1,1]
  stress[1,1,1] = Stress()
  stress[1,1,1]
  disp[1,1,1] = Displacement()
  disp[1,1,1]
end
