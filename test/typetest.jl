using JuliaBasicScheme
using Base.Test

@testset "Types" begin
  @test promote_type(Float32,BulkModulus) == Float64
  promote(2.0,BulkModulus(3.0))
  2.0+3BulkModulus(4.0)
  BulkModulus(4.0)
  convert(BulkModulus,LamesFirstParameter(),YoungsModulus())
  convert(BulkModulus,BulkModulus(),YoungsModulus())
  convert(BulkModulus,PoissonsRatio(),BulkModulus())

  K = BulkModulus(1.2345)
  E = YoungsModulus(3.4323)
  l = LamesFirstParameter(0.243443)
  nu = PoissonsRatio(0.223232)
  mu = ShearModulus(6.34331247)

  @test convert(BulkModulus,
                convert(LamesFirstParameter,K,E),
                convert(PoissonsRatio,K,E)).val ≈ K.val
  @test convert(YoungsModulus,
                convert(LamesFirstParameter,K,E),
                convert(PoissonsRatio,K,E)).val ≈ E.val

  CIso = IsotropicStiffnessTensor(l,mu)
  CTrans = convert(TransversalIsotropicZStiffnessTensor,CIso)

  stress  = Stress()
  strain  = Strain(Float64)

  strain.val = [1.234234,3.3247293,9.12398217,4.12837192,6.2394732,22.2374834]

  @test mult!(stress,CIso,strain).val ≈ mult!(stress,CTrans,strain).val

  @test convert(AnisotropicStiffnessTensor,CIso).C ≈ convert(AnisotropicStiffnessTensor,CTrans).C
  

end

@testset "Strain and Stress Fields" begin
  strain = Strain()
  strainF = StrainField(Float64,(2,2,2))
  stressF = StressField(Float64,(2,2,2))
  I = CartesianIndex((1,1,1))
  strainF[1,1,1] = strain
  strainF[I] = strain
  strain = strainF[1,1,1]
  strain = strainF[I]
  strainF+strainF
  strain.val = [1,2,3,4,5,6]
  init!(strainF,strain)
end

@testset "Coefficient Tensor Fields" begin
  l = LamesFirstParameter(0.243443)
  mu = ShearModulus(6.34331247)

  CIso = IsotropicStiffnessTensor(l,mu)
  CTrans = convert(TransversalIsotropicZStiffnessTensor,CIso)

  CIsoField = CoefficientTensorField(IsotropicStiffnessTensor,(3,3))
  
end
