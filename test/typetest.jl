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
  strain  = Strain()

  strain.val = [1.234234,3.3247293,9.12398217,4.12837192,6.2394732,22.2374834]

  @test mult!(stress,CIso,strain).val ≈ mult!(stress,CTrans,strain).val

  @test convert(AnisotropicStiffnessTensor,CIso).C ≈ convert(AnisotropicStiffnessTensor,CTrans).C
  

end

@testset "Strain and Stress Fields" begin
  strain = Strain()
  @show "a"
  strainF = StrainField{Float64}((2,2,2))
  @show "a"
  stressF = StressField{Float64}((2,2,2))
  @show "a"
  I = CartesianIndex((1,1,1))
  @show "a"
  strainF[1,1,1] = strain
  @show "a"
  strainF[I] = strain
  @show "a"
  strain = strainF[1,1,1]
  @show "a"
  strain = strainF[I]
  @show "a"
  strainF+strainF
  @show "a"
  strain.val = [1,2,3,4,5,6]
  @show "a"
  init!(strainF,strain)
  @show "a"

  b = copy(strainF)
  @show "a"
end

@testset "Coefficient Tensor Fields" begin
  l = LamesFirstParameter(0.243443)
  mu = ShearModulus(6.34331247)

  CIso = IsotropicStiffnessTensor(l,mu)
  CIso2 = IsotropicStiffnessTensor(l+1,mu-1)
  CTrans = convert(TransversalIsotropicZStiffnessTensor,CIso)
  CTrans2 = convert(TransversalIsotropicZStiffnessTensor,CIso2)

  CIsoField = CoefficientTensorField{IsotropicStiffnessTensor}((3,3))
  CTransField = CoefficientTensorField{TransversalIsotropicZStiffnessTensor}((3,3))
  for i in 1:3
    for j in 1:3
      if i == j
        CIsoField[i,j] = CIso
        CTransField[i,j] = CTrans
      else
        CIsoField[i,j] = CIso2
        CTransField[i,j] = CTrans2
      end

    end
  end
  solver = BasicScheme()
  C0 = getReferenceTensor(CIsoField,solver)
  C0Trans = getReferenceTensor(CTransField,solver)

  C0An = convert(AnisotropicStiffnessTensor,C0)
  C0An2 = convert(AnisotropicStiffnessTensor,C0Trans)
end
