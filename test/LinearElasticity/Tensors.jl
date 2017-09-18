using JuliaFFTHomogenization
using Base.Test

function f(i,j,k,l)
  2.0i^2 +
  2.3j^2 +
  2.6k^2 +
  2.7l^2 +
  (i*j-3.0k*l)*i^2
end

function f1(i,j,k,l)
  1.0
end

@testset "Tensors" begin
  T = toVoigtSSymFourthOrder(f)
  TI = toInverseVoigtSSymFourthOrder(f)

  @test issymmetric(T)
  @test issymmetric(TI)

  T = toVoigtSSymFourthOrder(f1)
  TI = toInverseVoigtSSymFourthOrder(f1)

  @test all(T .== 1.0)
  @test all(TI[1:3,1:3] .== 1.0)
  @test all(TI[4:6,4:6] .== 4.0)
  @test all(TI[1:3,4:6] .== 2.0)
  @test all(TI[4:6,1:3] .== 2.0)
end
