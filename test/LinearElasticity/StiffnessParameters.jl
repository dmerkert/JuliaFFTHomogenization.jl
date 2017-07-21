using JuliaFFTHomogenization
using Base.Test

@testset "StiffnessParameters" begin
  BulkModulus()
  YoungsModulus()
  LamesFirstParameter()
  ShearModulus()
  PoissonsRatio()

  BulkModulus(1.0)
  YoungsModulus(1.0)
  LamesFirstParameter(0.5)
  ShearModulus(0.5)
  PoissonsRatio(0.5)

  types = (
           :BulkModulus,
           :YoungsModulus,
           :LamesFirstParameter,
           :ShearModulus,
           :PoissonsRatio
          )

  for t1 in types
    for t2 in types
      for t3 in types
        @eval b = $t2()
        @eval c = $t3()
        if (t2 != t3)
          @eval a = convert($t1,b,c)

          if (t1 != t2) && (t1 != t3)
            @eval @test convert($t2,a,c).val ≈ b.val
          end
        end

        @test b+c ≈ b.val+c.val
        @test b-c ≈ b.val-c.val
        @test b*c ≈ b.val*c.val
        @test b/c ≈ b.val/c.val
      end
    end
  end

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

end
