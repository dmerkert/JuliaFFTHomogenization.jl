using JuliaFFTHomogenization
using Base.Test

@testset "Laminate" begin

C1 = IsotropicStiffnessTensor(
                              YoungsModulus(29.3333333333333),
                              PoissonsRatio(0.3333333333)
                             )

C2 = IsotropicStiffnessTensor(
                              YoungsModulus(2.666666666666),
                              PoissonsRatio(0.3333333333)
                             )

amtC1 = 0.1
amtC2 = 0.9

CEff1 = laminateStiffness(C1,C2,amtC1,amtC2,[1.0;0.0;0.0])
CEff2 = laminateStiffness(C1,C2,amtC1,amtC2,[0.0;1.0;0.0])
CEff3 = laminateStiffness(C1,C2,amtC1,amtC2,[0.0;0.0;1.0])

CEffNum1 = AnisotropicStiffnessTensor(
                                      [4.4 2.2 2.2 0 0 0;
                                       2.2 7.1 3.1 0 0 0;
                                       2.2 3.1 7.1 0 0 0;
                                       0   0   0   2 0 0;
                                       0   0   0   0 1.1 0;
                                       0   0   0   0 0 1.1
                                      ]
                                     )
CEffNum2 = AnisotropicStiffnessTensor(
                                      [7.1 2.2 3.1 0 0 0;
                                       2.2 4.4 2.2 0 0 0;
                                       3.1 2.2 7.1 0 0 0;
                                       0   0   0   1.1 0 0;
                                       0   0   0   0 2 0;
                                       0   0   0   0 0 1.1
                                      ]
                                     )
CEffNum3 = AnisotropicStiffnessTensor(
                                      [7.1 3.1 2.2 0 0 0;
                                       3.1 7.1 2.2 0 0 0;
                                       2.2 2.2 4.4 0 0 0;
                                       0   0   0   1.1 0 0;
                                       0   0   0   0 1.1 0;
                                       0   0   0   0 0 2
                                      ]
                                     )

@test CEff1.C ≈ CEffNum1.C
@test CEff2.C ≈ CEffNum2.C
@test CEff3.C ≈ CEffNum3.C
end
