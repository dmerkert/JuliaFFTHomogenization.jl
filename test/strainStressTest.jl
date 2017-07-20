using JuliaFFTHomogenization
using Base.Test

@testset "Strain and Stress" begin
  strain = Strain()
  for i in 1:6
    strain.val[i] = 1.0i
    @test strain[i] == strain.val[i]
  end
  @test strain[1,1] == strain.val[1]
  @test 2.0strain[1,2] == strain.val[6]
  @test 2.0strain[1,3] == strain.val[5]
  @test 2.0strain[2,1] == strain.val[6]
  @test strain[2,2] == strain.val[2]
  @test 2.0strain[2,3] == strain.val[4]
  @test 2.0strain[3,1] == strain.val[5]
  @test 2.0strain[3,2] == strain.val[4]
  @test strain[3,3] == strain.val[3]

  for i in 1:6
    strain[i] = 1.0i+10.0
    @test strain.val[i] == 1.0i+10.0
  end

  for i in 1:3
    for j in 1:3
      r = rand()
      strain[i,j] = r
      i == j && @test strain.val[i] == r
      i == 1 && j == 2 && @test strain.val[6] == 2.0r
      i == 1 && j == 3 && @test strain.val[5] == 2.0r
      i == 2 && j == 1 && @test strain.val[6] == 2.0r
      i == 2 && j == 3 && @test strain.val[4] == 2.0r
      i == 3 && j == 1 && @test strain.val[5] == 2.0r
      i == 3 && j == 2 && @test strain.val[4] == 2.0r
    end
  end


  stress = Stress()
  for i in 1:6
    stress.val[i] = 1.0i
    @test stress[i] == stress.val[i]
  end
  @test stress[1,1] == stress.val[1]
  @test stress[1,2] == stress.val[6]
  @test stress[1,3] == stress.val[5]
  @test stress[2,1] == stress.val[6]
  @test stress[2,2] == stress.val[2]
  @test stress[2,3] == stress.val[4]
  @test stress[3,1] == stress.val[5]
  @test stress[3,2] == stress.val[4]
  @test stress[3,3] == stress.val[3]

  for i in 1:6
    stress[i] = 1.0i+10.0
    @test stress.val[i] == 1.0i+10.0
  end

  for i in 1:3
    for j in 1:3
      r = rand()
      stress[i,j] = r
      i == j && @test stress.val[i] == r
      i == 1 && j == 2 && @test stress.val[6] == r
      i == 1 && j == 3 && @test stress.val[5] == r
      i == 2 && j == 1 && @test stress.val[6] == r
      i == 2 && j == 3 && @test stress.val[4] == r
      i == 3 && j == 1 && @test stress.val[5] == r
      i == 3 && j == 2 && @test stress.val[4] == r
    end
  end
end
