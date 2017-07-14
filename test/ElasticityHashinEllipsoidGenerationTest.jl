using JuliaFFTHomogenization
using Base.Test
using MAT
using MPAWL

@testset "Elasticity Hashin Ellipsoid generation diagonal" begin
  L = Lattice([64 0 0;0 64 0;0 0 1])

  problem = ElasticityHashinEllipsoid(L)

  file = matopen("ElasticityHashinEllipsoidDiag64.mat")
  epsilon0Matlab = read(file, "epsilon0")
  CMatlab = read(file, "C")
  tau0Matlab = read(file, "tau0")
  epsilonMatlab = read(file, "epsilonCorrect")
  close(file)

  @test problem.macroscopicStrain.val ≈ epsilon0Matlab
  @test get(problem.averageStress).val ≈ tau0Matlab
  for i in CartesianRange(size(problem.stiffness))
    CTmp = problem.stiffness[i]
    CTmpAnisotropic = convert(AnisotropicStiffnessTensor,CTmp)
    CTmpMatlab = CMatlab[:,:,1,i]
    @test CTmpMatlab ≈ CTmpAnisotropic.C
  end

  strainFieldMatlab = StrainField{Float64}(size(get(problem.strain)))
  for i in CartesianRange(size(get(problem.strain)))
    epsilonTmp = get(problem.strain)[i]
    epsilonMatlabTmp = epsilonMatlab[1,1,:,i]
    strainFieldMatlab[i] = Strain(epsilonMatlabTmp)
    @test epsilonTmp.val ≈ epsilonMatlabTmp
  end
end

@testset "Elasticity Hashin Ellipsoid generation Rank 1" begin
  L = Lattice([64 1 0;0 64 0;0 0 1])

  problem = ElasticityHashinEllipsoid(L)

  file = matopen("ElasticityHashinEllipsoidRank1.mat")
  epsilon0Matlab = read(file, "epsilon0")
  CMatlab = read(file, "C")
  tau0Matlab = read(file, "tau0")
  epsilonMatlab = read(file, "epsilonCorrect")
  close(file)

  @test problem.macroscopicStrain.val ≈ epsilon0Matlab
  @test get(problem.averageStress).val ≈ tau0Matlab
  for i in CartesianRange(size(problem.stiffness))
    CTmp = problem.stiffness[i]
    CTmpAnisotropic = convert(AnisotropicStiffnessTensor,CTmp)
    CTmpMatlab = CMatlab[:,:,1,i]
    @test CTmpMatlab ≈ CTmpAnisotropic.C
  end

  for i in CartesianRange(size(get(problem.strain)))
    epsilonTmp = get(problem.strain)[i]
    epsilonMatlabTmp = epsilonMatlab[1,1,:,i]
    @test epsilonTmp.val ≈ epsilonMatlabTmp
  end
end

