using JuliaBasicScheme
using Base.Test
using MAT

@testset "Elasticity Hashin Ellipsoid generation diagonal" begin
  L = Lattice([64 0 0;0 64 0;0 0 1])

  (C,epsilon0) = ElasticityHashinEllipsoid(L)
  (epsilon,tau0) = ElasticityHashinEllipsoidAnalytic(L)

  file = matopen("ElasticityHashinEllipsoidDiag64.mat")
  epsilon0Matlab = read(file, "epsilon0")
  CMatlab = read(file, "C")
  tau0Matlab = read(file, "tau0")
  epsilonMatlab = read(file, "epsilonCorrect")
  close(file)

  @test epsilon0.val ≈ epsilon0Matlab
  @test tau0.val ≈ tau0Matlab 
  for i in CartesianRange(size(C))
    CTmp = C[i]
    CTmpAnisotropic = convert(AnisotropicStiffnessTensor,CTmp)
    CTmpMatlab = CMatlab[:,:,1,i]
    @test CTmpMatlab ≈ CTmpAnisotropic.C
  end

  for i in CartesianRange(size(epsilon))
    epsilonTmp = epsilon[i]
    epsilonMatlabTmp = epsilonMatlab[1,1,:,i]
    @test epsilonTmp.val ≈ epsilonMatlabTmp
  end
end

@testset "Elasticity Hashin Ellipsoid generation Rank 1" begin
  L = Lattice([64 1 0;0 64 0;0 0 1])

  (C,epsilon0) = ElasticityHashinEllipsoid(L)
  (epsilon,tau0) = ElasticityHashinEllipsoidAnalytic(L)

  file = matopen("ElasticityHashinEllipsoidRank1.mat")
  epsilon0Matlab = read(file, "epsilon0")
  CMatlab = read(file, "C")
  tau0Matlab = read(file, "tau0")
  epsilonMatlab = read(file, "epsilonCorrect")
  close(file)

  @test epsilon0.val ≈ epsilon0Matlab
  @test tau0.val ≈ tau0Matlab 
  for i in CartesianRange(size(C))
    CTmp = C[i]
    CTmpAnisotropic = convert(AnisotropicStiffnessTensor,CTmp)
    CTmpMatlab = CMatlab[:,:,1,i]
    @test CTmpMatlab ≈ CTmpAnisotropic.C
  end

  for i in CartesianRange(size(epsilon))
    epsilonTmp = epsilon[i]
    epsilonMatlabTmp = epsilonMatlab[1,1,:,i]
    @test epsilonTmp.val ≈ epsilonMatlabTmp
  end
end

@testset "Basic scheme on MacroscopicGradientProblem" begin
  L = Lattice([64 0 0;0 64 0;0 0 1])
  (C,epsilon0) = ElasticityHashinEllipsoid(L)
  (epsilon,tau0) = ElasticityHashinEllipsoidAnalytic(L)

  problem = LinearElasticityProblem(L,C,epsilon0)
  approximationMethod = ApproximationMethod(FFT(),Gamma0())
  solver = BasicScheme(;printSkip=1,verbose=true,maxIter=1000)
  solve!(problem,approximationMethod,solver)
  @show norm(get(problem.strain).val[:]-epsilon.val[:])
end
