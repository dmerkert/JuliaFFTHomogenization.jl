export ElasticityHashinEllipsoid

# [geometryFunction] = Elasticity_Hashin_circle(M,...) generates the Hashin structure for the equations of linear elasticity
#
# OPTIONAL ARGUMENTS
#    c1: coordinate system of the ellipse in x direction (1)
#    c2: coordinate system of the ellipse in y direction (1)
#    c3: coordinate system of the ellipse in z direction (1)
#    pInner: elliptic radius of the core wrt. c1,c2,c3 (0.25)
#    pOuter: elliptic radius of the coating wrt. c1,c2,c3 (0.4)
#    nuInner: Poisson's ratio in the core (0.3)
#    nuOuter: Poisson's ratio in the coating (0.3)
#    EInner: Young's modulus in the core (1)
#    EOuter: Young's modulus in the coating (10)
#    normal: normal for rotation ([1;0;0])
#
# OUTPUT
#    geometryFunction: returns a function with
#    C =  geometryFunction(points) where C has the form [3, 3, quadrature point, patternArray] and point is a [component,|det(M)|] matrix
function ElasticityHashinEllipsoid{R <: Real}(lattice :: Lattice;
                                                c1::R=0.5,
                                                c2::R=0.6,
                                                c3::R=Inf,
                                                pInner::R=-0.25+0.2^2,
                                                pOuter::R=-0.25+0.3^2,
                                                nuInner::PoissonsRatio=PoissonsRatio(0.3),
                                                nuOuter::PoissonsRatio=PoissonsRatio(0.3),
                                                EInner::YoungsModulus=YoungsModulus(10.0),
                                                EOuter::YoungsModulus=YoungsModulus(1.0),
                                                normal::Array{R,1}=[1.0,0.0,0.0]
                                               )
  @argcheck c1 > 0
  @argcheck c2 > 0
  @argcheck c3 > 0
  @argcheck length(normal) == 3

  @argcheck pInner < pOuter
  @argcheck -c1^2 >= -c2^2
  @argcheck -c2^2 >= -c3^2
  @argcheck pOuter > -c1^2
  @argcheck pInner > -c1^2

  muInner     = convert(ShearModulus,EInner,nuInner)
  kappaInner  = convert(BulkModulus,EInner,nuInner)
  lambdaInner = convert(LamesFirstParameter,kappaInner,muInner)
  muOuter     = convert(ShearModulus,EOuter,nuOuter)
  kappaOuter  = convert(BulkModulus,EOuter,nuOuter)
  lambdaOuter = convert(LamesFirstParameter,kappaOuter,muOuter)

  #rotation matrix
  rotation = RotationMatrixNormal([1.0;0.0;0.0],normal)

  #semi-axes
  lc1 = sqrt(c1^2+pInner)
  lc2 = sqrt(c2^2+pInner)
  lc3 = sqrt(c3^2+pInner)
  le1 = sqrt(c1^2+pOuter)
  le2 = sqrt(c2^2+pOuter)
  le3 = sqrt(c3^2+pOuter)

  @assert le1 < 0.5
  @assert le2 < 0.5
  @assert le3 < 0.5 || le3 == Inf

  #volume fraction
  f1 = 0
  if isfinite(lc3)
    f1 = (lc1*lc2*lc3)/(le1*le2*le3)
  else
    f1 = (lc1*lc2)/(le1*le2)
  end
  f2 = 1.0-f1

  @assert 0 < f1 < 1
  @assert 0 < f2 < 1

  #depolarization factors
  (dc1,dc2,dc3) = DepolarizationFactors(lc1,lc2,lc3);
  (de1,de2,de3) = DepolarizationFactors(le1,le2,le3);
  M = (diagm([dc1,dc2,dc3])-f1*diagm([de1,de2,de3]))/f2;

  @assert isapprox(trace(M),1.0,rtol=0.0,atol=1e-12)

  epsilon0 = Strain()
  epsilonMatrix =
  (3.0kappaOuter+4.0muOuter)/
  (9.0(kappaInner-kappaOuter))*eye(3) +
  f2*M

  toVoigt!(epsilon0,epsilonMatrix)
  transform!(epsilon0,rotation)

  tau0 = Stress()
  tauMatrix =
  (
   kappaOuter*(kappaInner+4.0/3.0*muOuter)/
   (kappaInner-kappaOuter) +
   4.0muOuter*f1/3.0
  ) * eye(3) +
  2.0muOuter*f2*(M-eye(3)/3.0)

  toVoigt!(tau0,tauMatrix)
  transform!(tau0,rotation)

  EllipsoidMatrixInner = diagm([lc1^-2,lc2^-2,lc3^-2]);
  EllipsoidMatrixOuter = diagm([le1^-2,le2^-2,le3^-2]);
  EllipsoidMatrixInnerRotated = rotation*EllipsoidMatrixInner*rotation';
  EllipsoidMatrixOuterRotated = rotation*EllipsoidMatrixOuter*rotation';

  sigma2 = 1.0;
  sigma1 = 1.0+9.0(kappaInner-kappaOuter)/(3.0kappaOuter+4.0muOuter);

  C = CoefficientTensorField{StiffnessTensor}(lattice.size)
  epsilon = StrainField{Float64}(lattice.size)

  for index in getSamplingIterator(lattice)
    point = getSamplingPoint(lattice,index)/(2.0pi)
    if dot(point,EllipsoidMatrixInnerRotated*point) <= 1.0
      C[index] = IsotropicStiffnessTensor(lambdaInner,muInner)
    elseif dot(point,EllipsoidMatrixOuterRotated*point) <= 1.0
      C[index] = IsotropicStiffnessTensor(lambdaOuter,muOuter)
    else
      v = tau0.val./epsilon0.val
      for i in 1:6
        if epsilon0[i] == 0.0
          v[i] = 1.0
        end
      end
      C[index] = DiagonalStiffnessTensor(v)
    end


    point = rotation'*(getSamplingPoint(lattice,index)/(2.0pi))

    if dot(point,EllipsoidMatrixInner*point) <= 1.0
      s = Strain([1.0/(sigma1-1.0);1.0/(sigma1-1.0);1.0/(sigma1-1.0);0.0;0.0;0.0])
      epsilon[index] = s
    elseif dot(point,EllipsoidMatrixOuter*point) > 1.0
      epsilon[index] = epsilon0
    else
      p = CartesianToEllipsoidal(point,c1,c2,c3)

      c = [c1;c2;c3]
      dc = [dc1;dc2;dc3]

      dpdx = zeros(3)
      if c3 == Inf
        dpdx =
        [2.0point[1:2]./(c[1:2].^2+p).*1.0/(sum(point[1:2].^2./(c[1:2].^2+p).^2));0.0]
      else
        dpdx = 2.0point./(c.^2+p)*1.0/(sum(point.^2./(c.^2+p).^2))
      end

      g = 0
      if c3 == Inf
        g = t -> prod(c[1:2].^2+t)
      else
        g = t -> prod(c.^2+t)
      end

      integrand = 0
      if c3 == Inf
        integrand = [1.0./( (c[1:2].^2+p).*sqrt(prod(c[1:2].^2+p)) );0.0]
      else
        integrand = 1.0./( (c.^2+p).*sqrt(prod(c.^2+p)) )
      end

      s = Strain()
      s.val[4] = 2.0*point[2]*sqrt(g(pInner))./(2.0sigma2)*dpdx[3]*integrand[2]
      s.val[5] = 2.0*point[1]*sqrt(g(pInner))./(2.0sigma2)*dpdx[3]*integrand[1]
      s.val[6] = 2.0*point[1]*sqrt(g(pInner))./(2.0sigma2)*dpdx[2]*integrand[1]

      (dilp1,dilp2,dilp3) = DepolarizationFactors(sqrt(c1^2+p),sqrt(c2^2+p),sqrt(c3^2+p));

      s.val[1] = 1.0./(sigma1-sigma2) +
      sqrt(g(pInner))/(2.0sigma2) *
      (
       2.0dc[1]/sqrt(g(pInner)) -
       2.0dilp1/sqrt(g(p)) +
      point[1]*integrand[1]*dpdx[1]
     )
      s.val[2] = 1.0./(sigma1-sigma2) +
      sqrt(g(pInner))/(2.0sigma2) *
      (
       2.0dc[2]/sqrt(g(pInner)) -
       2.0dilp2/sqrt(g(p)) +
      point[2]*integrand[2]*dpdx[2]
     )
      s.val[3] = 1.0./(sigma1-sigma2) +
      sqrt(g(pInner))/(2.0sigma2) *
      (
       2.0dc[3]/sqrt(g(pInner)) -
       2.0dilp3/sqrt(g(p)) +
      point[3]*integrand[3]*dpdx[3]
     )

      transform!(s,rotation)
      epsilon[index] = s
    end



  end

  problem = LinearElasticityProblem(lattice,C,epsilon0)
  problem.strain = epsilon
  problem.averageStress = tau0
  problem
end

