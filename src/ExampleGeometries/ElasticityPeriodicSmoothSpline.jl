export ElasticityPeriodicSmoothSpline

function ElasticityPeriodicSmoothSpline(
                             lattice :: Lattice,
                             order :: Int
                            )

  C = CoefficientTensorField{IsotropicStiffnessTensor}(lattice.size)
  ν = PoissonsRatio(0.3)

  (p0,p1) = getPolys(order)

  for coord in getSamplingIterator(lattice)
    point = getSamplingPoint(lattice,coord)./(2.0*pi)
    E = YoungsModulus(5.0*evalPolys(point,p0,p1))
    C[coord] = IsotropicStiffnessTensor(ν,E)
  end

  strain0 = Strain([1.0;0.0;0.0;0.0;0.0;0.0])

  problem = LinearElasticityProblem(lattice, C,strain0)
end


