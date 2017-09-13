export ElasticitySingleFiber

function ElasticitySingleFiber{
                               R <: Real
                              }(
                                lattice :: Lattice,
                                epsilon0;
                                radius :: R = 0.1,
                                fiberLength :: R = 0.5,
                                nuInner::PoissonsRatio=PoissonsRatio(0.3),
                                nuOuter::PoissonsRatio=PoissonsRatio(0.3),
                                EInner::YoungsModulus=YoungsModulus(10.0),
                                EOuter::YoungsModulus=YoungsModulus(1.0),
                                direction::Array{R,1}=[1.0;1.0;0.0]
                               )
  @argcheck radius > 0
  @argcheck fiberLength > 0

  @argcheck length(direction) == 3

  #= muInner     = convert(ShearModulus,EInner,nuInner) =#
  #= kappaInner  = convert(BulkModulus,EInner,nuInner) =#
  #= lambdaInner = convert(LamesFirstParameter,kappaInner,muInner) =#
  #= muOuter     = convert(ShearModulus,EOuter,nuOuter) =#
  #= kappaOuter  = convert(BulkModulus,EOuter,nuOuter) =#
  #= lambdaOuter = convert(LamesFirstParameter,kappaOuter,muOuter) =#

  C = CoefficientTensorField{IsotropicStiffnessTensor}(lattice.size)

  for index in getSamplingIterator(lattice)
    point = getSamplingPoint(lattice,index)/(2.0pi)

    #line: f(t) = t*direction
    #projection onto line
    t = dot(point,direction)/dot(direction,direction)
    #radius
    r = norm(point-t*direction)
    #distance
    dist = norm(t*direction)



    if (r <= radius) && (dist <= 0.5fiberLength)
      C[index] = IsotropicStiffnessTensor(EInner,nuInner)
    else
      C[index] = IsotropicStiffnessTensor(EOuter,nuOuter)
    end
  end

  LinearElasticityProblem(lattice,C,epsilon0)
end

