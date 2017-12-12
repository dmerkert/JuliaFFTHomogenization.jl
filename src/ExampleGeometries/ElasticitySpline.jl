export ElasticityBoxSpline

function ElasticityBoxSpline(
                             lattice :: Lattice,
                             scaling :: Lattice,
                             Xi :: NTuple{N,Array{Float64,1}}
                            ) where {N}

  C = CoefficientTensorField{IsotropicStiffnessTensor}(lattice.size)
  ν = PoissonsRatio(0.3)

  YoungsModuli = Array{Float64}(lattice.size)

  for coord in getFrequencyIterator(lattice)
    freq = getFrequencyPoint(lattice,coord)
    YoungsModuli[coord] = ckPeriodicBoxSpline(scaling.M'\freq,Xi)
  end

  YoungsModuli[1] = 2.0*lattice.m

  YoungsModuli = real(patternifft(YoungsModuli,lattice))
  @assert all(YoungsModuli .> 0.0)

  for coord in getSamplingIterator(lattice)
    C[coord] = IsotropicStiffnessTensor(ν,YoungsModulus(YoungsModuli[coord]))
  end

  strain0 = Strain([1.0;0.0;0.0;0.0;0.0;0.0])

  problem = LinearElasticityProblem(lattice, C,strain0)
end


