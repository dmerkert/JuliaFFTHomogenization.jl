export ElasticityBoxSpline

function ElasticityBoxSpline(
                             lattice :: Lattice,
                             scaling :: Lattice,
                             Xi :: NTuple{N,Array{Float64,1}}
                            ) where {N}

  C = CoefficientTensorField{IsotropicStiffnessTensor}(lattice.size)
  λ = LamesFirstParameter(0.3)

  YoungsModuli = Array{Float64}(lattice.size)

  for coord in getFrequencyIterator(lattice)
    freq = getFrequencyPoint(lattice,coord)
    YoungsModuli[coord] = ckPeriodicBoxSpline(scaling.M'\freq,Xi)*lattice.m
  end

  YoungsModuli = real(patternifft(YoungsModuli,lattice)).+1.0
  @assert all(YoungsModuli .> 0.0)

  for coord in getSamplingIterator(lattice)
    C[coord] = IsotropicStiffnessTensor(λ,YoungsModulus(YoungsModuli[coord]))
  end

  strain0 = Strain([1.0;0.0;0.0;0.0;0.0;0.0])

  problem = LinearElasticityProblem(lattice, C,strain0)
end


