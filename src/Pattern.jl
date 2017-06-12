type TensorProductGrid{I <: Integer} <: Pattern
  N :: Array{I,1}

  function TensorProductGrid{I}(N :: Array{I,1})
    @argcheck all(N > 0.0)
    new(N)
  end
end
TensorProductGrid{I}(N :: Array{I,1}) = TensorProductGrid{I}(N)

function PatternIterator(pattern :: TensorProductGrid)
  CartesianRange(tuple(pattern.N...))
end

function PatternPointOfIndex!{R <: Real}(
                                         point :: Array{R,1},
                                         pattern :: TensorProductGrid,
                                         index :: CartesianIndex
                                        )

  point = collect(index.I)./pattern.N
end

function GeneratingSetPointOfIndex!{C <: Complex}(
                                         point :: Array{C,1},
                                         pattern :: TensorProductGrid,
                                         index :: CartesianIndex
                                        )

  point = collect(index.I)./pattern.N
end
