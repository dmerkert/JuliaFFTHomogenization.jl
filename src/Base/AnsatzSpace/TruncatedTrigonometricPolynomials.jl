export TruncatedTrigonometricPolynomials,
mult!,
kernel,
ck

type TruncatedTrigonometricPolynomials <: AnsatzSpace end

function coefficients2Values!(
                              pointField :: SolutionTensorField{R,G,N},
                              space :: TruncatedTrigonometricPolynomials,
                              L :: Lattice
                             ) where {
                                      R <: AbstractFloat,
                                      G,
                                      N
                                     }
  #nothing to do here
  pointField
end

function values2Coefficients!(
                              pointField :: SolutionTensorField{R,G,N},
                              space :: TruncatedTrigonometricPolynomials,
                              L :: Lattice
                             ) where {
                                      R <: AbstractFloat,
                                      G,
                                      N
                                     }
  #nothing to do here
  pointField
end

function mult!(c :: SolutionTensorField{R,F1,M},
               G :: GreenOperator,
               b :: SolutionTensorField{R,F2,M},
               referenceCoefficient :: C,
               L :: Lattice,
               space :: TruncatedTrigonometricPolynomials
              ) where {
                       R,
                       F1,
                       F2,
                       M,
                       C<:CoefficientTensor
                      }
  @argcheck size(b) == size(c)

  tmpGamma = initMult(G)
  tmpc = c[start(getFrequencyIterator(L))] :: F1{R}
  tmpb = b[start(getFrequencyIterator(L))] :: F2{R}
  tmpFrequency = getFrequencyPoint(L,start(getFrequencyIterator(L))) :: Array{Int,1}
  tmpFrequencyFloat = zeros(Float64,size(tmpFrequency)) :: Array{Float64,1}

  _mult!(c,G,b,referenceCoefficient,L,tmpGamma,tmpc,tmpb,tmpFrequency,tmpFrequencyFloat,getFrequencyIterator(L))
end

function _mult!(c :: SolutionTensorField{R,F1,M},
               G :: GreenOperator,
               b :: SolutionTensorField{R,F2,M},
               referenceCoefficient :: C,
               L :: Lattice,
               tmpGamma :: TMPG,
               tmpc :: F3,
               tmpb :: F4,
               tmpFrequency :: Array{I,1},
               tmpFrequencyFloat :: Array{Float64,1},
               range :: CartesianRange{CartesianIndex{N}}
              ) where {R,F1,F2,M,C<:CoefficientTensor,TMPG,N,F3,F4,I}

  for i in range
    #= getFrequencyPoint!(L,i,tmpFrequency,tmpFrequencyFloat) =#
    tmpFrequency = getFrequencyPoint(L,i)
    tmpb.val .= b.val[1:6,i]
    mult!(tmpc,G,tmpb,referenceCoefficient,tmpFrequency,tmpGamma)
    c.val[1:6,i] .= tmpc.val
  end
  c
end

@inline function kernel(
                x :: Array{R,1},
                space :: TruncatedTrigonometricPolynomials
               ) where {
                        R <: Real
                       }

  if all(x .< 1.0/2.0) && all(x .>= -1.0/2.0)
    return 1.0
  else 
    return 0.0
  end
end

@inline function ck(
            k :: Array{I,1},
            L :: Lattice{LI,MF,MF2},
            space :: TruncatedTrigonometricPolynomials
           ) where {
                    I <: Integer,
                    LI,MF,MF2
                   }
  kernel(L.MTFactorize\k,space)
end
