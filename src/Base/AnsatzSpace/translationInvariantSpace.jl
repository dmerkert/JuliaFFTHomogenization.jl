abstract type TranslationInvariantSpace <: AnsatzSpace end
abstract type OrthonormalTranslationInvariantSpace <: TranslationInvariantSpace end

#= function coefficients2Values!( =#
#=                               pointField :: SolutionTensorField{R,G,N}, =#
#=                               space :: Space, =#
#=                               L :: Lattice =#
#=                              ) where { =#
#=                                       R <: AbstractFloat, =#
#=                                       G, =#
#=                                       N, =#
#=                                       Space <: TranslationInvariantSpace =#
#=                                      } =#

#=   dims = [2:(L.rank+1)...] =#
#=   bSq = ansatzSpaceBracketSums(space,Lattice) =#
#=   changeBasis!( =#
#=                pointField.val, =#
#=                pointField.val, =#
#=                L, =#
#=                bSq, =#
#=                dims, =#
#=                inputDomain="Space", =#
#=                outputDomain="Space" =#
#=               ) =#
#=   pointField =#
#= end =#

#= function values2Coefficients!( =#
#=                               pointField :: SolutionTensorField{R,G,N}, =#
#=                               space :: Space, =#
#=                               L :: Lattice =#
#=                              ) where { =#
#=                                       R <: AbstractFloat, =#
#=                                       G, =#
#=                                       N, =#
#=                                       Space <: TranslationInvariantSpace =#
#=                                      } =#
#=   dims = [2:(L.rank+1)...] =#
#=   bSq = ansatzSpaceBracketSums(space,Lattice) =#
#=   changeBasis!( =#
#=                pointField.val, =#
#=                pointField.val, =#
#=                L, =#
#=                1.0./bSq, =#
#=                dims, =#
#=                inputDomain="Space", =#
#=                outputDomain="Space" =#
#=               ) =#
#=   pointField =#
#= end =#

function coefficients2Values!(
                              pointField :: SolutionTensorField{R,G,N},
                              space :: Space,
                              L :: Lattice
                             ) where {
                                      R <: AbstractFloat,
                                      G,
                                      N,
                                      Space <: OrthonormalTranslationInvariantSpace
                                     }
  #nothing to do here
  pointField
end

function values2Coefficients!(
                              pointField :: SolutionTensorField{R,G,N},
                              space :: Space,
                              L :: Lattice
                             ) where {
                                      R <: AbstractFloat,
                                      G,
                                      N,
                                      Space <: OrthonormalTranslationInvariantSpace
                                     }
  #nothing to do here
  pointField
end

function bracketSum(
                    frequency :: Array{I,1},
                    space :: Space,
                    L :: Lattice
                   ) where {
                            I <: Integer,
                            Space <: TranslationInvariantSpace
                           }
  error("bracketSum not implemented!")
end

function ck(
            frequency :: Array{I,1},
            space :: Space,
            L :: Lattice
           ) where {
                    I <: Integer,
                    Space <: TranslationInvariantSpace
                   }

  error("ck not implemented!")
end

function mult!(c :: SolutionTensorField{R,F1,M},
               G :: GreenOperator,
               b :: SolutionTensorField{R,F2,M},
               referenceCoefficient :: C,
               L :: Lattice,
               space :: Space
              ) where {
                       R,
                       F1,
                       F2,
                       M,
                       C<:CoefficientTensor,
                       Space <: TranslationInvariantSpace
                      }
  @argcheck size(b) == size(c)

  tmpGamma = initMult(G)
  tmpc = c[start(getFrequencyIterator(L))] :: F1{R}
  tmpcTmp = c[start(getFrequencyIterator(L))] :: F1{R}
  tmpb = b[start(getFrequencyIterator(L))] :: F2{R}
  tmpFrequency = getFrequencyPoint(L,start(getFrequencyIterator(L))) :: Array{Int,1}
  tmpFrequencyFloat = zeros(Float64,size(tmpFrequency)) :: Array{Float64,1}

  _mult!(c,G,b,referenceCoefficient,L,space,tmpGamma,tmpc,tmpcTmp,tmpb,tmpFrequency,tmpFrequencyFloat,getFrequencyIterator(L))
end

function _mult!(c :: SolutionTensorField{R,F1,M},
               G :: GreenOperator,
               b :: SolutionTensorField{R,F2,M},
               referenceCoefficient :: C,
               L :: Lattice,
               space :: Space,
               tmpGamma :: TMPG,
               tmpc :: F3,
               tmpcTmp :: F3,
               tmpb :: F4,
               tmpFrequency :: Array{I,1},
               tmpFrequencyFloat :: Array{Float64,1},
               range :: CartesianRange{CartesianIndex{N}}
              ) where {
                       R,
                       F1,
                       F2,
                       M,
                       C<:CoefficientTensor,
                       TMPG,
                       N,
                       F3,
                       F4,
                       I,
                       Space <: TranslationInvariantSpace
  }

  for i in range
    getFrequencyPoint!(L,i,tmpFrequency,tmpFrequencyFloat)
    tmpc.val .= 0.0
    for bs in BracketSumIterator(tmpFrequency,-2,2,L)
      tmpb.val .= b.val[1:6,i]
      abs2ck = abs2(ck(bs,space,L))
      if abs2ck > 0.0
        mult!(tmpcTmp,G,tmpb,referenceCoefficient,bs,tmpGamma)
        tmpc.val .+= tmpcTmp.val * abs2ck
      end
    end
    c.val[1:6,i] .= tmpc.val
  end
  c
end



