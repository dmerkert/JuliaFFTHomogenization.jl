abstract StiffnessParameter <: Real

type BulkModulus <: StiffnessParameter
  val :: Float64
  function BulkModulus(val)
    @argcheck val > 0.0
    new(val)
  end
end
BulkModulus() = BulkModulus(1.0)

type YoungsModulus <: StiffnessParameter
  val :: Float64
  function YoungsModulus(val)
    new(convert(Float64,val))
  end
end
YoungsModulus() = YoungsModulus(1.0)

type LamesFirstParameter <: StiffnessParameter
  val :: Float64
  function LamesFirstParameter(val)
    convert(Float64,val) >= 0.0 || warning("Lame's first parameter is usually positive")
    new(convert(Float64,val))
  end
end
LamesFirstParameter() = LamesFirstParameter(1.0)

type ShearModulus <: StiffnessParameter
  val :: Float64
  function ShearModulus(val)
    @argcheck val >= convert(Float64,val)
    new(convert(Float64,val))
  end
end
ShearModulus() = ShearModulus(1.0)

type PoissonsRatio <: StiffnessParameter
  val :: Float64
  function PoissonsRatio(val)
    @argcheck -1.0 <= convert(Float64,val) <= 0.5
    new(convert(Float64,val))
  end
end
PoissonsRatio() = PoissonsRatio(0.25)


promote_rule{S <: StiffnessParameter,F <: Number}(::Type{F},::Type{S}) =
promote_type(F,Float64)
promote_rule{S <: StiffnessParameter,T <: StiffnessParameter}(::Type{S},::Type{T}) =
Float64
convert{T <: Real}(::Type{T},P :: StiffnessParameter) =  convert(T,P.val)
convert{S <: StiffnessParameter, T<:AbstractFloat}(::Type{S},a::T) =
S(convert(Float64,a))
+{S <: StiffnessParameter}(a::S,b::S) = convert(Float64,a)+convert(Float64,b)
-{S <: StiffnessParameter}(a::S,b::S) = convert(Float64,a)-convert(Float64,b)
*{S <: StiffnessParameter}(a::S,b::S) = convert(Float64,a)*convert(Float64,b)
/{S <: StiffnessParameter}(a::S,b::S) = convert(Float64,a)/convert(Float64,b)
function convert{S <: StiffnessParameter,
                 T <: StiffnessParameter,
                 U <: StiffnessParameter}(::Type{S}, a :: T, b :: U)
  if S == T
    return a
  elseif S == U
    return b
  else
    convert(S,b,a)
  end
end

convert(::Type{LamesFirstParameter}, K :: BulkModulus, E :: YoungsModulus) =
(3K*(3K-E))/(9K-E) |> LamesFirstParameter

convert(::Type{ShearModulus}, K :: BulkModulus, E :: YoungsModulus) =
3K*E/(9K-E) |> ShearModulus

convert(::Type{PoissonsRatio}, K :: BulkModulus, E :: YoungsModulus) =
(3K-E)/(6K) |> PoissonsRatio

convert(::Type{YoungsModulus}, K :: BulkModulus, l :: LamesFirstParameter) =
(9K*(K-l))/(3K-l) |> YoungsModulus

convert(::Type{ShearModulus}, K :: BulkModulus, l :: LamesFirstParameter) =
3(K-l)/2 |> ShearModulus

convert(::Type{PoissonsRatio}, K :: BulkModulus, l :: LamesFirstParameter) =
l/(3K-l) |> PoissonsRatio

convert(::Type{YoungsModulus}, K :: BulkModulus, G :: ShearModulus) =
9K*G/(3K+G) |> YoungsModulus

convert(::Type{LamesFirstParameter}, K :: BulkModulus, G :: ShearModulus) =
K-2G/3 |> LamesFirstParameter

convert(::Type{PoissonsRatio}, K :: BulkModulus, G :: ShearModulus) =
(3K-2G)/(2(3K+G)) |> PoissonsRatio

convert(::Type{YoungsModulus}, K :: BulkModulus, nu :: PoissonsRatio) =
3K*(1-2nu) |> YoungsModulus

convert(::Type{LamesFirstParameter}, K :: BulkModulus, nu :: PoissonsRatio) =
3K*nu/(1+nu) |> LamesFirstParameter

convert(::Type{ShearModulus}, K :: BulkModulus, nu :: PoissonsRatio) =
(3K*(1-2nu))/(2(1+nu)) |> ShearModulus

convert(::Type{BulkModulus}, E :: YoungsModulus, l :: LamesFirstParameter) =
(E+3l+sqrt(E^2+9l^2+2E*l))/6 |> BulkModulus

convert(::Type{ShearModulus}, E :: YoungsModulus, l :: LamesFirstParameter) =
(E-3l+sqrt(E^2+9l^2+2E*l))/4 |> ShearModulus

convert(::Type{PoissonsRatio}, E :: YoungsModulus, l :: LamesFirstParameter) =
2l/(E+l+sqrt(E^2+9l^2+2E*l)) |> PoissonsRatio

convert(::Type{BulkModulus}, E :: YoungsModulus, G :: ShearModulus) =
E*G/(3(3G-E)) |> BulkModulus

convert(::Type{LamesFirstParameter}, E :: YoungsModulus, G :: ShearModulus) =
G*(E-2G)/(3G-E) |> LamesFirstParameter

convert(::Type{PoissonsRatio}, E :: YoungsModulus, G :: ShearModulus) =
E/(2G)-1 |> PoissonsRatio

convert(::Type{BulkModulus}, E :: YoungsModulus, nu :: PoissonsRatio) =
E/(3(1-2nu)) |> BulkModulus

convert(::Type{LamesFirstParameter}, E :: YoungsModulus, nu :: PoissonsRatio) =
E*nu/((1+nu)*(1-2nu)) |> LamesFirstParameter

convert(::Type{ShearModulus}, E :: YoungsModulus, nu :: PoissonsRatio) =
E/(2(1+nu)) |> ShearModulus

convert(::Type{BulkModulus}, l :: LamesFirstParameter, G :: ShearModulus) =
l+2G/3 |> BulkModulus

convert(::Type{YoungsModulus}, l :: LamesFirstParameter, G :: ShearModulus) =
G*(3l+2G)/(l+G) |> YoungsModulus

convert(::Type{PoissonsRatio}, l :: LamesFirstParameter, G :: ShearModulus) =
l/(2(l+G)) |> PoissonsRatio

convert(::Type{BulkModulus}, l :: LamesFirstParameter, nu :: PoissonsRatio) =
l*(1+nu)/(3nu) |> BulkModulus

convert(::Type{YoungsModulus}, l :: LamesFirstParameter, nu :: PoissonsRatio) =
l*(1+nu)*(1-2nu)/nu |> YoungsModulus

convert(::Type{ShearModulus}, l :: LamesFirstParameter, nu :: PoissonsRatio) =
l*(1-2nu)/(2nu) |> ShearModulus

convert(::Type{BulkModulus}, G :: ShearModulus, nu :: PoissonsRatio) =
2G*(1+nu)/(3(1-2nu)) |> BulkModulus

convert(::Type{YoungsModulus}, G :: ShearModulus, nu :: PoissonsRatio) =
2G*(1+nu) |> YoungsModulus

convert(::Type{LamesFirstParameter}, G :: ShearModulus, nu :: PoissonsRatio) =
2G*nu/(1-2nu) |> LamesFirstParameter

