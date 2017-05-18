abstract StiffnessTensor <: CoefficientTensor
abstract StiffnessParameter


type BulkModulus{R <: AbstractFloat} <: StiffnessParameter
  val :: R
  function BulkModulus{R}(val :: R)
    val > 0 || error("Bulk modulus must be > 0")
    new(val)
  end
end
BulkModulus{R <: AbstractFloat}(val :: R) = BulkModulus{R}(val)

type YoungsModulus{R <: AbstractFloat} <: StiffnessParameter
  val :: R
end

type LamesFirstParameter{R <: AbstractFloat} <: StiffnessParameter
  val :: R
  function LamesFirstParameter{R}(val :: R)
    val >= 0.0 || warning("Lame's first parameter is usually positive")
    new(val)
  end
end
LamesFirstParameter{R <: AbstractFloat}(val :: R) = LamesFirstParameter{R}(val)

type ShearModulus{R <: AbstractFloat} <: StiffnessParameter
  val :: R
  function ShearModulus{R}(val :: R)
    val >= 0.0 || error("Shear modulus must be >= 0")
    new(val)
  end
end
ShearModulus{R <: AbstractFloat}(val :: R) = ShearModulus{R}(val)

type PoissonsRatio{R <: AbstractFloat} <: StiffnessParameter
  val :: R
  function PoissonsRatio{R}(val :: R)
    -1.0 <= val <= 0.5 || error("Poisson's ratio must be => -1 and <= 0.5")
    new(val)
  end
end
PoissonsRatio{R <: AbstractFloat}(val :: R) = PoissonsRatio{R}(val)

function convert!{R <: StiffnessParameter}(a :: R,b :: R,c :: StiffnessParameter)
  a.val = b.val
end
function convert!{R <: StiffnessParameter}(a :: R,b :: StiffnessParameter,c :: R)
  a.val = c.val
end
function convert!(K :: BulkModulus,E :: YoungsModulus,l :: LamesFirstParameter)
  K.val = (E.val+3.0l.val+sqrt(E.val^2+9.0l.val^2+2.0E.val*l.val))/6.0
end
function convert!(K :: BulkModulus,E :: YoungsModulus,mu :: ShearModulus)
  K.val = E.val*mu.val/(3.0*(3.0mu.val-E.val))
end
function convert!(K :: BulkModulus,E :: YoungsModulus,nu :: PoissonsRatio)
  K.val = E.val/(3.0*(1.0-2.0nu.val))
end




type Strain{R <: Number} <: GradientSolutionTensor
  val :: Array{R,1}

  function Strain{R}(val :: Array{R,1})
    length(val) == 6 || error("Strain vector has wrong size")
    new(val)
  end
end
Strain{R}(val :: Array{R,1}) = Strain{R}(val)

function Strain(R)
  Strain(Array(R,6))
end

function Strain{R <: Number}(val :: Array{R,2})
  size(val) == (3,3) || error("Wrong size")
  issymmetric(val) || error("Only symmetric strains supported")

  Strain([val[[1,5,9]];2*val[[8,7,4]]])
end

type Stress{R <: Number} <: FluxSolutionTensor
  val :: Array{R,1}

  function Stress{R}(val :: Array{R,1})
    length(val) == 6 || error("Stress vector has wrong size")
    new(val)
  end

end
Stress{R}(val :: Array{R,1}) = Stress{R}(val)

function Stress(R)
  Stress(Array(R,6))
end

function Stress{R <: Number}(val :: Array{R,2})
  size(val) == (3,3) || error("Wrong size")
  issymmetric(val) || error("Only symmetric stresses supported")

  Stress([val[[1,5,9]];val[[8,7,4]]])
end

type Displacement{R <: Number} <: PrimarySolutionTensor
  val :: Array{R,1}

  function Displacement{R}(val :: Array{R,1})
    length(val) == 3 || error("Displacement vector has wrong size")
    new(val)
  end
end
Displacement{R}(val :: Array{R,1}) = Displacement{R}(val)

function Displacement(R)
  Displacement(Array(R,3))
end

immutable IsotropicStiffnessTensor{R <: AbstractFloat} <: StiffnessTensor
  lambda :: LamesFirstParameter{R}
  mu     :: ShearModulus{R}
end

immutable TransversalIsotropicZ{R <: AbstractFloat} <: StiffnessTensor
  E_p   :: YoungsModulus{R}
  E_t   :: YoungsModulus{R}
  nu_p  :: PoissonsRatio{R}
  nu_pt :: PoissonsRatio{R}
  nu_tp :: PoissonsRatio{R}
  mu_t  :: ShearModulus{R}
end

type StiffnessTensorField{R <: StiffnessTensor} <: CoefficientTensorField
  val :: Array{R}
end
StiffnessTensorField(dims) = 
StiffnessTensorField{StiffnessTensor}(Array{StiffnessTensor}(dims))
StiffnessTensorField(R,dims) = 
StiffnessTensorField{R}(Array{R}(dims))
get(stiffnessTensorField :: StiffnessTensorField,Index) =
stiffnessTensorField.val[Index]
size(stiffnessTensorField :: StiffnessTensorField) =
size(stiffnessTensorField.val)

type StrainField{R <: Number} <: GradientSolutionTensorField
  val :: Array{R}
end
StrainField(R,dims) = StrainField(Array{R}(tuple(dims...,6)))
function get!{R}(strain :: Strain{R},strainField :: StrainField{R},Index)
  strain.val = strainField.val[Index,:]
end
function set!{R}(strainField :: StrainField{R},strain :: Strain{R},Index)
  strainField.val[Index,:] = strain.val
end
function add!{R}(strain :: Strain{R},strainField :: StrainField{R},Index)
  strain.val += strainField.val[Index,:]
end
size(strainField :: StrainField) = size(strainField.val)[1:end-1]

type StressField{R <: Number} <: FluxSolutionTensorField
  val :: Array{R}
end
StressField(R,dims) = StressField(Array{R}(tuple(dims...,6)))
function get!{R}(stress :: Stress{R},stressField :: StressField{R},Index)
  stress.val = stressField.val[Index,:]
end
function set!{R}(stressField :: StressField{R},stress :: Stress{R},Index)
  stressField.val[Index,:] = stress.val
end
function add!{R}(stress :: Stress{R},stressField :: StressField{R},Index)
  stress.val += stressField.val[Index,:]
end
size(stressField :: StressField) = size(stressField.val)[1:end-1]


type DisplacementField{R <: Number} <: PrimarySolutionTensorField
  val :: Array{R}
end
DisplacementField(R,dims) =
DisplacementField(Array{R}(tuple(dims...,6)))
