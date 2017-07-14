module JuliaFFTHomogenization

import Base.*
import Base./
import Base.+
import Base.-
import Base.getindex
import Base.norm
import Base.promote_rule
import Base.convert
import Base.length
import Base.copy
import Base.copy!
import Base.setindex!
import Base.getindex
using ArgCheck
using MPAWL
import Base.size
using Polynomials

FFTW.set_num_threads(Sys.CPU_CORES)

include("Base/Types.jl")
include("Base/Transformations.jl")
include("Base/ApproximationMethod.jl")
include("Base/ConvergenceCriterion.jl")
include("Base/Solvers.jl")
include("Base/Solvers/BasicScheme.jl")

include("LinearElasticity/StiffnessParameters.jl")
include("LinearElasticity/StrainStressDisplacement.jl")
include("LinearElasticity/StrainStressDisplacementFields.jl")
include("LinearElasticity/StiffnessTensors.jl")
include("LinearElasticity/GreenOperators.jl")
include("LinearElasticity/Problem.jl")
include("LinearElasticity/ReferenceTensor.jl")

include("Misc/Misc.jl")
include("Misc/Ellipsoids.jl")

include("ExampleGeometries/ElasticityHashinEllipsoid.jl")
include("ExampleGeometries/Elasticity1DLaminate.jl")

end # module
