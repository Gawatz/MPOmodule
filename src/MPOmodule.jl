module MPOmodule
using LinearAlgebra, TensorOperations

export MPO, stringMPO, localOp
export getBondDim, eachLocalOp
export eachoperation, eachOp, get_MPOTensor
export (*), (+), (∥), depara!
export MPO_to_StringArray
export MPOsparseT

include("localOp.jl")
include("mpoStruct.jl")
include("operatorString.jl")

end # module
