struct OpString 
	idx::String

	function OpString(a::String)
		return new(a)
	end

	function OpString(a::Vector{String})

		if length(a)>1
			res = OpString(a[1]) + OpString(a[2:end])
		else
			res = OpString(a[1])
		end

		return new(res.idx)

	end
end

#
# overwritting base operators for OpString matrix multiplication 
#

# in case either a or b is empty "" creat new string by just combining otherwise add a "→" between them 
import Base.+
(+)(a::OpString, b::OpString) = a.idx=="" || b.idx=="" ? OpString(string(a.idx, b.idx)) : OpString(string(a.idx,"→ ",b.idx))

# in case either a or b is empty "" the new string is also empty otherwise split operator string and add a ⋅ inbetween 
# note that a an b OpStrings have to have same length
import Base.*
(*)(a::OpString, b::OpString) = a.idx == "" || b.idx == "" ? OpString("") : OpString(string.(split(a.idx,"→ "),"⋅", split(b.idx,"→ ")))

import Base.zero
zero(a::OpString) = OpString("")

function StringArray_to_OpStringArray(A::AbstractArray{<:String})
	
	B = Array{OpString}
	B = fill(OpString(""),size(A)...)
	B[:,:] = OpString.(A[:,:])
	
	return B
end

function MPO_to_StringArray(mpo_A::MPO)
	bDim = getBondDim(mpo_A)
	bin = 0
	bout = 0
	try 
		bin = bDim[1]
		bout = bDim[2]
	catch
		bin = bDim
		bout = bDim
	end
	
	B = Array{String}
	B = fill("", bin, bout)
	Op_index = mpo_A.Op_index
	
	Idx = LinearIndices(fill(1,bin,bout))
	mpo_idx = [findfirst(x-> x == idx, Idx) for idx in Op_index]


	B[mpo_idx] .= [string(i) for i in 1:size(mpo_idx)[1]]

	return B
end

MPO_to_OpString(mpoA::MPO) = StringArray_to_OpStringArray(MPO_to_StringArray(mpoA))

# matrix multiplication

function stringM_mul(A::AbstractArray{OpString}, B::AbstractArray{OpString})	
	return A*B
end

stringM_mul(A::AbstractArray{String}, B::AbstractArray{String}) = stringM_mul(StringArray_to_OpStringArray(A), StringArray_to_OpStringArray(B))

"""
   convert_OpStringArrayToTuple(A::AbstractArray{<:String})

converts and Array respresenting the MPO product as strings to an Array containing
tuples where each int corresponds to the operator index of the local MPO

#


"""
function convert_OpStringArrayToTuple(A::AbstractArray{OpString})
	bin, bout = size(A)
	B = Array{Any,2}
	B = fill([],bin,bout)
	for i in 1:bin
		row = []
		for j in 1:bout			
			try
				#decompose string
				terms = split(A[i,j].idx,"→ ")
				terms = [Tuple(parse.(Int, split(x,"⋅"))) for x in terms]
				
				B[i,j] = terms
			catch
				continue
			end

		end
	end

	return B
end

"""
    operatorString(MPOvec; op_String = nothing)

Returns all operator strings contained in an mpo chain in form of a list of integer tuples.
Example:


# Arguments:
 

"""
function stringMPO(MPOvec::Vector{MPO{localOp{DT,D}}}; op_String::Union{Array{OpString,2},Nothing} = nothing) where {DT<:Number, D}
	mpo = MPOvec[1]
	B = MPO_to_OpString(mpo)

	op_String = op_String == nothing ? B : stringM_mul(op_String,B)
	
	if size(MPOvec[2:end])[1] > 0
		stringMPO(MPOvec[2:end], op_String = op_String)
	else
		return convert_OpStringArrayToTuple(op_String) 
	
	end
end
