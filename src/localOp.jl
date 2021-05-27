#
#	Local Operator is a sparse operator structure only storing in and outgoing indices 
#	and the coefficient associated with eacht in-out index pair:
#
#
#	σʸ = (1,2),(2,1) index pairs with coeffcients -1.0im and 1.0im
#
"""
    localOp

Local Operator is a sparse operator structure only storing in and outgoing indices 
and the coefficient associated with each in-out index pair:

σʸ = (1,2),(2,1) index pairs with coeffcients -1.0im and 1.0im

To Do: !!! Essentially could be replaced by sparse array with coo indexing.
"""
struct localOp{DT<:Number, D}
	c::Vector{DT}
	in::Vector{Int}
 	out::Vector{Int}

 	function localOp(A::AbstractArray{<:Number,2})
		idx = findall(x->x!=0, A)
		DT = eltype(A)

		in = [x[1] for x in idx]
		out = [x[2] for x in idx]
		c = A[idx]

		return new{DT, size(A)[1]}(c, in, out)
	end

	function localOp(c::Vector{<:Number}, in::Vector{Int}, out::Vector{Int}, D::Int)
		# check dimensions
		D == maximum(in) || throw(ArgumentError("ingoing index exceeds the local operator dimension"))
		D == maximum(out) || throw(ArgumentError("outgoing index exceeds the local operator dimension"))
		length(c) <= D*D || throw(ArgumentError("specified entires of local operator exceeds dimension"))
		# check if each element is only referenced ones 
		allunique(in) || throw(ArgumentError("ingoing indicies contains duplicate"))
		allunique(out) || throw(ArgumentError("outgoing indicies contains duplicate"))


		return new{eltype(c), D}(c, in, out)
	end
end



physdim(op::localOp{DT, D}) where{DT<:Number, D} = D

import Base.convert
function convert(::Type{Array{DT, 2}}, op::localOp{DT, D}) where {DT<:Number, D}
	OP = zeros(DT, D, D)
	@inbounds for x in eachoperation(op)
		OP[x[2],x[3]] = x[1]
	end

	return OP
end

import Base.Array
Array(op::localOp{DT, D}) where {DT<:Number, D} = convert(Array{DT,2}, op)

import Base.show
function Base.show(io::IO, op::localOp)	
	print(io, "local operator: ", Array(op))
end

import Base.*
*(a::Number, b::localOp{DT,D}) where {DT<:Number, D}= localOp(a.*b.c, b.in, b.out, D)
*(b::localOp{DT,D}, a::Number) where {DT<:Number, D}= localOp(a.*b.c, b.in, b.out, D)


import Base.+
function +(a::localOp{DT1, D}, b::localOp{DT2, D}) where {DT1<:Number, DT2<:Number, D} 

end

import Base.eltype
eltype(op::localOp{DT,D}) where {DT<:Number, D} = DT


# get list with tuples (coef, in_index, out_index)
function eachoperation(a::localOp)
	return zip(a.c, a.in, a.out)
end

# check if two localOP are parallel 
function (∥)(a::localOp, b::localOp)	
	parallel = false
	scaling = 0
	if a.in == b.in && a.out == b.out	
		if iszero(diff(a.c./b.c))
			parallel = true
			scaling = (a.c./b.c)[1]
		end
	end

	return parallel, scaling
end
