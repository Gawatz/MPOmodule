#
#	MPO structure
#
struct MPO{DT<:Union{localOp, AbstractArray}}
	bDim::Union{Int, Tuple{Int,Int}}
	Operator::Vector{DT}
	Op_index::Vector{Int}
	
	function MPO(bDim::Union{Tuple{Int,Int}, Int}, Operator::Union{Vector{<:localOp}, Vector{<:AbstractArray}}, Op_index::Vector{Int})
		allunique(Op_index) || throw(ArgumentError("operator indecies contains duplicate"))
		dim = typeof(bDim) == Int ? bDim*bDim : *(bDim...)
		dim >= maximum(Op_index) || throw(ArgumentError("operator indecies exceed MPO bond dimension"))


		return new{eltype(Operator)}(bDim, Operator, Op_index)
	end
end

# define MPOsparseT
const MPOsparseT = MPO{<:localOp}


import Base.eltype
eltype(mpo::MPOsparseT) = promote_type(eltype.(mpo.Operator)...)

getBondDim(mpo::MPO{DT}) where {DT<:Union{localOp, AbstractArray}} = typeof(mpo.bDim) == Int ? (mpo.bDim, mpo.bDim) : mpo.bDim

physDim(Op::localOp{DT, D}) where {DT<:Number, D} = D
physDim(mpo::MPO) = physDim(mpo.Operator[1])


import Base.getindex
function getindex(mpo::MPO{DT}, i::Int) where {DT<:Union{localOp, AbstractArray}}
	idx = findfirst(x->x == i, mpo.Op_index)


	if idx == nothing
		return nothing
	else
		return mpo.Operator[idx]
	end
end

getindex(mpo::MPO{DT}, i::Int, j::Int) where {DT<:Union{localOp, AbstractArray}} = getindex(mpo, LinearIndices((getBondDim(mpo)[1], getBondDim(mpo)[2]))[i,j])
getindex(mpo::MPO{DT}, i::Int, j::Colon) where {DT<:Union{localOp, AbstractArray}} = [getindex(mpo, x) for x in LinearIndices((getBondDim(mpo)[1], getBondDim(mpo)[2]))[i,:]]
getindex(mpo::MPO{DT}, i::Colon, j::Int) where {DT<:Union{localOp, AbstractArray}} = [getindex(mpo, x) for x in LinearIndices((getBondDim(mpo)[1], getBondDim(mpo)[2]))[:,j]]

"""
    eachLocalOp(mpo)

Iterates over all localOp and returns them in a tuple with the index position of the MPO.

# Arguments:
    - mpo:

return:
	
	!!!ADD example!!!

"""
function eachLocalOp(mpo::MPO{DT}) where {DT<:Union{localOp, AbstractArray}} 
	idx = CartesianIndices(getBondDim(mpo))[mpo.Op_index]
	idx = [(x[1],x[2]) for x in idx]

	return zip(mpo.Operator, idx)
end

#
#	this should be changed
#
"""
	iteration over all terms in MPO
"""
function eachoperation(mpo::MPOsparseT)
	tmp = []
	@inbounds for x in eachLocalOp(mpo)
		for y in eachoperation(x[1])
			push!(tmp, [y[1], Tuple(y[2:3]), x[2]])
		end
	
	end
	return tmp
end


import Base.show
function Base.show(io::IO, mpo::MPOsparseT)
	
	bDim = getBondDim(mpo)
	idx = [CartesianIndices(bDim)[x] for x in mpo.Op_index]
	B = zeros(Int, bDim)
	B[idx] .= [1:length(mpo.Op_index)...]

	print(io, "operators in MPO:\n $B \n")

	for (i,x) in enumerate(mpo.Operator)
		
		print(io, " operator: $i --> $(Array(x))\n")
	
	end
end

function getMPOTensor(mpo::MPO{DT}) where {DT<:Union{localOp, AbstractArray}}
	dims = getBondDim(mpo)
	phys_dim = physDim(mpo)
	data_type = eltype(mpo)
	MPOTensor = zeros(data_type, (phys_dim, dims..., phys_dim))

	for (op, (i,j)) in eachLocalOp(mpo)
		MPOTensor[:, i, j, :] = Array(op)
	end

	return MPOTensor
end

function getMPOFromTensor(mpoTensor::Array{DT, 4}) where {DT<:Number}
	bond_dim = size(mpoTensor)[2:3]

	operator = localOp[]
	idx = Int[]
	for col in 1:bond_dim[2]
		for row in 1:bond_dim[1]

			op = mpoTensor[:,row,col,:]

			if iszero(op)!=true
				push!(operator, localOp(op))
				push!(idx, LinearIndices(bond_dim)[row, col])
			end
		end
	end

	mpo = MPO(bond_dim, operator, idx)
	return mpo
end

#
#	some MPO algebra
#
import Base.(*)
import Base.(+)

(*)(a::Number, MPOvec::Vector{<:MPO{DT}}) where {DT<:Union{localOp, AbstractArray}} = Vector{MPO{DT}}([MPO(MPOvec[1].bDim, a .* MPOvec[1].Operator, MPOvec[1].Op_index),MPOvec[2:end]...])
(+)(mpo_a::Vector{<:MPO{DT}}, mpo_b::Vector{<:MPO{DT}}) where {DT<:Union{localOp, AbstractArray{<:Number}}} = add_MPO(mpo_a, mpo_b)

function add_MPO(MPOvecA::Vector{<:MPO{DT}}, MPOvecB::Vector{<:MPO{DT}}) where {DT<:Union{localOp, AbstractArray{<:Number}}}
	@assert length(MPOvecA) == length(MPOvecB) "you are trying to add
					to MPO with different length"

	newMPOvec = Vector{MPO{DT}}([])				
	
	for i in 1:size(MPOvecA)[1]
		
		MPOA = MPOvecA[i]
		MPOB = MPOvecB[i]
		
		bDim_A = MPOA.bDim
		bDim_B = MPOB.bDim
		
		bDim_A = typeof(bDim_A) == Int ? (bDim_A, bDim_A) : bDim_A
		bDim_B = typeof(bDim_B) == Int ? (bDim_B, bDim_B) : bDim_B

		# take care of last and first site 
		new_bin = i != 1 ? bDim_A[1]+bDim_B[1] : 1
		new_bout = i != size(MPOvecA)[1] ? bDim_A[2]+bDim_B[2] : 1

		OP_vec_A = MPOA.Operator
		OP_idx_A = MPOA.Op_index
		
		OP_vec_B = MPOB.Operator
		OP_idx_B = MPOB.Op_index

		#
		#	bring into block diagonal form
		#
		
		# adjust OP_idx_A
		idx = CartesianIndices((bDim_A[2],bDim_A[1]))[OP_idx_A]
		OP_idx_A = LinearIndices((new_bout,bDim_A[1]))[idx]


		# adjust Op_idx_B
		idx = CartesianIndices((bDim_B[2],bDim_B[1]))[OP_idx_B]
		idx = [Array([Tuple(x)...]).+[new_bout-bDim_B[2],new_bin-bDim_B[1]] for x in idx]	
		OP_idx_B = [LinearIndices((new_bout, new_bin))[x...] for x in idx]
		new_Op_vec = Vector{DT}([OP_vec_A..., OP_vec_B...])
		new_Op_idx = Vector{Int}([OP_idx_A..., OP_idx_B...])
		
		newMPO = MPO((new_bin, new_bout), new_Op_vec, new_Op_idx)
		push!(newMPOvec, newMPO)
	end

	return newMPOvec
end

function (*)(MPOvecA::Vector{<:MPOsparseT}, MPOvecB::Vector{<:MPOsparseT})

	@assert length(MPOvecA) == length(MPOvecB) "you are trying to mulitply
						    to MPO with different length"

	new_MPOvec = MPO{localOp}[]
	for i in 1:length(MPOvecA)

		mpoA = getMPOTensor(MPOvecA[i])
		mpoB = getMPOTensor(MPOvecB[i])

		@tensor new_mpo[α, α′, d, d′, β, β′] := mpoA[d, α, β, γ]*mpoB[γ, α′, β′, d′]
		
		new_mpo = reshape(new_mpo, *(size(new_mpo)[1:2]...), size(new_mpo)[3], size(new_mpo)[4], *(size(new_mpo)[5:end]...))
		new_mpo = permutedims(new_mpo, [2,1,4,3])
		new_mpo = getMPOFromTensor(new_mpo)
		push!(new_MPOvec, new_mpo)
	end
	
	return new_MPOvec
end

#
#	multiplication of MPO with Array 
#	(multiplication acts as matrix multiplication only on the bond dimension)
#


function (*)(mpo::MPO{localOp{DT,D}}, a::AbstractArray{DT, 2}) where {DT<:Number, D}
	tmp = []

	#multiplication from the right
	@inbounds for col in 1:size(a)[2]
		for row in 1:getBondDim(mpo)[1]
			array_col = @view a[:,col]  # get column
			idx = findall(x->x!=0, array_col) # find nz elements in column
			if !isempty(idx) 
				mpo_row = mpo[row,:]  # eventually replace with a view
				idx_nothing = findall(x->typeof(x) == Nothing, mpo_row) #check if there is a localOp
				idx = setdiff!(idx,idx_nothing) #set difference of idx: make sure that element of idx is not in idx_nothing
				if !isempty(idx)  
					mpo_row = mpo_row[idx] #[mpo[row,i] for i in idx] # get elements which are at nz column index
					coef = array_col[idx] # [array_col[i] for i in idx]
					new_localOp = sum(coef.*mpo_row) #obtain actual new localOp
					push!(tmp, (new_localOp, (row,col)))
				end
			end
		end

	end

	# construct new OP
	Operator = [x[1] for x in tmp]
	Op_index = [LinearIndices((getBondDim(mpo)[1], size(a)[2]))[x[2]...] for x in tmp]


	idx = findall(x->x!=nothing, Operator)
	Operator = Vector{localOp{DT,D}}(Operator[idx])
	Op_index = Op_index[idx]
	# this would just be necessary if we would loop first over row
	#idx = sortperm(Op_index)
	#Op_index = Op_index[idx]
	#Operator = Operator[idx]


	return MPO((getBondDim(mpo)[1], size(a)[2]), Operator, Op_index)
end

function (*)(a::AbstractArray{DT, 2}, mpo::MPO{localOp{DT,D}}) where {DT<:Number, D}
	tmp = []

	#multiplication from the left
	@inbounds for col in 1:getBondDim(mpo)[2]
		for row in 1:size(a)[1]

			array_row = @view a[row,:]  # get row
			idx = findall(x->x!=0, array_row) # find nz elements in row
			if !isempty(idx) 
				mpo_col = mpo[:,col]
				idx_nothing = findall(x->typeof(x) == Nothing, mpo_col) #check if there is a localOp
				idx = setdiff!(idx,idx_nothing)
				if !isempty(idx) 
					mpo_col = mpo_col[idx] #[mpo[i, col] for i in idx] # get elements which are at nz column index
					coef = array_row[idx] #[array_row[i] for i in idx]
					new_localOp = sum(coef.*mpo_col) #obtain actual new localOp
					push!(tmp, (new_localOp, (row,col)))
				end
			end
		end

	end

	# construct new OP
	Operator = [x[1] for x in tmp]
	Op_index = [LinearIndices((size(a)[1], getBondDim(mpo)[2]))[x[2]...] for x in tmp]
	
	idx = findall(x->x!=nothing, Operator)
	Operator = Vector{localOp{DT,D}}(Operator[idx])
	Op_index = Op_index[idx]
	
	# this would just be necessary if we would loop first over row
	#idx = sortperm(Op_index)
	#Op_index = Op_index[idx]
	#Operator = Operator[idx]


	return MPO((size(a)[1], getBondDim(mpo)[2]), Operator, Op_index)
end

#
#	function for deparalization of MPO (shrinks complexity)
#
#	so far only for localOp structure not defined for Array operator
#

function checkPara(a, b)
	@assert length(a) == length(b)

	parallel = false
	c = 0
	if typeof.(a) != typeof.(b)
		return (false, 0.0)
	end

	idx = findall(x->typeof(x) != Nothing, a)
	#tmp = [parallelOp(a[x], b[x]) for x in idx]
	tmp = [a[x]∥b[x] for x in idx]  # check if operator are parallel  
	if length(unique([x[2] for x in tmp])) != 1 # check if the scale between them are identical
		return (false, 0.0)
	end


	return (tmp[1][1], tmp[1][2])
end


function findParaRow(mpo::MPO{localOp{DT,D}}) where {DT<:Number, D}
	a = zeros(ComplexF64, getBondDim(mpo)[1])
	a[1] = 1.0
	T = [a]

	mpo_row = mpo[1,:]
	K = Any[mpo_row]

	for row_i in length(K)+1:getBondDim(mpo)[1]
		para = false

		for k_idx in 1:length(K)

			mpo_row = mpo[row_i, :]
			k = K[k_idx]

			para, scale = checkPara(mpo_row, k)

			if para == true

				T[k_idx][row_i] = scale
				break
			end
		end

		if para != true
			push!(K, mpo_row)
			a = zeros(ComplexF64, getBondDim(mpo)[1])
			a[row_i] = 1.0
			push!(T, a)
		end

	end

	# recat K into mpo
	
	K = hcat(K...)
	K = K[Transpose(LinearIndices(size(K)))]
	idx = findall(x->typeof(x) != Nothing, K)
	op_idx = LinearIndices(size(K))[idx]
	K = MPO(size(K), Vector{localOp{DT,D}}(K[idx]), op_idx)

	
	return K, hcat(T...)
end

function findParaCol(mpo::MPO{localOp{DT,D}}) where {DT<:Number, D}	
	a = zeros(ComplexF64, getBondDim(mpo)[2])
	a[1] = 1.0
	T = [a]
	mpo_col = mpo[:,1]
	K = Any[mpo_col]

	for col_i in length(K)+1:getBondDim(mpo)[2]
		para = false

		for k_idx in 1:length(K)
			mpo_col = mpo[:,col_i]
			k = K[k_idx]

			para, scale = checkPara(mpo_col, k)
			if para == true

				T[k_idx][col_i] = scale
				break
			end

		end

		if para != true
			
			push!(K, mpo_col)
			a = zeros(ComplexF64, getBondDim(mpo)[2])
			a[col_i] = 1.0
			push!(T,a)



		end

	end


	# recat K into mpo
	K = hcat(K...)
	idx = findall(x->typeof(x) != Nothing, K)
	op_idx = LinearIndices(size(K))[idx]
	K = MPO(size(K), Vector{localOp{DT,D}}(K[idx]), op_idx)

	return K, Transpose(hcat(T...))
end

#
#	deparallizes mpo-chain
#
function depara!(MPOvec::Vector{MPO{localOp{DT,D}}}) where {DT<:Number, D}
	for site in 1:length(MPOvec)
		#@show "dePara col", site
		mpo = MPOvec[site]

		mpo, T = findParaCol(mpo)
		MPOvec[site] = mpo

		if site+1 <= length(MPOvec)
			mpo_next = MPOvec[site+1]
			mpo_next = T*mpo_next
			MPOvec[site+1] = mpo_next
		end
	end
	
	for site in length(MPOvec):-1:1
		#@show "dePara row", site
		mpo = MPOvec[site]
		mpo, T = findParaRow(mpo)
		MPOvec[site] = mpo
		if site-1 > 1
			mpo_pre = MPOvec[site-1]
			mpo_pre = mpo_pre*T
			MPOvec[site-1] = mpo_pre
		end
	end


	for site in 1:length(MPOvec)
		#@show "dePara col", site	
		mpo = MPOvec[site]

		mpo, T = findParaCol(mpo)
		MPOvec[site] = mpo

		if site+1 <= length(MPOvec)
			mpo_next = MPOvec[site+1]
			mpo_next = T*mpo_next
			MPOvec[site+1] = mpo_next
		end
	end

	return MPOvec
end
