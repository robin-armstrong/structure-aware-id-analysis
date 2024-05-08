using LinearAlgebra
using Random

include("../utilities/returnstructs.jl")

"""
	rsvd([rng], A, k; oversamp = 0, power = 0, minimal = false)

Compute an approximate factorization `A = U*diagm(S)*Vt` where `U` consists of `k` orthonormal columns (left singular vector estimates), `Vt` consists of `k` orthonormal rows (right singular vector estimates), and `S` is a vector of length `k` (singular value estimates) using the algorithm of Halko, Martinsson, and Tropp (2011) with oversampling `oversamp` and `power` steps of power iteration. If `minimal == true` then `U` is not computed and only the singular value estmiates are returned. If  `minimal == false` then `U`, `S`, and `Vt` are returned. All internal randomness is drawn from `rng`.
"""
function rsvd(rng::AbstractRNG, A::Matrix, k::Integer ; 
				oversamp::Integer = 0,
				power::Integer = 0,
				minimal::Bool = false)
	
	if(power < 0)
		throw(ErrorException("power iteration parameter ("*string(power)*") must be nonnegative"))
		
	elseif((k < 0) || (k > min(size(A, 1), size(A, 2))))
		throw(ArgumentError("the target rank ("*string(k)*") must be at least 1 and at most "*string(min(size(A, 1), size(A, 2)))))
	
	elseif(oversamp < 0)
		throw(ArgumentError("the oversampling parameter ("*string(oversamp)*") must be nonnegative"))
	end
	
	A_sk = A*randn(rng, size(A, 2), k + oversamp)
	Q = Matrix(qr(A_sk).Q)
	
	for i = 1:power
		Q = Matrix(qr(A*(A'*Q)).Q)
	end
	
	svdobj = svd(Q'*A, alg = LinearAlgebra.QRIteration())
	
	if(minimal)
		return svdobj.S[1:k]
	end
	
	return RSVD(Q*svdobj.U[:, 1:k], svdobj.S[1:k], svdobj.Vt[1:k, :])
end

function rsvd(A::Matrix, k::Integer ; 
				oversamp::Integer = 0,
				power::Integer = 0,
				minimal::Bool = false)
	
	return rsvd(Random.default_rng(), A, k, oversamp = oversamp, power = power, minimal = minimal)
end
