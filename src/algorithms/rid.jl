using LinearAlgebra
using Random

include("../utilities/returnstructs.jl")

"""
	rid([rng], A, k; oversamp = 0, minimal = false, orthonormal = true)

Compute an approximate factorization `A = C*X` where `C` consists of `k` skeleton columns from `A`. Choose columns using Businger-Golub QRCP on a sketch of `A` computed using oversampling `oversamp`. If `minimal == true` then only the indices of the skeleton columns are computed. If `minimal == false` then `C` and `X` are returned along with the indices of the skeleton columns. If `orthonormal == true` then the columns of `C` are orthonormalized; otherwise its columns are the corresponding columns of `A` unmodified. All internal randomness is drawn from `rng`.

NOTE: running `rid` with `orthonormal = false` may be numerically unstable when `k` exceeds the numerical rank of `A`.
"""
function rid(rng::AbstractRNG, A::Matrix, k::Integer ;
				oversamp::Integer = 0,
				minimal::Bool = false,
				orthonormal::Bool = true)
	
	if((k < 0) || (k > min(size(A, 1), size(A, 2))))
		throw(ArgumentError("the target rank ("*string(k)*") must be at least 1 and at most "*string(min(size(A, 1), size(A, 2)))))
	
	elseif(oversamp < 0)
		throw(ArgumentError("the oversampling parameter ("*string(oversamp)*") must be nonnegative"))
	end
	
	A_sk = randn(rng, k + oversamp, size(A, 1))*A
	qrobj = qr(A_sk, ColumnNorm())
	perm = qrobj.p[1:k]
	
	if(minimal)
		return perm
	end
	
	C = orthonormal ? Matrix(qr(A[:, perm]).Q) : A[:, perm]
	Cp = orthonormal ? C' : pinv(C)
	R = Cp*A
	
	return orthonormal ? OrthoSkeletalDecomp(perm, C, R) : SkeletalDecomp(perm, C, R)
end

function rid(A::Matrix, k::Integer ;
				oversamp::Integer = 0,
				minimal::Bool = false,
				orthonormal::Bool = true)
	
	return rid(Random.default_rng(), A, k, oversamp = oversamp, minimal = minimal, orthonormal = orthonormal)
end
