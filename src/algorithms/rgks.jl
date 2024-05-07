using LinearAlgebra
using Random

include("../utilities/returnstructs.jl")

"""
	rgks([rng], A, k; oversamp = 0, minimal = false, orthonormal = true)

Compute an approximate factorization `A = C*R` where `C` consists of `k` skeleton columns from `A`. Choose columns using the procedure of Golub, Klemma, and Stewart (1976) on a sketch of `A` computed using oversampling `oversamp`. If `minimal == true` then only the indices of the skeleton columns are computed. If `minimal == false` then `C` and `R` are returned along with the indices of the skeleton columns. If `orthonormal == true` then the columns of `C` are orthonormalized; otherwise its columns are the corresponding columns of `A` unmodified. All internal randomness is drawn from `rng`.

NOTE: running `rgks` with `orthonormal = false` may be numerically unstable when `k` exceeds the numerical rank of `A`.
"""
function rgks(rng::AbstractRNG, A::Matrix, k::Integer ;
				oversamp::Integer = 0,
				minimal::Bool = false,
				orthonormal::Bool = true)
	
	if((k < 0) || (k > min(size(A, 1), size(A, 2))))
		throw(ArgumentError("the target rank ("*string(k)*") must be at least 1 and at most "*string(min(size(A, 1), size(A, 2)))))
	
	elseif(oversamp < 0)
		throw(ArgumentError("the oversampling parameter ("*string(oversamp)*") must be nonnegative"))
	end
	
	A_sk = A*randn(rng, size(A, 2), k + oversamp)
	Q = Matrix(qr(A_sk).Q)
	
	svdobj = svd(Q'*A, alg = LinearAlgebra.QRIteration())	
	Vt = svdobj.Vt[1:k, :]
	qrobj = qr(Vt, ColumnNorm())
	perm = qrobj.p[1:k]
	
	if(minimal)
		return perm
	end
	
	C = orthonormal ? Matrix(qr(A[:, perm]).Q) : A[:, perm]
	Cp = orthonormal ? C' : pinv(C)
	R = Cp*A
	
	return orthonormal ? OrthoSkeletalDecomp(perm, C, R) : SkeletalDecomp(perm, C, R)
end

function rgks(A::Matrix, k::Integer ;
				oversamp::Integer = 0,
				minimal::Bool = false,
				orthonormal::Bool = true)
				
	return rgks(Random.default_rng(), A, k, oversamp = oversamp, minimal = minimal, orthonormal = orthonormal)
end
