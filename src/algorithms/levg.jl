using LinearAlgebra
using StatsBase
using Random

include("../utilities/returnstructs.jl")

"""
	levg([rng], A, k; oversamp = 0, minimal = false, orthonormal = true)

Compute an approximate factorization `A = C*R` where `C` consists of `k + oversamp` skeleton columns from `A`. Choose columns randomly using leverage score sampling with replacement. If `minimal == true` then only the indices of the skeleton columns are computed. If `minimal == false` then `C` and `R` are returned along with the indices of the skeleton columns. If `orthonormal == true` then the columns of `C` are orthonormalized; otherwise its columns are the corresponding columns of `A` unmodified. All internal randomness is drawn from `rng`.

NOTE: running `levg` with `orthonormal = false` may be numerically unstable when `k + oversamp` exceeds the numerical rank of `A`.
"""
function levg(rng::AbstractRNG, A::Matrix, k::Integer ;
				oversamp::Integer = 0,
				minimal::Bool = false,
				orthonormal::Bool = true)
	if(oversamp < 0)
		throw(ArgumentError("the oversampling parameter ("*string(oversamp)*") must be nonnegative"))
	end
	
	Vt = svd(A).Vt[1:k, :]
	levg = Vector{Float64}(undef, size(A, 2))
	
	for j = 1:size(A, 2)
		levg[j] = norm(Vt[:, j])^2
	end
	
	cols = sample(rng, 1:size(A, 2), Weights(levg), k + oversamp)
	
	if(minimal)
		return cols
	end
	
	C = orthonormal ? Matrix(qr(A[:, cols]).Q) : A[:, cols]
	Cp = orthonormal ? C' : pinv(C)
	R = Cp*A
	
	return orthonormal ? OrthoSkeletalDecomp(cols, C, R) : SkeletalDecomp(cols, C, R)
end

function levg(A::Matrix, k::Integer ;
				oversamp::Integer = 0,
				minimal::Bool = false,
				orthonormal::Bool = true)
	
	return levg(Random.default_rng(), A, k, oversamp = oversamp, minimal = minimal, orthonormal = orthonormal)
end
