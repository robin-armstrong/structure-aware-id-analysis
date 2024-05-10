using LinearAlgebra
using StatsBase
using Random

include("../utilities/returnstructs.jl")

"""
	levg([rng], A, k; oversamp = 0, minimal = false, orthonormal = true, leverage_scores)

Compute an approximate factorization `A = C*X` where `C` consists of `k + oversamp` skeleton columns from `A`. Choose columns randomly using leverage score sampling with replacement. Optionally, leverage scores can be precomputed and provided through the `leverage_score` argument. If `minimal == true` then only the indices of the skeleton columns are computed. If `minimal == false` then `C` and `X` are returned along with the indices of the skeleton columns. If `orthonormal == true` then the columns of `C` are orthonormalized; otherwise its columns are the corresponding columns of `A` unmodified. All internal randomness is drawn from `rng`.

NOTE: running `levg` with `orthonormal = false` may be numerically unstable when `k + oversamp` exceeds the numerical rank of `A`.
"""
function levg(rng::AbstractRNG, A::Matrix, k::Integer ;
				oversamp::Integer = 0,
				minimal::Bool = false,
				orthonormal::Bool = true,
				leverage_scores::Union{Vector{T}, Nothing} = nothing) where T <: Real
	if(oversamp < 0)
		throw(ArgumentError("the oversampling parameter ("*string(oversamp)*") must be nonnegative"))
	end
	
	lscores_provided = !isnothing(leverage_scores)

	if(lscores_provided && (minimum(leverage_scores.*(1 .- leverage_scores)) < 0))
		throw(ArgumentError("all leverage scores must be between 0 and 1"))
	end
	
	lscores = lscores_provided ? leverage_scores : zeros(size(A, 2))
	
	if(!lscores_provided)
		Vt = svd(A).Vt[1:k, :]
		
		for j = 1:size(A, 2)
			lscores[j] = norm(Vt[:, j])^2
		end
	end
	
	cols = sample(rng, 1:size(A, 2), Weights(lscores), k + oversamp)
	
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
				orthonormal::Bool = true,
				leverage_scores::Union{Vector{T}, Nothing} = nothing) where T <: Real
	
	return levg(Random.default_rng(), A, k, oversamp = oversamp, minimal = minimal, orthonormal = orthonormal, leverage_scores = leverage_scores)
end
