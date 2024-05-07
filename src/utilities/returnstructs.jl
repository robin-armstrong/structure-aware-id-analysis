# struct, iterator, and show method for randomized SVD

"""
	RSVD(U::Matrix, S::Vector, Vt::Matrix)

Return type for `rsvd`, the randomized singular value decomposition. The `S` vector contains estimated singular values in decreasing order, while `U` and `Vt` are the estimated orthogonal factors.
"""
struct RSVD
	U::Matrix
	S::Vector
	Vt::Matrix
end

Base.iterate(F::RSVD) = (F.U, Val(:S))
Base.iterate(F::RSVD, ::Val{:S}) = (F.S, Val(:Vt))
Base.iterate(F::RSVD, ::Val{:Vt}) = (F.Vt, Val(:done))
Base.iterate(F::RSVD, ::Val{:done}) = nothing

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, F::RSVD)
	summary(io, F); println(io)
	println(io, "U factor:")
	show(io, mime, F.U)
	println(io, "\nsingular values:")
	show(io, mime, F.S)
	println(io, "\nVt factor:")
	show(io, mime, F.Vt)
end

# struct, iterator, and show method for randomized skeletal decompositions

"""
	SkeletalDecomp(skel::Vector, C::Matrix, R::Matrix)

Return type for low-rank methods that compute interpolative decompositions. The vector `skel` contains the indices of the skeleton columns chosen used for the decomposition, while `C` and `R` are the factors of the decomposition.
"""
struct SkeletalDecomp
	skel::Vector
	C::Matrix
	R::Matrix
end

Base.iterate(F::SkeletalDecomp) = (F.skel, Val(:C))
Base.iterate(F::SkeletalDecomp, ::Val{:C}) = (F.C, Val(:R))
Base.iterate(F::SkeletalDecomp, ::Val{:R}) = (F.R, Val(:done))
Base.iterate(F::SkeletalDecomp, ::Val{:done}) = nothing

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, F::SkeletalDecomp)
	summary(io, F); println(io)
	println(io, "skeleton column indices:")
	show(io, mime, F.skel)
	println(io, "\nskeleton columns:")
	show(io, mime, F.C)
	println(io, "\ninterpolation matrix:")
	show(io, mime, F.R)
end

# struct, iterator, and show method for orthonormalized skeletal decompositions

"""
	OrthoSkeletalDecomp(skel::Vector, C::Matrix, R::Matrix)

Return type for low-rank methods that compute interpolative decompositions with orthonormalized skeleton columns. The vector `skel` contains the indices of the skeleton columns chosen used for the decomposition, while `C` and `R` are the factors of the decomposition.
"""
struct OrthoSkeletalDecomp
	skel::Vector
	C::Matrix
	R::Matrix
end

Base.iterate(F::OrthoSkeletalDecomp) = (F.skel, Val(:C))
Base.iterate(F::OrthoSkeletalDecomp, ::Val{:C}) = (F.C, Val(:R))
Base.iterate(F::OrthoSkeletalDecomp, ::Val{:R}) = (F.R, Val(:done))
Base.iterate(F::OrthoSkeletalDecomp, ::Val{:done}) = nothing

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, F::OrthoSkeletalDecomp)
	summary(io, F); println(io)
	println(io, "skeleton column indices:")
	show(io, mime, F.skel)
	println(io, "\northonormalized skeleton columns:")
	show(io, mime, F.C)
	println(io, "\ninterpolation matrix:")
	show(io, mime, F.R)
end
