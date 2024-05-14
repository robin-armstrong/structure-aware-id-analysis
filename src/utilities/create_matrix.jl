using LinearAlgebra
using Random
using SpecialFunctions
using Hadamard

"""
	SpectrumObject(levels::Vector{Float64}, decayRanges::Vector{Int64}, decayTypes::Vector{String})

Create a new `SpectrumObject` to be passed to the `createMatrix` method, specifying the singular spectrum of the matrix to be created. Entries of `levels` specify the levels of the singular spectrum in decreasing order. Transitions between levels occur in the index ranges given by `decayRanges`, a vector of increasing integers. Each consecutive pair of integers in `decayRanges` specifies a range of indices in which decay occurs. Consequently, `length(decayRanges)` must equal `2*(length(levels) - 1)`. The type of decay in each decay window is specified by the corresponding entry in `decayTypes`, a vector of strings which are either 'smoothgap' or 'exponential'. Consequently, `2*length(decayTypes)` must equal `length(decayRanges)`.
"""
struct SpectrumObject
	levels::Vector{Float64}
	decayRanges::Vector{Int64}
	decayTypes::Vector{String}
	
	function SpectrumObject(levels::Vector{Float64}, decayRanges::Vector{Int64}, decayTypes::Vector{String})
		if(minimum(levels) <= 0)
			throw(ArgumentError("entries of levels vector must be strictly positive"))
		end
		
		for i = 2:length(levels)
			if(levels[i] >= levels[i - 1])
				throw(ArgumentError("entries of levels vector must be strictly decreasing"))
			end
		end
		
		if(minimum(decayRanges) <= 0)
			throw(ArgumentError("entries of decayRanges must be positive"))
		end
		
		for i = 2:length(decayRanges)
			if(decayRanges[i] <= decayRanges[i - 1])
				throw(ArgumentError("entries of decayRanges must be strictly increasing"))
			end
		end
		
		if(length(decayRanges) != 2*(length(levels) - 1))
			throw(ArgumentError("length(decayRanges) must equal 2*(length(levels) - 1)"))
		end
		
		for i = 1:length(decayTypes)
			if((decayTypes[i] != "smoothgap") && (decayTypes[i] != "exponential"))
				throw(ArgumentError("unrecognized decay type '"*decayTypes[i]*"', decay types must be either 'smoothgap' or 'exponential'"))
			end
		end
		
		if(2*length(decayTypes) != length(decayRanges))
			throw(ArgumentError("2*length(decayTypes) must equal length(decayRanges)"))
		end
		
		return new(levels, decayRanges, decayTypes)
	end
end

"""
	create_matrix([rng], dim, spectrum; coher, minimal = false)

Randomly generate a square matrix of dimension `dim` whose singular values decay as specified by `spectrum`, a variable of type `SpectrumObject`. Left and right singular vectors are orthogonalizations of Gaussian random matrices, unless the value of `coher` is set. In this case, `coher` is a real number between 0 and 1 controlling the level of coherence, and right singular vectors are given as

	V = orth(coher*P + (1 - coher)*H)

where `P` is a permutation matrix, `H` is a normalized Hadamard matrix, and `orth` denotes orthogonalization via a polar decomposition. Factors `U`, `S`, and `Vt` of a singular value decomposition are returned. If `minimal = true` then only `S` is computed and returned. Randomness is drawn from `rng`.
"""
function create_matrix(rng::AbstractRNG, dim::Int, spectrum::SpectrumObject; coher::Union{Real, Nothing} = nothing, minimal::Bool = false)
	# checking compatibility of arguments
	if(spectrum.decayRanges[end] > dim)
		throw(ArgumentError("spectrum object has decay indices that exceed the matrix dimension"))

	elseif(coher != nothing)
		if((coher > 1.0) || (coher < 0.0))
			throw(ArgumentError("parameter coher must be between 0.0 and 1.0"))
		end
	end

	# constructing the singular spectrum
	
	levels = spectrum.levels
	ranges = spectrum.decayRanges
	types = spectrum.decayTypes
	
	S = Vector{Float64}(undef, dim)
	lastIndexSet = 0
	
	for typeIndex = 1:length(types)
		levelIndex = typeIndex
		rangeIndex = 2*(typeIndex - 1) + 1
		
		# filling out a flat portion of the spectrum before a decay region
		
		S[lastIndexSet + 1:ranges[rangeIndex] - 1] = levels[levelIndex]*ones(ranges[rangeIndex] - lastIndexSet - 1)
		
		# populating a decay region from ranges[rangeIndex] to ranges[rangeIndex + 1]
		
		numDecayPoints = ranges[rangeIndex + 1] - ranges[rangeIndex] + 1
		
		if(types[typeIndex] == "smoothgap")
			t = Array(range(-2, 2.5, numDecayPoints))
			decayCurve = .5*(erfc.(t) .+ 1)                               # set x to be an erfc decay curve
			beta = log(levels[levelIndex]/levels[levelIndex + 1])/log(3)  # correcting the height of the decay curve
			decayCurve = decayCurve.^beta
			decayCurve .*= levels[levelIndex]/decayCurve[1]               # making sure the curve starts at the right value
		else
			rho = (levels[levelIndex + 1]/levels[levelIndex])^(1/(numDecayPoints - 1))
			decayCurve = levels[levelIndex]*rho.^(0:numDecayPoints - 1)
		end
		
		S[ranges[rangeIndex]:ranges[rangeIndex + 1]] = decayCurve
		lastIndexSet = ranges[rangeIndex + 1]
	end
	
	# filling out a remaining flat portion of the spectrum
	
	S[lastIndexSet + 1:end] = levels[end]*ones(dim - lastIndexSet)

	if(minimal) return S end
	
	# constructing the right singular subspaces
	
	if(coher == nothing)
		Vt = Matrix(qr(randn(rng, dim, dim)).Q)
	else
		H = hadamard(dim)
		H = H[randperm(rng, dim), randperm(rng, dim)]/sqrt(dim)
		P = I(dim)
		P = P[:, randperm(rng, dim)]

		s = svd(coher*P + (1 - coher)*H)
		Vt = s.U*s.Vt
	end
	
	# constructing the left singular subspaces
	
	U = Matrix(qr(randn(rng, dim, dim)).Q)
	
	return U, S, Vt
end

function create_matrix(dim::Int, spectrum::SpectrumObject; coher::Union{Real, Nothing} = nothing, minimal::Bool = false)
	create_matrix(Random.default_rng(), dim, spectrum, coher = coher, minimal = minimal)
end
