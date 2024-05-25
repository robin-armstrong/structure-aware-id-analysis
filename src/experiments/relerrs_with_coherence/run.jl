using LinearAlgebra
using Distributions
using StatsBase
using Random
using PyPlot
using JLD2

include("../../algorithms/rgks.jl")
include("../../algorithms/rcpqr.jl")
include("../../algorithms/rsvd.jl")
include("../../algorithms/levg.jl")
include("../../utilities/create_matrix.jl")
include("../../utilities/plot_config.jl")

####################################################################
##################### SCRIPT PARAMETERS ############################
####################################################################

# path prefix for all output generated by this script
destination = "src/experiments/relerrs_with_coherence/alg_comparison"
readme      = "Testing low-rank approximation algorithms on two matrices with different subspace geometries."

rng        = MersenneTwister(1)
n          = 512
spectrum   = SpectrumObject([1., 1e-4], [1, 250], ["smoothgap"])
krange     = 1:250
coherences = [.2, .5]
numtrials  = 100

plot_only = false

####################################################################
##################### DATA GENERATION ##############################
####################################################################

if(!plot_only)
    function fprintln(s)
        println(s)
        flush(stdout)
    end

    function run_relerrs_with_coherence(destination, readme, rng, n, spectrum, krange, coherences, numtrials)
        logstr  = readme*"\n\n"
        logstr *= "rng        = "*string(rng)*"\n"
        logstr *= "n          = "*string(n)*"\n"
        logstr *= "spectrum   = "*string(spectrum)*"\n"
        logstr *= "krange     = "*string(krange)*"\n"
        logstr *= "coherences = "*string(coherences)*"\n"
        logstr *= "numtrials  = "*string(numtrials)*"\n"

        fprintln("\n"*logstr)
        
        logfile = destination*"_log.txt"
        touch(logfile)
        io = open(logfile, "w")
        write(logfile, logstr)
        close(io)

        A  = zeros(2, n, n)
        U  = zeros(2, n, n)
        Vt = zeros(2, n, n)
        S  = zeros(2, n)
        
        U1, S1, Vt1 = create_matrix(rng, n, spectrum, coher = coherences[1])
	    U[1, :, :]  = U1
        Vt[1, :, :] = Vt1
        S[1, :]     = S1
        A[1, :, :]  = U1*Diagonal(S1)*Vt1

        U2, S2, Vt2 = create_matrix(rng, n, spectrum, coher = coherences[2])
	    U[2, :, :]  = U2
        Vt[2, :, :] = Vt2
        S[2, :]     = S2
        A[2, :, :]  = U2*Diagonal(S2)*Vt2

        err_data   = Dict()
        means  = Dict()
        stds   = Dict()
        quants = Dict()
        
        for alg in ["rsvd_q0", "rcpqr", "rgks"]
            err_data[alg]   = zeros(2, length(krange), numtrials)
            means[alg]  = zeros(2, length(krange))
            stds[alg]   = zeros(2, length(krange))
            quants[alg] = zeros(2, length(krange), 2)
        end

        coher_data = zeros(2, length(krange))

        totalNumMatrices = 2*length(krange)*numtrials
        matrixCounter    = 0

        for matrix_id = 1:2
            for i = 1:length(krange)
                k           = krange[i]
                levg_scores = zeros(n)

                for j = 1:n
                    levg_scores[j] = norm(Vt[matrix_id, 1:k, j])^2
                end

                coher_data[matrix_id, i] = maximum(levg_scores)

                for t = 1:numtrials
                    matrixCounter += 1

                    if(matrixCounter % 10 == 0)
                        fprintln("   testing matrix "*string(matrixCounter)*" out of "*string(totalNumMatrices))
                    end

                    ov = ceil(Int64, .1*k)

                    U_appx, S_appx, Vt_appx = rsvd(rng, A[matrix_id, :, :], k, oversamp = ov, power = 0)
                    err_data["rsvd_q0"][matrix_id, i, t] = norm(A[matrix_id, :, :] - U_appx*Diagonal(S_appx)*Vt_appx)
                    
                    _, Q, X = rcpqr(rng, A[matrix_id, :, :], k, oversamp = ov)
                    err_data["rcpqr"][matrix_id, i, t] = norm(A[matrix_id, :, :] - Q*X)

                    _, Q, X = rgks(rng, A[matrix_id, :, :], k, oversamp = ov)
                    err_data["rgks"][matrix_id, i, t] = norm(A[matrix_id, :, :] - Q*X)
                end

                @save destination*"_data.jld2" krange numtrials U S Vt coher_data err_data means stds quants
            end
        end

        for alg in ["rsvd_q0", "rcpqr", "rgks"]
            for matrix_id = 1:2
                means[alg][matrix_id, :] = vec(mean(err_data[alg][matrix_id, :, :], dims = 2))
                stds[alg][matrix_id, :] = vec(std(err_data[alg][matrix_id, :, :], dims = 2))
                
                for i = 1:length(krange)
                    quants[alg][matrix_id, i, 1] = quantile(err_data[alg][matrix_id, i, :], .05)
                    quants[alg][matrix_id, i, 2] = quantile(err_data[alg][matrix_id, i, :], .95)
                end
            end
        end

        @save destination*"_data.jld2" krange numtrials U S Vt coher_data err_data means stds quants
    end

    run_relerrs_with_coherence(destination, readme, rng, n, spectrum, krange, coherences, numtrials)
end

####################################################################
##################### PLOTTING #####################################
####################################################################

@load destination*"_data.jld2" krange numtrials U S Vt coher_data err_data means stds quants

optimal     = zeros(2, length(krange))
matrix_norm = zeros(2)

for matrix_id = 1:2
    matrix_norm[matrix_id] = norm(S[matrix_id, :])
    optimal[matrix_id, :]  = [norm(S[matrix_id, (k + 1):end]) for k in krange]
end

fig, (err1, coher1, err2, coher2) = subplots(2, 2, height_ratios = [1, .5], figsize = (8, 6))
err   = [err1, err2]
coher = [coher1, coher2]

mfreq = 20
errbar     = "confidence"   # either "confidence" or "quantile"
confidence = .95

alpha = quantile(Normal(0, 1), 1 - .5*(1 - confidence))

for matrix_id = 1:2
    if(matrix_id == 1)
        err[matrix_id].set_ylabel("Relative Frobenius Error")
        coher[matrix_id].set_ylabel(L"Coherence ($c_k$)")
    end
    
    coher[matrix_id].set_xlabel(L"Approximation Rank ($k$)")
    coher[matrix_id].set_ylim([0, 1])
    coher[matrix_id].plot(krange, sqrt.(coher_data[matrix_id, :]), color = "black")
    
    err[matrix_id].set_yscale("log")
    err[matrix_id].plot(krange, optimal[matrix_id, :]/matrix_norm[matrix_id], color = "black", linestyle = "dashed", label = "Optimal")

    for alg in ["rsvd_q0", "rcpqr", "rgks"]
        (alg == "rsvd_q1") && continue
        (alg == "levg") && continue

        err[matrix_id].plot(krange, means[alg][matrix_id, :]/matrix_norm[matrix_id], color = algcolors[alg], marker = algmarkers[alg], markevery = mfreq, markerfacecolor = "none", label = alglabels[alg])

        if(errbar == "confidence")
            bar1 = (means[alg][matrix_id, :] .- alpha*stds[alg][matrix_id, :]/sqrt(numtrials))/matrix_norm[matrix_id]
            bar2 = (means[alg][matrix_id, :] .+ alpha*stds[alg][matrix_id, :]/sqrt(numtrials))/matrix_norm[matrix_id]
            err[matrix_id].fill_between(krange, bar1, bar2, color = algcolors[alg], alpha = .2)
        
        elseif(errbar == "quantile")
            err[matrix_id].fill_between(krange, quants[alg][matrix_id, :, 1]/matrix_norm[matrix_id], quants[alg][matrix_id, :, 2], color = algcolors[alg], alpha = .2)
        else
            throw(ArgumentError("unrecognized error bar type, '"*errbar*"'"))
        end
    end    
end

err1.legend()
savefig(destination*"_plot.pdf", bbox_inches = "tight")
close(fig)
