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
include("../../utilities/plot_config.jl")

####################################################################
##################### SCRIPT PARAMETERS ############################
####################################################################

# path prefix for all output generated by this script
destination = "src/experiments/kahan_example/oneblock"
readme      = "Testing low-rank approximation algorithms on a Kahan matrix."

rng       = MersenneTwister(1)      # random souce, set for reproducibility
n         = 100                     # dimension of Kahan blocks
blocks    = 1                       # number of Kahan blocks
alpha     = .9                      # parameterizes the Kahan blocks
pert      = 1e3                     # controls a diagonal perturbation on the test matrix
krange    = 1:2:100                 # range of approximation ranks to test
numtrials = 100                     # trials per approximation rank

plot_only = false

####################################################################
##################### DATA GENERATION ##############################
####################################################################

if(!plot_only)
    function fprintln(s)
        println(s)
        flush(stdout)
    end

    function kahan(n, alpha, pert)
        beta = sqrt(max(1. - alpha^2, 0.))
        K    = triu(-beta*ones(n, n))
        d    = 1
        
        for i = 1:n
            K[i, i]  = 1. + pert*eps()*(n - i + 1)
            K[i, :] *= d
            d       *= alpha
        end

        return K
    end
    
    function run_kahan_example(destination, readme, rng, n, blocks, alpha, pert, krange, numtrials, plot_only)
        logstr  = readme*"\n\n"
        logstr *= "rng       = "*string(rng)*"\n"
        logstr *= "n         = "*string(n)*"\n"
        logstr *= "blocks    = "*string(blocks)*"\n"
        logstr *= "alpha     = "*string(alpha)*"\n"
        logstr *= "pert      = "*string(pert)*"\n"
        logstr *= "krange    = "*string(krange)*"\n"
        logstr *= "numtrials = "*string(numtrials)*"\n"

        fprintln("\n"*logstr)
        
        logfile = destination*"_log.txt"
        touch(logfile)
        io = open(logfile, "w")
        write(logfile, logstr)
        close(io)

        fprintln("generating block Kahan matrix...")

        K  = zeros(n*blocks, n*blocks)
        
        for b = 1:blocks
            id1 = (b - 1)*n + 1
            id2 = b*n

            K[id1:id2, id1:id2] = kahan(n, alpha, pert)
        end

        fprintln("computing SVD of test matrix...")

        block_kahan_svd = svd(K)

        data   = Dict()
        means  = Dict()
        stds   = Dict()
        quants = Dict()

        for alg in ["dgeqp3", "levg", "rcpqr", "rgks"]
            data[alg]   = zeros(length(krange), numtrials)
            means[alg]  = zeros(length(krange))
            stds[alg]   = zeros(length(krange))
            quants[alg] = zeros(length(krange), 2)
        end

        fprintln("running tests...\n")

        trialcounter = 0

        for t = 1:numtrials
            for i = 1:length(krange)
                trialcounter += 1
                
                k       = krange[i]
                lscores = [norm(block_kahan_svd.Vt[1:k, j])^2 for j = 1:size(K, 2)]
                lscores = min.(lscores, 1.)
                lscores = max.(lscores, 0.)

                if(trialcounter % 10 == 0)
                    fprintln("   trial "*string(trialcounter)*" out of "*string(numtrials*length(krange)))
                end

                r1 = rgks(rng, K, k, oversamp = ceil(Int64, .1*k))
                r2 = rcpqr(rng, K, k, oversamp = ceil(Int64, .1*k))
                r3 = levg(rng, K, k, oversamp = ceil(Int64, .1*k), leverage_scores = lscores)
                r4 = qr(K, ColumnNorm())

                # optimally reducing the levg approximation to rank k

                U = svd(r3.X).U
                Q = r3.Q*U[:, 1:k]
                X = Q'*K

                # finding the approximation subspace used by DGEQP3
                W = r4.Q*Matrix{Float64}(I(n)[:, 1:k])

                data["rgks"][i, t]    = norm(K - (r1.Q)*(r1.X))
                data["rcpqr"][i, t]     = norm(K - (r2.Q)*(r2.X))
                data["levg"][i, t]    = norm(K - Q*X)
                data["dgeqp3"][i, t]  = norm(K - W*(W'*K))
            end

            @save destination*"_data.jld2" krange numtrials block_kahan_svd data means stds quants
        end

        fprintln("\ncalculating approximation error statistics...")
        
        for alg in ["dgeqp3", "levg", "rcpqr", "rgks"]
            means[alg]  = vec(mean(data[alg], dims = 2))
            stds[alg]   = vec(std(data[alg], dims = 2))

            for i = 1:length(krange)
                quants[alg][i, 1] = quantile(data[alg][i, :], .05)
                quants[alg][i, 2] = quantile(data[alg][i, :], .95)
            end
        end

        @save destination*"_data.jld2" krange numtrials block_kahan_svd data means stds quants
    end

    run_kahan_example(destination, readme, rng, n, blocks, alpha, pert, krange, numtrials, plot_only)
end

####################################################################
##################### PLOTTING #####################################
####################################################################

@load destination*"_data.jld2" krange numtrials block_kahan_svd data means stds quants

matrix_norm = norm(block_kahan_svd.S)
optimal     = [norm(block_kahan_svd.S[(k + 1):end]) for k in krange]

ioff()
fig, (norm_rel, opt_rel) = subplots(1, 2, figsize = (10, 4))

mfreq      = 5
errbar     = "confidence"   # either "confidence" or "quantile"
confidence = .95

alpha = quantile(Normal(0, 1), 1 - .5*(1 - confidence))

norm_rel.set_xlabel(L"Approximation Rank ($k$)")
norm_rel.set_ylabel("Relative Frobenius Error")
norm_rel.set_yscale("log")

opt_rel.set_xlabel(L"Approximation Rank ($k$)")
opt_rel.set_ylabel("Frobenius Error Suboptimality")
opt_rel.set_ylim([1., 3.])

for alg in ["dgeqp3", "levg", "rcpqr", "rgks"]
    norm_rel.plot(krange, means[alg]/matrix_norm, color = algcolors[alg], marker = algmarkers[alg], markevery = mfreq, markerfacecolor = "none", label = alglabels[alg])
    opt_rel.plot(krange, means[alg]./optimal, color = algcolors[alg], marker = algmarkers[alg], markevery = mfreq, markerfacecolor = "none")

    if(errbar == "quantile")
        norm_rel.fill_between(krange, quants[alg][:, 1]/matrix_norm, quants[alg][:, 2]/matrix_norm, color = algcolors[alg], alpha = .2)
        opt_rel.fill_between(krange, quants[alg][:, 1]./optimal, quants[alg][:, 2]./optimal, color = algcolors[alg], alpha = .2)
    elseif(errbar == "confidence")
        norm_rel.fill_between(krange, (means[alg] .+ alpha*stds[alg]/sqrt(numtrials))/matrix_norm, (means[alg] .- alpha*stds[alg]/sqrt(numtrials))/matrix_norm, color = algcolors[alg], alpha = .2)
        opt_rel.fill_between(krange, (means[alg] .+ alpha*stds[alg]/sqrt(numtrials))./optimal, (means[alg] .- alpha*stds[alg]/sqrt(numtrials))./optimal, color = algcolors[alg], alpha = .2)
    else
        throw(ArgumentError("unrecognized error bar type, '"*errbar*"'"))
    end
end

norm_rel.plot(krange, optimal/matrix_norm, color = "black", linestyle = "dashed", label = "Optimal")
norm_rel.legend()

savefig(destination*"_plot.pdf", bbox_inches = "tight")
close(fig)
