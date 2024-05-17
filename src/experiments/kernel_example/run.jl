using LinearAlgebra
using Distributions
using StatsBase
using Random
using PyPlot
import Plots
using JLD2

include("../../algorithms/rgks.jl")
include("../../algorithms/rid.jl")
include("../../algorithms/rsvd.jl")
include("../../algorithms/levg.jl")
include("../../utilities/plot_config.jl")

####################################################################
##################### SCRIPT PARAMETERS ############################
####################################################################

# path prefix for all output generated by this script
destination = "src/experiments/kernel_example/onecomponent"
readme      = "Testing low-rank approximation algorithms on a kernel evaluation matrix, generated from a Gaussian cloud."

rng          = MersenneTwister(1)   # random souce, set for reproducibility
num_clusters = 1                    # number of clusters in mixture model
num_points   = 2000                 # number of points per cluster
radius       = 0                    # clusters are equispaced around a circle of this radius
noise        = 5                    # standard deviation of Gaussian mixture components
bandwidth    = 2                    # bandwidth of Gaussian kernel
normalize    = false
krange       = 1:2:99               # range of approximation ranks to test
k_anim       = 80
numtrials    = 50                   # trials per approximation rank

plot_only = true

####################################################################
##################### DATA GENERATION ##############################
####################################################################

function fprintln(s)
    println(s)
    flush(stdout)
end

if(!plot_only)
    function run_kernel_example(destination, readme, rng, num_clusters, num_points, radius, noise, bandwidth, krange, numtrials, plot_only)
        logstr  = readme*"\n\n"
        logstr *= "rng          = "*string(rng)*"\n"
        logstr *= "num_clusters = "*string(num_clusters)*"\n"
        logstr *= "num_points   = "*string(num_points)*"\n"
        logstr *= "radius       = "*string(radius)*"\n"
        logstr *= "noise        = "*string(noise)*"\n"
        logstr *= "bandwidth    = "*string(bandwidth)*"\n"
        logstr *= "normalize    = "*string(normalize)*"\n"
        logstr *= "krange       = "*string(krange)*"\n"
        logstr *= "numtrials    = "*string(numtrials)*"\n"

        fprintln("\n"*logstr)
        
        logfile = destination*"_log.txt"
        touch(logfile)
        io = open(logfile, "w")
        write(logfile, logstr)
        close(io)

        angles        = range(0, 2*pi, num_clusters + 1)
        centers       = zeros(2, num_clusters)
        centers[1, :] = radius*cos.(angles[1:num_clusters])
        centers[2, :] = radius*sin.(angles[1:num_clusters])

        N       = num_clusters*num_points
        points  = zeros(2, N)
        sqnorms = zeros(N)
        
        fprintln("generating data points...")

        for i = 1:N
            c            = 1 + floor(Int64, num_clusters*rand(rng))
            points[:, i] = centers[:, c] .+ noise*randn(rng, 2)
            sqnorms[i]   = norm(points[:, i]).^2
        end

        fprintln("evaluating kernel matrix...")

        K  = zeros(N, N)
        K += sqnorms*ones(1, N)
        K += ones(N)*sqnorms'
        K -= 2*points'*points
        K /= -2*bandwidth^2
        broadcast!(exp, K, K)

        if(normalize)
            D = inv(Diagonal(sqrt.(K*ones(N))))
            lmul!(D, K)
            rmul!(K, D)
        end

        fprintln("computing SVD of kernel matrix...")

        kernel_svd = svd(K)

        data   = Dict()
        means  = Dict()
        stds   = Dict()
        quants = Dict()

        for alg in ["dgeqp3", "levg", "rid", "rgks"]
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
                lscores = [norm(kernel_svd.Vt[1:k, j])^2 for j = 1:size(K, 2)]
                lscores = min.(lscores, 1.)
                lscores = max.(lscores, 0.)

                if(trialcounter % 10 == 0)
                    fprintln("   trial "*string(trialcounter)*" out of "*string(numtrials*length(krange)))
                end

                r1 = rgks(rng, K, k, oversamp = ceil(Int64, .1*k))
                r2 = rid(rng, K, k, oversamp = ceil(Int64, .1*k))
                r3 = levg(rng, K, k, oversamp = ceil(Int64, .1*k), leverage_scores = lscores)
                r4 = qr(K, ColumnNorm())

                # optimally reducing the levg approximation to rank k

                U = svd(r3.X).U
                Q = r3.Q*U[:, 1:k]
                X = Q'*K

                # finding the approximation subspace used by DGEQP3
                W = r4.Q*Matrix{Float64}(I(N)[:, 1:k])

                data["rgks"][i, t]    = norm(K - (r1.Q)*(r1.X))
                data["rid"][i, t]     = norm(K - (r2.Q)*(r2.X))
                data["levg"][i, t]    = norm(K - Q*X)
                data["dgeqp3"][i, t]  = norm(K - W*(W'*K))
            end

            @save destination*"_data.jld2" krange numtrials points kernel_svd data means stds quants
        end

        fprintln("\ncalculating approximation error statistics...")
        
        for alg in ["dgeqp3", "levg", "rid", "rgks"]
            means[alg]  = vec(mean(data[alg], dims = 2))
            stds[alg]   = vec(std(data[alg], dims = 2))

            for i = 1:length(krange)
                quants[alg][i, 1] = quantile(data[alg][i, :], .05)
                quants[alg][i, 2] = quantile(data[alg][i, :], .95)
            end
        end

        @save destination*"_data.jld2" krange numtrials points kernel_svd data means stds quants
    end

    run_kernel_example(destination, readme, rng, num_clusters, num_points, radius, noise, bandwidth, krange, numtrials, plot_only)
end

####################################################################
##################### PLOTTING #####################################
####################################################################

@load destination*"_data.jld2" krange numtrials points kernel_svd data means stds quants

fprintln("plotting error statistics...")

kernel_norm = norm(kernel_svd.S)
optimal     = [norm(kernel_svd.S[(k + 1):end]) for k in krange]

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

for alg in ["dgeqp3", "levg", "rid", "rgks"]
    norm_rel.plot(krange, means[alg]/kernel_norm, color = algcolors[alg], marker = algmarkers[alg], markevery = mfreq, markerfacecolor = "none", label = alglabels[alg])
    opt_rel.plot(krange, means[alg]./optimal, color = algcolors[alg], marker = algmarkers[alg], markevery = mfreq, markerfacecolor = "none")

    if(errbar == "quantile")
        norm_rel.fill_between(krange, quants[alg][:, 1]/kernel_norm, quants[alg][:, 2]/kernel_norm, color = algcolors[alg], alpha = .2)
        opt_rel.fill_between(krange, quants[alg][:, 1]./optimal, quants[alg][:, 2]./optimal, color = algcolors[alg], alpha = .2)
    elseif(errbar == "confidence")
        norm_rel.fill_between(krange, (means[alg] .+ alpha*stds[alg]/sqrt(numtrials))/kernel_norm, (means[alg] .- alpha*stds[alg]/sqrt(numtrials))/kernel_norm, color = algcolors[alg], alpha = .2)
        opt_rel.fill_between(krange, (means[alg] .+ alpha*stds[alg]/sqrt(numtrials))./optimal, (means[alg] .- alpha*stds[alg]/sqrt(numtrials))./optimal, color = algcolors[alg], alpha = .2)
    else
        throw(ArgumentError("unrecognized error bar type, '"*errbar*"'"))
    end
end

norm_rel.plot(krange, optimal/kernel_norm, color = "black", linestyle = "dashed", label = "Optimal")
norm_rel.legend()

savefig(destination*"_plot.pdf", bbox_inches = "tight")
close(fig)

fprintln("creating point selection animation...")

K      = kernel_svd.U*Diagonal(kernel_svd.S)*kernel_svd.Vt
p_rid  = rid(K, k_anim, oversamp = ceil(Int64, .1*k_anim)).p
p_rgks = rgks(K, k_anim, oversamp = ceil(Int64, .1*k_anim)).p

colnorms = [norm(K[:, j]) for j = 1:size(K, 2)]
l_scores = [norm(kernel_svd.Vt[1:k_anim, j]) for j = 1:size(K, 2)]

rid_cloud  = Plots.plot()
rid_cols   = Plots.plot()
rgks_cloud = Plots.plot()
rgks_cols  = Plots.plot()

Plots.title!(rid_cloud, "RID")
Plots.title!(rgks_cloud, "RGKS")
Plots.xlabel!(rid_cols, "Leverage Score")
Plots.xlabel!(rgks_cols, "Leverage Score")
Plots.xlims!(rgks_cloud, minimum(points[1, :]), maximum(points[1, :]))
Plots.xlims!(rid_cloud, minimum(points[1, :]), maximum(points[1, :]))
Plots.xlims!(rid_cols, 0., 1.1*maximum(l_scores))
Plots.xlims!(rgks_cols, 0., 1.1*maximum(l_scores))
Plots.ylabel!(rid_cols, "Column Norm")
Plots.ylims!(rgks_cloud, minimum(points[2, :]), maximum(points[2, :]))
Plots.ylims!(rid_cloud, minimum(points[2, :]), maximum(points[2, :]))
Plots.ylims!(rid_cols, 0., 1.1*maximum(colnorms))
Plots.ylims!(rgks_cols, 0., 1.1*maximum(colnorms))

tiles = Plots.plot(rid_cloud, rgks_cloud, rid_cols, rgks_cols, layout = (2, 2))
Plots.scatter!(tiles[1], points[1, :], points[2, :], markercolor = :gray, markersize = 2, markeralpha = .1, legend = false)
Plots.scatter!(tiles[2], points[1, :], points[2, :], markercolor = :gray, markersize = 2, markeralpha = .1, legend = false)
Plots.scatter!(tiles[3], l_scores, colnorms, markercolor = :gray, markersize = 2, markeralpha = .1, legend = false)
Plots.scatter!(tiles[4], l_scores, colnorms, markercolor = :gray, markersize = 2, markeralpha = .1, legend = false)

anim = Plots.@animate for j = 1:k_anim
    fprintln("  frame "*string(j)*" of "*string(k_anim)*"...")
    Plots.scatter!(tiles[1], points[1, p_rid[j]], points[1, p_rid[j]], markercolor = :green, markersize = 4, legend = false)
    Plots.scatter!(tiles[2], points[1, p_rgks[j]], points[1, p_rgks[j]], markercolor = :red, markersize = 4, legend = false)
    Plots.scatter!(tiles[3], l_scores[p_rid[j]], colnorms[p_rid[j]], markercolor = :green, markersize = 4, legend = false)
    Plots.scatter!(tiles[4], l_scores[p_rgks[j]], colnorms[p_rgks[j]], markercolor = :red, markersize = 4, legend = false)
end

Plots.fig(anim, destination*"_animation.gif", fps = 3)