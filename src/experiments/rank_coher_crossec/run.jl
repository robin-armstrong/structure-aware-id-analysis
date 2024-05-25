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
include("../../utilities/create_matrix.jl")

####################################################################
##################### SCRIPT PARAMETERS ############################
####################################################################

# path prefix for all output generated by this script
destination = "src/experiments/rank_coher_crossec/crossections"
readme      = "Low-rank approximation algorithms compared over cross-sections of approximation rank and coherence."

rng          = MersenneTwister(1)   # random souce, set for reproducibility
n            = 4096
decay_start  = 1
decay_end    = 200
sigma_start  = 100.
sigma_end    = 1.
krange       = 1:5:300          # range of approximation ranks to test
num_cohers   = 100
k_cross      = 100
c_cross      = .2               # a number between 0 and 1
numtrials    = 100              # trials per approximation rank

plot_only = false

####################################################################
##################### DATA GENERATION ##############################
####################################################################

function fprintln(s)
    println(s)
    flush(stdout)
end

if(!plot_only)
    function run_crossecs(destination, readme, rng, n, decay_Start, decay_end, sigma_start, sigma_end, krange, num_cohers, k_cross, c_cross, numtrials)
        logstr  = readme*"\n\n"
        logstr *= "rng         = "*string(rng)*"\n"
        logstr *= "n           = "*string(n)*"\n"
        logstr *= "decay_start = "*string(decay_start)*"\n"
        logstr *= "decay_end   = "*string(decay_end)*"\n"
        logstr *= "sigma_start = "*string(sigma_start)*"\n"
        logstr *= "sigma_end   = "*string(sigma_end)*"\n"
        logstr *= "krange      = "*string(krange)*"\n"
        logstr *= "num_cohers  = "*string(num_cohers)*"\n"
        logstr *= "k_cross     = "*string(k_cross)*"\n"
        logstr *= "c_cross     = "*string(c_cross)*"\n"
        logstr *= "numtrials   = "*string(numtrials)*"\n"

        fprintln("\n"*logstr)
        
        logfile = destination*"_log.txt"
        touch(logfile)
        io = open(logfile, "w")
        write(logfile, logstr)
        close(io)

        data_vsrank    = Dict()
        means_vsrank   = Dict()
        stds_vsrank    = Dict()
        quants_vsrank  = Dict()
        data_vscoher   = Dict()
        means_vscoher  = Dict()
        stds_vscoher   = Dict()
        quants_vscoher = Dict()
        
        for alg in ["levg", "rsvd_q0", "rgks"]
            data_vsrank[alg]   = zeros(length(krange), numtrials)
            means_vsrank[alg]  = zeros(length(krange))
            stds_vsrank[alg]   = zeros(length(krange))
            quants_vsrank[alg] = zeros(length(krange), 2)
            
            data_vscoher[alg]   = zeros(num_cohers, numtrials)
            means_vscoher[alg]  = zeros(num_cohers)
            stds_vscoher[alg]   = zeros(num_cohers)
            quants_vscoher[alg] = zeros(num_cohers, 2)
        end
        
        cohers   = zeros(num_cohers)
        spct     = SpectrumObject([sigma_start, sigma_end], [decay_start, decay_end], ["smoothgap"])
        U, S, _  = create_matrix(rng, n, spct)

        totaltrials_rank  = numtrials*length(krange)

        H = hadamard(n)
		H = H[randperm(rng, n), randperm(rng, n)]/sqrt(n)
		P = I(n)
		P = P[:, randperm(rng, n)]

        svdobj = svd(c_cross*P + (1 - c_cross)*H)
		Vt     = svdobj.U*svdobj.Vt
        A      = U*Diagonal(S)*Vt

        fprintln("performing tests over approximation rank...\n")

        trialcounter = 0

        for i = 1:length(krange)
            k       = krange[i]
            ov      = ceil(Int64, .1*k)
            lscores = [norm(Vt[1:k, :])^2 for j = 1:n]
            lscores = min.(lscores, 1.)
            lscores = max.(lscores, 0.)

            for t = 1:numtrials
                trialcounter += 1

                if(trialcounter % 10 == 0)
                    fprintln("  testing matrix "*string(trialcounter)*" out of "*string(totaltrials_rank))
                end

                r1 = levg(rng, A, k, oversamp = ov, leverage_scores = lscores)
                data_vsrank["levg"][i, t] = norm(A - r1.Q*r1.X)

                r2 = rsvd(rng, A, k, oversamp = ov)
                data_vsrank["rsvd_q0"][i, t] = norm(A - r2.U*Diagonal(r2.S)*r2.Vt)

                r3 = rgks(rng, A, k, oversamp = ov)
                data_vsrank["rgks"][i, t] = norm(A - r3.Q*r3.X)
            end

            @save destination*"_data.jld2" krange cohers k_cross c_cross numtrials S data_vsrank means_vsrank stds_vsrank quants_vsrank data_vscoher means_vscoher stds_vscoher quants_vscoher
        end

        fprintln("\nperforming tests over coherence...\n")

        crange = range(0, 1, num_cohers)
        ov     = ceil(Int64, .1*k_cross)
        
        trialcounter = 0

        for i = 1:num_cohers
            c         = crange[i]
            svdobj    = svd(c*P + (1 - c)*H)
		    Vt        = svdobj.U*svdobj.Vt
            A         = U*Diagonal(S)*Vt
            lscores   = [norm(Vt[1:k_cross, j])^2 for j = 1:n]
            lscores   = min.(lscores, 1.)
            lscores   = max.(lscores, 0.)
            cohers[i] = maximum(sqrt.(lscores))

            for t = 1:numtrials
                trialcounter += 1

                if(trialcounter % 10 == 0)
                    fprintln("  testing matrix "*string(trialcounter)*" out of "*string(totaltrials_rank))
                end

                r1 = levg(rng, A, k_cross, oversamp = ov, leverage_scores = lscores)
                data_vscoher["levg"][i, t] = norm(A - r1.Q*r1.X)

                r2 = rsvd(rng, A, k_cross, oversamp = ov)
                data_vscoher["rsvd_q0"][i, t] = norm(A - r2.U*Diagonal(r2.S)*r2.Vt)

                r3 = rgks(rng, A, k_cross, oversamp = ov)
                data_vscoher["rgks"][i, t] = norm(A - r3.Q*r3.X)
            end

            @save destination*"_data.jld2" krange cohers k_cross c_cross numtrials S data_vsrank means_vsrank stds_vsrank quants_vsrank data_vscoher means_vscoher stds_vscoher quants_vscoher
        end

        println("\nsorting data based on coherence...")

        sp     = sortperm(cohers)
        cohers = cohers[sp]
        
        for alg in ["levg", "rsvd_q0", "rgks"]
            data_vscoher[alg] = data_vscoher[alg][sp, :]
        end

        fprintln("calculating approximation error statistics...")

        for alg in ["levg", "rsvd_q0", "rgks"]
            means_vsrank[alg]  = vec(mean(data_vsrank[alg], dims = 2))
            stds_vsrank[alg]   = vec(std(data_vsrank[alg], dims = 2))

            for i = 1:length(krange)
                quants_vsrank[alg][i, 1] = quantile(data_vsrank[alg][i, :], .05)
                quants_vsrank[alg][i, 2] = quantile(data_vsrank[alg][i, :], .95)
            end

            means_vscoher[alg]  = vec(mean(data_vscoher[alg], dims = 2))
            stds_vscoher[alg]   = vec(std(data_vscoher[alg], dims = 2))

            for i = 1:num_cohers
                quants_vscoher[alg][i, 1] = quantile(data_vscoher[alg][i, :], .05)
                quants_vscoher[alg][i, 2] = quantile(data_vscoher[alg][i, :], .95)
            end
        end

        @save destination*"_data.jld2" krange cohers k_cross c_cross numtrials S data_vsrank means_vsrank stds_vsrank quants_vsrank data_vscoher means_vscoher stds_vscoher quants_vscoher
    end

    run_crossecs(destination, readme, rng, n, decay_start, decay_end, sigma_start, sigma_end, krange, num_cohers, k_cross, c_cross, numtrials)
end

####################################################################
##################### PLOTTING #####################################
####################################################################

@load destination*"_data.jld2" krange cohers k_cross c_cross numtrials S data_vsrank means_vsrank stds_vsrank quants_vsrank data_vscoher means_vscoher stds_vscoher quants_vscoher

fprintln("plotting error statistics...")

matrixnorm    = norm(S)
optimal       = [norm(S[(k + 1):end]) for k in krange]
optimal_cross = norm(S[(k_cross + 1):end])

mfreq      = 5
errbar     = "quantile"   # either "confidence" or "quantile"
confidence = .95

alpha = quantile(Normal(0, 1), 1 - .5*(1 - confidence))

ioff()
fig, (vsrank_rel, vscoher, vsrank_opt, stablerank) = subplots(2, 2, figsize = (8, 8))

vsrank_rel.set_xlabel(L"Approximation Rank ($k$)")
vsrank_rel.set_ylabel("Relative Frobenius Error")
vsrank_rel.set_yscale("log")
stablerank.set_xlabel(L"Approximation Rank ($k$)")
stablerank.set_ylabel(L"Log Stable Rank ($\log_{10}(r_k)$)")
stablerank.set_ylim([1, 3.7])
vsrank_opt.set_xlabel(L"Approximation Rank ($k$)")
vsrank_opt.set_ylabel("Frobenius Error Suboptimality")
vscoher.set_xlabel(L"Coherence ($c_{100}$)")
vscoher.set_ylabel(L"Frobenius Suboptimality ($k = 100$)")

vsrank_opt.set_yticks([1, 2])
stablerank.set_yticks([1, 2, 3])

for alg in ["levg", "rsvd_q0", "rgks"]
    vsrank_rel.plot(krange, means_vsrank[alg]/matrixnorm, color = algcolors[alg], marker = algmarkers[alg], markevery = mfreq, markerfacecolor = "none", label = alglabels[alg])
    vsrank_opt.plot(krange, means_vsrank[alg]./optimal, color = algcolors[alg], marker = algmarkers[alg], markevery = mfreq, markerfacecolor = "none", label = alglabels[alg])
    vscoher.plot(cohers, means_vscoher[alg]/optimal_cross, color = algcolors[alg], marker = algmarkers[alg], markevery = mfreq, markerfacecolor = "none")

    if(errbar == "quantile")
        vsrank_rel.fill_between(krange, quants_vsrank[alg][:, 1]/matrixnorm, quants_vsrank[alg][:, 2]/matrixnorm, color = algcolors[alg], alpha = .2)
        vsrank_opt.fill_between(krange, quants_vsrank[alg][:, 1]./optimal, quants_vsrank[alg][:, 2]./optimal, color = algcolors[alg], alpha = .2)
        vscoher.fill_between(cohers, quants_vscoher[alg][:, 1]/optimal_cross, quants_vscoher[alg][:, 2]/optimal_cross, color = algcolors[alg], alpha = .2)
    elseif(errbar == "confidence")
        vsrank_rel.fill_between(krange, (means_vsrank[alg] .+ alpha*stds_vsrank[alg]/sqrt(numtrials))/matrixnorm, (means_vsrank[alg] .- alpha*stds_vsrank[alg]/sqrt(numtrials))/matrixnorm, color = algcolors[alg], alpha = .2)
        vsrank_opt.fill_between(krange, (means_vsrank[alg] .+ alpha*stds_vsrank[alg]/sqrt(numtrials))./optimal, (means_vsrank[alg] .- alpha*stds_vsrank[alg]/sqrt(numtrials))./optimal, color = algcolors[alg], alpha = .2)
        vscoher.fill_between(cohers, (means_vscoher[alg] .+ alpha*stds_vscoher[alg]/sqrt(numtrials))/optimal_cross, (means_vscoher[alg] .- alpha*stds_vscoher[alg]/sqrt(numtrials))/optimal_cross, color = algcolors[alg], alpha = .2)
    else
        throw(ArgumentError("unrecognized error bar type, '"*errbar*"'"))
    end
end

vsrank_rel.plot(krange, optimal/matrixnorm, color = "black", linestyle = "dashed", label = "Optimal")
vsrank_rel.legend()

sr = [norm(S[(k + 1):end])^2/S[k + 1]^2 for k in krange]
stablerank.plot(krange, log10.(sr), color = "black")

savefig(destination*"_plot.pdf", bbox_inches = "tight")
close(fig)