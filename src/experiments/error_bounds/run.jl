using LinearAlgebra
using Distributions
using StatsBase
using Random
using PyPlot
using JLD2

include("../../utilities/plot_config.jl")
include("../../utilities/create_matrix.jl")

####################################################################
##################### SCRIPT PARAMETERS ############################
####################################################################

# DEFAULT PARAMETER SETTINGS:
#
# rng         = MersenneTwister(1)
# n           = 512
# decay_start = 20
# decay_end   = 80
# sigma_start = 10.
# sigma_end   = .001
# krange      = 1:1:100
# num_cohers  = 50
# k_cross     = 45
# c_cross     = .2

destination = "src/experiments/error_bounds/errorbounds"
readme      = "Trying to get the stupid script to work."

rng         = MersenneTwister(1)
n           = 512
decay_start = 20
decay_end   = 80
sigma_start = 10.
sigma_end   = .001
krange      = 1:1:100
num_cohers  = 50
k_cross     = 45
c_cross     = .2

plot_only = false

####################################################################
##################### DATA GENERATION ##############################
####################################################################

function fprintln(s)
    println(s)
    flush(stdout)
end

if(!plot_only)
    function run_errbound_plots(rng, n, decay_start, decay_end, sigma_start, sigma_end, krange, num_cohers, k_cross, c_cross)
        logstr  = readme*"\n\n"
        logstr *= "rng         = "*string(rng)*"\n"
        logstr *= "n           = "*string(n)*"\n"
        logstr *= "decay_start = "*string(decay_start)*"\n"
        logstr *= "decay_end   = "*string(decay_end)*"\n"
        logstr *= "sigma_start = "*string(sigma_start)*"\n"
        logstr *= "sigma_end   = "*string(sigma_end)*"\n"

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
        
        for key in ["spec_err", "frob_err", "spec_bound", "frob_bound", "cond_bound"]
            data_vsrank[key]   = zeros(length(krange))
            data_vscoher[key]   = zeros(num_cohers)
        end
        
        cohers   = zeros(num_cohers)
        spct     = SpectrumObject([sigma_start, sigma_end], [decay_start, decay_end], ["smoothgap"])
        U, S, _  = create_matrix(rng, n, spct)

        H = hadamard(n)
		H = H[randperm(rng, n), randperm(rng, n)]/sqrt(n)
		P = I(n)
		P = P[:, randperm(rng, n)]

        svdobj = svd(c_cross*P + (1 - c_cross)*H)
		Vt     = svdobj.U*svdobj.Vt
        A      = U*Diagonal(S)*Vt

        fprintln("performing tests over approximation rank...\n")

        trialcounter    = 0
        cos_to_sqtan(x) = x.^(-2) .- 1

        for i = 1:length(krange)
            k = krange[i]

            if(trialcounter % 10 == 0)
                fprintln("  testing matrix "*string(i)*" out of "*string(length(krange)))
            end

            J     = qr(Vt[1:k, :], ColumnNorm()).p[1:k]
            Q     = Matrix{Float64}(qr(A[:,J]).Q)
            E     = A - Q*(Q'*A)
            E_svd = svd(E)
            V_svd = svd(Vt[1:k, J])

            sr = norm(S[(k + 1):end])^2/S[k + 1]^2
            
            data_vsrank["spec_err"][i]   = E_svd.S[1]
            data_vsrank["frob_err"][i]   = norm(E_svd.S)
            data_vsrank["spec_bound"][i] = S[k + 1]/V_svd.S[k]
            data_vsrank["frob_bound"][i] = norm(S[(k + 1):end])*sqrt(1 + sum(cos_to_sqtan.(V_svd.S))/sr)
            data_vsrank["cond_bound"][i] = S[k + 1]*E_svd.S[1]/E_svd.S[k + 1]
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
            cohers[i] = maximum(sqrt.(lscores))


            if(trialcounter % 10 == 0)
                fprintln("  testing matrix "*string(i)*" out of "*string(num_cohers))
            end

            J     = qr(Vt[1:k_cross, :], ColumnNorm()).p[1:k_cross]
            Q     = Matrix{Float64}(qr(A[:,J]).Q)
            E     = A - Q*(Q'*A)
            E_svd = svd(E)
            V_svd = svd(Vt[1:k_cross, J])

            sr = norm(S[(k_cross + 1):end])^2/S[k_cross + 1]^2
            
            data_vscoher["spec_err"][i]   = E_svd.S[1]
            data_vscoher["frob_err"][i]   = norm(E_svd.S)
            data_vscoher["spec_bound"][i] = S[k_cross + 1]/V_svd.S[k_cross]
            data_vscoher["frob_bound"][i] = norm(S[(k_cross + 1):end])*sqrt(1 + sum(cos_to_sqtan.(V_svd.S))/sr)
            data_vscoher["cond_bound"][i] = S[k_cross + 1]*E_svd.S[1]/E_svd.S[k_cross + 1]
        end

        println("\nsorting data based on coherence...")

        sp     = sortperm(cohers)
        cohers = cohers[sp]
        
        for key in ["spec_err", "frob_err", "spec_bound", "frob_bound", "cond_bound"]
            data_vscoher[key] = data_vscoher[key][sp]
        end

        @save destination*"_data.jld2" krange cohers k_cross c_cross S data_vsrank data_vscoher
    end

    run_errbound_plots(rng, n, decay_start, decay_end, sigma_start, sigma_end, krange, num_cohers, k_cross, c_cross)
end

####################################################################
##################### PLOTTING #####################################
####################################################################

@load destination*"_data.jld2" krange cohers k_cross c_cross S data_vsrank data_vscoher

spec_optimal = [S[k + 1] for k in krange]
frob_optimal = [norm(S[(k + 1):end]) for k in krange]
gamma        = [S[k + 1]/S[k] for k in krange]

ioff()
fig, (spec_rank, frob_rank, cond_rank, spec_coher, frob_coher, cond_coher) = subplots(3, 2, figsize = (8, 12))

for plot in [spec_rank, cond_rank]
    plot.set_ylabel("Spectral Error Suboptimality")
end

frob_rank.set_ylabel("Frobenius Error Suboptimality")

for plot in [spec_rank, frob_rank, cond_rank]
    plot.set_xlabel(L"Approximation Rank ($k$)")
end

for plot in [spec_coher, frob_coher, cond_coher]
    plot.set_xlabel(L"Coherence ($c_{45}$)")
end

mfreq = 10

spec_rank.set_yscale("log")
spec_rank.plot(krange, data_vsrank["spec_err"]./spec_optimal, color = "black", label = "GKS")
spec_rank.plot(krange, data_vsrank["spec_bound"]./spec_optimal, color = "brown", marker = "D", markerfacecolor = "none", markevery = mfreq, label = "Thm 5.1")
spec_rank.plot(krange, gamma, color = "black", linestyle = "dotted", label = L"\sigma_{k + 1}/\sigma_k")
spec_rank.legend()

frob_rank.set_yscale("log")
frob_rank.plot(krange, data_vsrank["frob_err"]./frob_optimal, color = "black", label = "GKS")
frob_rank.plot(krange, data_vsrank["frob_bound"]./frob_optimal, color = "brown", marker = "D", markerfacecolor = "none", markevery = mfreq, label = "Thm 5.3")
frob_rank.plot(krange, gamma, color = "black", linestyle = "dotted", label = L"\sigma_{k + 1}/\sigma_k")
frob_rank.legend()

cond_rank.set_yscale("log")
cond_rank.plot(krange, data_vsrank["spec_err"]./spec_optimal, color = "black", label = "GKS")
cond_rank.plot(krange, data_vsrank["cond_bound"]./spec_optimal, color = "brown", marker = "D", markerfacecolor = "none", markevery = mfreq, label = "Thm 5.4")
cond_rank.plot(krange, gamma, color = "black", linestyle = "dotted", label = L"\sigma_{k + 1}/\sigma_k")
cond_rank.legend()

spec_coher.set_yscale("log")
spec_coher.plot(cohers, data_vscoher["spec_err"]/S[k_cross + 1], color = "black", label = "GKS")
spec_coher.plot(cohers, data_vscoher["spec_bound"]/S[k_cross + 1], color = "brown", marker = "D", markerfacecolor = "none", markevery = mfreq, label = "Thm 5.1")
spec_coher.legend()

frob_coher.set_yscale("log")
frob_coher.plot(cohers, data_vscoher["frob_err"]/norm(S[(k_cross + 1):end]), color = "black", label = "GKS")
frob_coher.plot(cohers, data_vscoher["frob_bound"]/norm(S[(k_cross + 1):end]), color = "brown", marker = "D", markerfacecolor = "none", markevery = mfreq, label = "Thm 5.3")
frob_coher.legend()

cond_coher.set_yscale("log")
cond_coher.plot(cohers, data_vscoher["spec_err"]/S[k_cross + 1], color = "black", label = "GKS")
cond_coher.plot(cohers, data_vscoher["cond_bound"]/S[k_cross + 1], color = "brown", marker = "D", markerfacecolor = "none", markevery = mfreq, label = "Thm 5.4")
cond_coher.legend()

savefig(destination*"_plot.pdf", bbox_inches = "tight")
close(fig)
