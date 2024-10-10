using LinearAlgebra
using Distributions
using StatsBase
using Random
using PyPlot
using JLD2

include("../../algorithms/rsvd.jl")
include("../../utilities/plot_config.jl")
include("../../utilities/create_matrix.jl")

####################################################################
##################### SCRIPT PARAMETERS ############################
####################################################################

# DEFAULT PARAMETER SETTINGS:
#
# rng          = MersenneTwister(1)
# n            = 512
# decay_start  = 10
# decay_end    = 90
# sigma_start  = 100.
# sigma_end    = .00001
# krange       = 10:1:90
# coherence    = .2
# numtrials    = 100

# path prefix for all output generated by this script
destination = "src/experiments/angle_perturbations/angles"
readme      = "Principal angles relevant to the error analysis of RGKS."

rng          = MersenneTwister(1)   # random souce, set for reproducibility
n            = 512
decay_start  = 10
decay_end    = 90
sigma_start  = 100.
sigma_end    = .00001
krange       = 10:1:90             # range of approximation ranks to test
coherence    = .2
numtrials    = 100

plot_only = false

####################################################################
##################### DATA GENERATION ##############################
####################################################################

function fprintln(s)
    println(s)
    flush(stdout)
end

if(!plot_only)
    function run_angle_perturbations(destination, readme, rng, n, decay_start, decay_end, sigma_start, sigma_end, krange, coherence, numtrials)
        logstr  = readme*"\n\n"
        logstr *= "rng         = "*string(rng)*"\n"
        logstr *= "n           = "*string(n)*"\n"
        logstr *= "decay_start = "*string(decay_start)*"\n"
        logstr *= "decay_end   = "*string(decay_end)*"\n"
        logstr *= "sigma_start = "*string(sigma_start)*"\n"
        logstr *= "sigma_end   = "*string(sigma_end)*"\n"
        logstr *= "krange      = "*string(krange)*"\n"
        logstr *= "coherence   = "*string(coherence)*"\n"
        logstr *= "numtrials   = "*string(numtrials)*"\n"

        fprintln("\n"*logstr)
        
        logfile = destination*"_log.txt"
        touch(logfile)
        io = open(logfile, "w")
        write(logfile, logstr)
        close(io)

        spct     = SpectrumObject([sigma_start, sigma_end], [decay_start, decay_end], ["smoothgap"])
        U, S, Vt = create_matrix(rng, n, spct, coher = coherence)
        A        = U*Diagonal(S)*Vt

        data   = Dict()
        means  = Dict()
        stds   = Dict()
        quants = Dict()

        for key in ["rgks", "rgks_angle", "true_angle", "rsvd_angle", "angle_bound"]
            data[key]   = zeros(length(krange), numtrials)
            means[key]  = zeros(length(krange))
            quants[key] = zeros(length(krange), 2)
        end

        totaltrials  = length(krange)*numtrials
        trialcounter = 0

        safe_acos(x) = acos(max(min(x, 1.), 0.))

        for i = 1:length(krange)
            k = krange[i]
            
            for t = 1:numtrials
                trialcounter += 1

                if(trialcounter % 10 == 0)
                    fprintln("  testing matrix "*string(trialcounter)*" out of "*string(totaltrials))
                end

                ov                      = ceil(Int64, .1*k)
                U_appx, S_appx, Vt_appx = rsvd(rng, A, k, oversamp = ov)
                p_rgks                  = qr(Vt_appx[1:k, :], ColumnNorm()).p[1:k]
                Q                       = Matrix{Float64}(qr(A[:, p_rgks]).Q)

                data["rgks"][i, t]        = norm(A - Q*(Q'*A))
                data["rgks_angle"][i, t]  = safe_acos(svd(Vt_appx[1:k, p_rgks]).S[k])
                data["true_angle"][i, t]  = safe_acos(svd(Vt[1:k, p_rgks]).S[k])
                data["rsvd_angle"][i, t]  = safe_acos(svd(Vt[1:k, :]*Vt_appx[1:k, :]').S[k])
                data["angle_bound"][i, t] = data["rgks_angle"][i, t] + data["rsvd_angle"][i, t]
            end

            @save destination*"_data.jld2" krange S data means stds quants
        end

        fprintln("\ncalculating statistics...")
        
        for key in ["rgks", "rgks_angle", "true_angle", "rsvd_angle", "angle_bound"]
            means[key] = vec(mean(data[key], dims = 2))
            stds[key]  = vec(std(data[key], dims = 2))

            for i = 1:length(krange)
                quants[key][i, 1] = quantile(data[key][i, :], .05)
                quants[key][i, 2] = quantile(data[key][i, :], .95)
            end
        end

        @save destination*"_data.jld2" krange S data means stds quants
    end

    run_angle_perturbations(destination, readme, rng, n, decay_start, decay_end, sigma_start, sigma_end, krange, coherence, numtrials)
end

####################################################################
##################### PLOTTING #####################################
####################################################################

@load destination*"_data.jld2" krange S data means stds quants

optimal = [norm(S[(k + 1):end]) for k in krange]
gamma   = [S[k + 1]/S[k] for k in krange]

mfreq      = 10
errbar     = "quantile"   # either "confidence" or "quantile"
confidence = .95

alpha = quantile(Normal(0, 1), 1 - .5*(1 - confidence))

ioff()
fig, (angles, errors) = subplots(1, 2, figsize = (8.5, 4.3))
angles.set_xlabel(L"Approximation Rank ($k$)")
angles.axhline(.5*pi, color = "black", linestyle = "dashed", label = L"\pi / 2")
angles.set_ylim([1.45, 1.01*pi/2])
errors.set_xlabel(L"Approximation Rank ($k$)")

angles.plot(krange, means["rgks_angle"], color = "darkslategray", marker = "o", markevery = mfreq, markerfacecolor = "none", label = L"$\widehat{\varphi}_\mathrm{max}$")
angles.plot(krange, means["true_angle"], color = "forestgreen", marker = "s", markevery = mfreq, markerfacecolor = "none", label = L"$\varphi_\mathrm{max}$")
angles.plot(krange, means["angle_bound"], color = "violet", marker = "D", markevery = mfreq, markerfacecolor = "none", label = L"$\widehat{\varphi}_\mathrm{max} + \theta_\mathrm{max}$")
angles.legend()

errors.plot(krange, means["rgks"]./optimal, color = algcolors["rgks"], marker = algmarkers["rgks"], markevery = mfreq, markerfacecolor = "none", label = "RGKS Subopt.")
errors.plot(krange, gamma, color = "black", linestyle = "dotted", label = L"\sigma_{k + 1}/\sigma_k")
errors.legend()

if(errbar == "quantile")
    angles.fill_between(krange, quants["rgks_angle"][:, 1], quants["rgks_angle"][:, 2], color = "darkslategray", alpha = .2)
    angles.fill_between(krange, quants["true_angle"][:, 1], quants["true_angle"][:, 2], color = "forestgreen", alpha = .2)
    angles.fill_between(krange, quants["angle_bound"][:, 1], quants["angle_bound"][:, 2], color = "violet", alpha = .2)
    errors.fill_between(krange, quants["rgks"][:, 1]./optimal, quants["rgks"][:, 2]./optimal, color = algcolors["rgks"], alpha = .2)
elseif(errbar = "confidence")
    angles.fill_between(krange, (means["rgks_angle"] .- alpha*stds["rgks_angle"])/numtrials, (means["rgks_angle"] .+ alpha*stds["rgks_angle"])/numtrials, color = "darkslategray", alpha = .2)
    angles.fill_between(krange, (means["true_angle"] .- alpha*stds["true_angle"])/numtrials, (means["true_angle"] .+ alpha*stds["true_angle"])/numtrials, color = "forestgreen", alpha = .2)
    angles.fill_between(krange, (means["angle_bound"] .- alpha*stds["angle_bound"])/numtrials, (means["angle_bound"] .+ alpha*stds["angle_bound"])/numtrials, color = "violet", alpha = .2)
    errors.fill_between(krange, (means["rgks"] .- alpha*stds["rgks"])./(optimal*numtrials), (means["rgks"] .+ alpha*stds["rgks"])./(optimal*numtrials), color = algcolors["rgks"], alpha = .2)
else
    throw(ArgumentError("unrecognized error bar type, '"*errbar*"'"))
end

savefig(destination*"_plot.pdf")
close(fig)