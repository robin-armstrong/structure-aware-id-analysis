using LinearAlgebra
using SpecialFunctions
using StatsBase
using Hadamard
using Random
using PyPlot
using JLD2

include("../../algorithms/rsvd.jl")
include("../../utilities/plot_config.jl")

####################################################################
##################### SCRIPT PARAMETERS ############################
####################################################################

# DEFAULT PARAMETER SETTINGS:
#
# rng          = MersenneTwister(1)
# n            = 1024
# k            = 70
# decay_start  = 5
# decay_end    = 80
# sigma_start  = 1.
# sigma_end    = 1e-12
# noise        = 0.01
# numtrials    = 250

# path prefix for all output generated by this script
destination = "src/experiments/subspace_stats/substats"
readme      = "Different subspace error metrics for the randomized SVD."

rng          = MersenneTwister(1)   # random souce, set for reproducibility
n            = 1024
k            = 70
decay_start  = 5
decay_end    = 80
sigma_start  = 1.
sigma_end    = 1e-12
noise        = 0.01
numtrials    = 250

plot_only = false

####################################################################
##################### DATA GENERATION ##############################
####################################################################

function fprintln(s)
    println(s)
    flush(stdout)
end

if(!plot_only)
    function run_subspace_stats(destination, readme, rng, n, k, decay_start, decay_end, sigma_start, sigma_end, noise, numtrials)
        logstr  = readme*"\n\n"
        logstr *= "rng         = "*string(rng)*"\n"
        logstr *= "n           = "*string(n)*"\n"
        logstr *= "k           = "*string(k)*"\n"
        logstr *= "decay_start = "*string(decay_start)*"\n"
        logstr *= "decay_end   = "*string(decay_end)*"\n"
        logstr *= "sigma_start = "*string(sigma_start)*"\n"
        logstr *= "sigma_end   = "*string(sigma_end)*"\n"
        logstr *= "noise       = "*string(noise)*"\n"
        logstr *= "numtrials   = "*string(numtrials)*"\n"

        fprintln("\n"*logstr)
        
        logfile = destination*"_log.txt"
        touch(logfile)
        io = open(logfile, "w")
        write(logfile, logstr)
        close(io)
        
        # setting up the singular spectrum of the test matrix
        S = zeros(n)
        S[1:(decay_start-1)] .= sigma_start
        S[(decay_end+1):end] .= sigma_end

        decay_points  = decay_end - decay_start + 1
        t             = range(-2, 2, decay_points)
        decay_curve   = .5*(erfc.(t) .+ 1)
        beta          = log(sigma_start/sigma_end)/log(2.981)
        decay_curve .^= beta
        decay_curve  *= sigma_start/decay_curve[1]
        
        S[decay_start:decay_end] = decay_curve

        # other components of the test matrix SVD
        U  = Matrix{Float64}(qr(randn(rng, n, n)).Q)
        Vt = Dict()

        data = Dict()

        for coher in ["incoherent", "random", "coherent"]
            data[coher] = Dict()

            for key in ["nrm", "max", "med", "avg"]
                data[coher][key] = zeros(numtrials)
            end
        end

        totaltrials  = 3*numtrials
        trialcounter = 0

        for coher in ["incoherent", "random", "coherent"]
            if(coher == "incoherent")
                W       = hadamard(n)/sqrt(n) + noise*randn(rng, n, n)
                svdobj  = svd(W')
                Vt_true = svdobj.U*svdobj.Vt
            elseif(coher == "coherent")
                W       = I(n)[:, randperm(rng, n)] + noise*randn(rng, n, n)
                svdobj  = svd(W')
                Vt_true = svdobj.U*svdobj.Vt
            else
                Vt_true = Matrix{Float64}(qr(randn(rng, n, n)).Q)
            end

            Vt[coher] = Vt_true

            A      = U*Diagonal(S)*Vt_true
            P_true = Vt_true[1:k, :]'*Vt_true[1:k, :]

            for t = 1:numtrials
                trialcounter += 1

                if(trialcounter % 10 == 0)
                    fprintln("  testing matrix "*string(trialcounter)*" out of "*string(totaltrials))
                end

                Vt_appx = rsvd(rng, A, k, oversamp = 5).Vt
                P_appx  = Vt_appx[1:k, :]'*Vt_appx[1:k, :]
                E       = P_true - P_appx

                data[coher]["nrm"][t] = sqrt(1 - svd(Vt_appx*Vt_true[1:k, :]').S[k]^2)
                data[coher]["max"][t] = maximum(abs.(E))
                data[coher]["med"][t] = median(abs.(E))
                data[coher]["avg"][t] = mean(abs.(E))

                @save destination*"_data.jld2" U S Vt k data
            end
        end
    end

    run_subspace_stats(destination, readme, rng, n, k, decay_start, decay_end, sigma_start, sigma_end, noise, numtrials)
end

####################################################################
##################### PLOTTING #####################################
####################################################################

@load destination*"_data.jld2" U S Vt k data

ioff()
fig, (coherent, incoherent, random) = subplots(1, 3, figsize = (15, 4))
coherent.set_title("Coherent")
incoherent.set_title("Incoherent")
random.set_title("Random")

tile_list = [coherent, incoherent, random]
tile_keys = ["coherent", "incoherent", "random"]
bins      = 20
alpha     = .7
colors    = ["C0", "C1", "C2", "C3"]

for tile in tile_list
    tile.set_xscale("log")
    tile.set_xlim([1e-4, 2])
    tile.set_ylim([0, 80])
end

for i = 1:3
    tile_list[i].hist(data[tile_keys[i]]["nrm"], bins, alpha = alpha, color = "C0", label = L"$|| P - \widehat{P}_t ||_2$")
    tile_list[i].hist(data[tile_keys[i]]["max"], bins, alpha = alpha, color = "C1", label = L"$\mathrm{max} | P - \widehat{P} |_{i,\, j}$")
    tile_list[i].hist(data[tile_keys[i]]["med"], bins, alpha = alpha, color = "C2", label = L"$\mathrm{med} | P - \widehat{P} |_{i,\, j}$")
    tile_list[i].hist(data[tile_keys[i]]["avg"], bins, alpha = alpha, color = "C3", label = L"$\mathrm{avg} | P - \widehat{P} |_{i,\, j}$")
end

random.legend()

savefig(destination*"_plot.pdf")
close(fig)