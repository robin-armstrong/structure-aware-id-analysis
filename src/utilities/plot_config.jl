using PyPlot
using PyCall

rcParams = PyDict(matplotlib["rcParams"])
rcParams["font.family"] = "serif"

algnames   = ["rsvd_q0", "rsvd_q1", "rgks", "rid", "levg"]
alglabels  = Dict("rsvd_q0" => "RSVD (0 iterations)", "rsvd_q1" => "RSVD (1 iteration)", "rgks" => "RGKS", "rid" => "RID", "levg" => "LEVG")
algcolors  = Dict("rsvd_q0" => "blue", "rsvd_q1" => "magenta", "rgks" => "brown", "rid" => "darkslategrey", "levg" => "forestgreen")
algmarkers = Dict("rsvd_q0" => "o", "rsvd_q1" => "x", "rgks" => "D", "rid" => "s", "levg" => "v")
