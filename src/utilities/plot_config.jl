using PyPlot
using PyCall

rcParams = PyDict(matplotlib["rcParams"])

rcParams["font.family"] = "serif"
rcParams["font.size"]   = 12

algnames   = ["rsvd_q0", "rsvd_q1", "rgks", "rcpqr", "levg", "dgeqp3"]
alglabels  = Dict("rsvd_q0" => "RSVD (0 PI)", "rsvd_q1" => "RSVD (1 PI)", "rgks" => "RGKS", "rcpqr" => "GRID", "levg" => "LEVG", "dgeqp3" => "DGEQP3")
algcolors  = Dict("rsvd_q0" => "blue", "rsvd_q1" => "magenta", "rgks" => "brown", "rcpqr" => "darkslategrey", "levg" => "forestgreen", "dgeqp3" => "orange")
algmarkers = Dict("rsvd_q0" => "o", "rsvd_q1" => "x", "rgks" => "D", "rcpqr" => "s", "levg" => "v", "dgeqp3" => "+")
