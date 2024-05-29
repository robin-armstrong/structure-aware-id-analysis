This repository contains Julia code to replicate numerical experiments in "Structure Aware Analyses and Algorithms for Interpolative Decompositions" by R. Armstrong, A. Buzali, and A. Damle. To begin working with this code, run
```
git clone https://github.com/robin-armstrong/structure-aware-id-analysis.git
cd structure-aware-id-analysis
julia
```
from a terminal. Then enter the following commands in Julia:
```
import Pkg
Pkg.activate(".")
Pkg.instantiate()
```
Experiments are located within the `src/experiments` directory, and are intented to be run from the top level. Thus, to run a particular experiment, enter
```
include("src/experiments/<experiment-name>/run.jl")
```
in Julia, or enter
```
julia --project=. src/experiments/<experiment-name>/run.jl
```
from a shell at the top level.