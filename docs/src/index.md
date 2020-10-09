# Welcome to pfem4ec

`pfem4ec` is a fast [Julia](https://julialang.org/) in-core solver to compute effective electrical conductivity of heterogeneous materials from raw images using Pixel-Based Finite Element Method (PFEM). The software was implemented in Julia because the resulting code is as simple as it would be if it were written using a scripting language (such as Python or Matlab/Octave). On the other hand, as Julia takes advantage of a just-in-time (JIT) technology, the codes can also be as efficient as codes written using compiled programming languages (such as in C/C++ or Fortran).

# Performance Tips

## Requirements

The following libraries are required, which are likely already installed on your system:
+ [JSON.jl](https://github.com/JuliaIO/JSON.jl)
+ [SparseArrays.jl](https://github.com/JuliaLang/julia/blob/master/stdlib/SparseArrays/src/SparseArrays.jl)
+ [LinearAlgebra.jl](https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/LinearAlgebra.jl)


## Installation

Download source code in the terminal:
```bash
git clone https://github.com/LCC-UFF/pfem4ec
```

## How to run

### Linux or Windows

Firstly, it is important to state that, even though two input files are needed, when running the code, `pfem4ec.jl` takes a single input, a string containing the name of the model (without file extensions), hence the need to give the same name to correlated RAW and JSON files.  

Users can run `pfem4ec.jl` directly from the terminal, but they need to bear in mind that Julia compiles the code prior to running, meaning that each call will take some extra time to compile (for relatively small models, this usually takes more time than the analysis itself). Certainly, if no changes are made to `pfem4ec.jl`, there’s no need to compile the same code over and over again. In light of this, the authors suggest that users run `pfem4ec.jl` from within the Julia environment, as follows.

a. Start the Julia REPL:

```bash
julia
```
P.S.: Make sure that you’re on the directory where `pfem4ec.jl` and the input files are placed.

b. Include `pfem4ec.jl`:

```julia
include(“pfem4ec.jl”)
```
P.S.: The first call will compile and run, all subsequent calls will just run, considerably 	reducing overall elapsed time to obtain results.

c. Run the `pfem4ec` function:
```julia
pfem4ec(“model_filename”)
```

d. To obtain time metrics for each time `pfem4ec` is run, simply add the Julia macro @time before the call presented above.
```julia
@time pfem4ec(“model_filename”)
```

## Example usage

For now, `pfem4ec.jl` output consists on feedback of the analysis that is printed on the terminal that called the program. Users are provided with an estimation of allocated memory, convergence metrics for the PCG numerical method, a stamp for each function called by `pfem4ec` along the process and, finally, the resulting effective homogenized conductivity matrix. <a href="#figure1">Figure 1</a> shows the output for the 5x5 example model adopted in our [user guide](https://github.com/LCC-UFF/pfem4ec/blob/master/docs/user_guide_pfem.pdf).

<a name="figure1"><div id="figure1"></div></a>
<p align="center">
  <img src="https://github.com/LCC-UFF/pfem4ec/blob/master/images/julia_pfem4ec_output.png?raw=true">
</p>
<p align="center">Figure 1: pfem4ec output for analysis on 5x5 model.</p>