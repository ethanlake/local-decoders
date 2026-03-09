local-decoders
=====
Code relating to the construction and analysis of 1) Tsirelson's 1d memory [1] and its generalizations, and 2) a local Tsirelson-inspired decoder for the 2d toric code. All work done jointly with Shankar Balasubramanian and Margarita Davydova. The main paper is Ref. [2] below.

project structure
=====
```
src/                          # 
  LocalDecoders.jl            # module definition
  gate_builder.jl             # builds circuits for below-listed codes
  gate_applier.jl             # applies circuits to states in the presence of various types of noise
  helper_functions.jl         # miscellaneous auxillary functions
  simulation_parameters.jl    # repository of parameter values (thresholds, default samples, etc.)
scripts/                      # 
  circuit_simulator.jl        # main simulation code (see below)
  nilpotence_tester.jl        # verifies nilpotence of toric code automaton gadgets
  build_gates.jl              # standalone gate building and saving
plotters/                     # 
  circuit_plotter.py          # plots .jld2 output from circuit_simulator.jl
  1d_circuit_visualizer.py    # draws 1d circuit / spacetime diagrams
  2d_circuit_visualizer.py    # interactive toric code automaton visualizer
  scaling_plotter.py          # interactive scaling collapse plots
data/                         # simulation output (.jld2 files)
visuals/                      # visualizations and animations
```

`scripts/circuit_simulator.jl` can compute the following quantities for the currently supported codes listed below:
- Relaxation times (mode: `trel`)
- Time evolution of logical fidelity (mode: `Ft`)
- Spacetime trajectories under custom noise realizations (mode: `hist`)
- Searches over low-weight noise realizations that lead to logical errors (mode: `noise_search`)
- Logical error rate in the code capacity scenario (mode: `offline`)

setup
=====
```bash
# first-time setup: install dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'
```

guide to reproducing results in [2]
=====

For the full list of adjustable parameters in `scripts/circuit_simulator.jl`, see the `main()` block of that file. 
### visualize action of Tsirelson's automaton (Fig. 1 of [2])
```bash
julia --project scripts/circuit_simulator.jl --mode hist --model rep --gate decoder -l 2 # build a level-2 decoder (the minimal circuit on 3^3 bits that cleans up any input error), and simulate its action on a random initial state 
python3 plotters/1d_circuit_visualizer.py -fin data/rep2_gates_decoder.jld2 -history data/rep_l2_hist_decoder_test_unbiased.jld2 # plot the resulting spacetime history
```

### compute logical fidelity for Tsirelson's automaton (Fig. 3 of [2])
```bash 
julia --project scripts/circuit_simulator.jl --mode Ft --model rep -l 2 --pmin .012 --pmax .018 --nps 8 # compute the fidelity 
python3 plotters/circuit_plotter.py -fin data/rep_l2_Ft_R_test_unbiased.jld2 # plot the fidelity 
``` 

### compute relaxation time in a similar setting (Fig. 4 of [2])
```bash 
julia --project scripts/circuit_simulator.jl --mode trel --model rep -l 2 --pmin .008 --pmax .018 --nps 6 # generate the relaxation time data 
python3 plotters/circuit_plotter.py -fin data/rep_l2_trel_R_test_unbiased.jld2 # plot the result 
```

### run nilpotence tester (used to prove Lemma 5.5 of [2])
```bash 
julia --project scripts/nilpotence_tester.jl --clean-input-search --maxlevel 3
``` 

### visualize the action of the toric code automaton
```bash 
julia --project scripts/circuit_simulator.jl --mode hist --model tc --gate R -l 1 # build a level-1 R gate and compute how a level-1 corner R gets acted on by it (to change the initial condition, consult the code block on line 501)
python3 plotters/2d_circuit_visualizer.py -hist data/tc_l1_hist_R_test.jld2 -gates data/tc_l1_gates_R.jld2 # build an interactive visualization of the error correcting dynamics
``` 

### compute the relaxation time for the toric code automaton (Fig. 10 of [2])
```bash
julia --project scripts/circuit_simulator.jl --mode trel --model tc -l 0 --pmin 5.5e-5 --pmax 3e-4 # estimate the relaxation time at the lowest level (set measurementnoise = true or gadgetnoise = true to reproduce Fig. 11) 
python3 plotters/circuit_plotter.py -fin /Users/elake/Research/GitHub/local-decoders/scripts/../data/tc_l0_trel_R_test.jld2 # plot the result 
```

### perform a brute-force search for small-weight errors (see the discussion in Sec. 5.5 of [2])
```bash 
julia --project scripts/circuit_simulator.jl --mode noise_search --model tc --gate R -l 0 # check that no level-0 errors fail R0 (to check different failure conditions, modify the code in the block starting on line 546)
```

### see all options
```bash 
julia --project scripts/circuit_simulator.jl --help
```

### REPL example
```julia
julia --project
julia> using LocalDecoders
julia> include("scripts/circuit_simulator.jl")
julia> main(mode="Ft", model="rep", l=2, test=true)
```

currently supported codes
=====
- 1d Tsirelson's automaton with base codes:
    - 3- and 5-bit repetition codes
    - a particular [5, 2, 3] code
- A 2d version of Tsirelson's automaton based on 9-bit reptition codes
- our 2d toric code automaton

noise models
======
- "wire noise" (default): Gates are implimented reliably, and iid noise is applied to each spin after the gates in a given layer have acted. For all but the toric code, the noise replaces a given spin by $\pm1$ with probability $p(1\pm \eta)/2$, and does nothing with probability $p$. for the toric code, a given spin is flipped with probability $p$.
- "measurement noise": Each syndrome measurement independently gives the wrong outcome with probability $p$.
- "gadget noise": Each gadget fails independently with probability $p$. When a gadget fails, its spacetime support incurrs a randomly chosen error from a pre-defined set (toric code), or the spins in the support are randomized (other codes). For the toric code, the supports of each gadget and the errors that occur when they fail can be modified in `scripts/circuit_simulator.jl`.
- Each type of noise model also allows for a specific spacetime history of noise to be applied by appropriately modifying `noise_hist` in `scripts/circuit_simulator.jl`.


references
======
1. [Tsirelson's original paper](https://link.springer.com/book/10.1007/BFb0070079) (see page 26)
2. [Our paper](https://arxiv.org/pdf/2412.19803)
