local-decoders
=====
code relating to the construction and analysis of 1) Tsirelson's 1d memory [1] and its generalizations, and 2) a local Tsirelson-inspired decoder for the 2d toric code.

main files
=====
- `circuit_simulator.jl` main simulation code capable of computing the following quantities for the codes listed below: 
    - time evolution of logical fidelity
    - relaxation times
    - spacetime trajectories under custom noise realizations
    - searches over low-weight noise realizations that lead to logical errors
    - statistical properties in long-time steady states (various moments of magnetization, statistics of domain wall sizes, etc.)  
- `gate_builder.jl` builds circuits for below-listed codes 
- `gate_applier.jl` applies circuits to states in the presence of various types of noise 
- `nilpotence_tester.jl` used to verify the "nilpotence" property of gadgets used in the toric code automaton
- `helper_functions.jl` contains miscellaneous auxillary functions 
- `circuit_plotter.py` plots data output (as `.jld2` files) from `circuit_simulator.jl`
- `1d_circuit_visualizer.py` draws diagrams of circuits (and spacetime spin configurations) for 1d codes 
- `2d_circuit_visualizer.py` an interactive visualizer for the toric code automaton 

currently supported codes 
=====
- 1d Tsirelson's automaton with base codes
    - 3- and 5-bit repetition codes 
    - a [5, 2, 3] code
- A 2d version of Tsirelson's automaton based on 9-bit reptition codes
- 2d toric code

noise models
====== 
- "wire noise" (default): gates are implimented reliably, and iid noise is applied to each spin after the gates in a given layer have acted. for all but the toric code, the noise replaces a given spin by $\pm1$ with probability $p(1\pm \eta)/2$, and does nothing with probability $p$. for the toric code, a given spin is flipped with probability $p$. 
- "measurement noise": each syndrome measurement independently gives the wrong outcome with probability $p$.
- "gadget noise": each gadget fails independently with probability $p$. when a gadget fails, its spacetime support incurrs a randomly chosen error from a pre-defined set. the supports of each gadget and the errors that occur when they fail can be modified in `circuit_simulator.jl`.
- each type of noise model also allows for a specific spacetime history of noise to be applied by appropriately modifying `noise_hist` in `circuit_simulator.jl`. 

usage notes
====== 
- as written, `circuit_simulator.jl` is designed to be run in the REPL; changes to parameter values should be made by making direct edits to the relevant lines.
- the function `parameter_repository` in `circuit_simulator` contains numerical values for the thresholds of various codes.
- the fps and length of animations produced by `2d_circuit_visualizer.py` needed to be specified by hand by editing the relevant section of code.

references
======
1. [Tsirelson's original paper](https://link.springer.com/book/10.1007/BFb0070079) (see page 26)



  
