Paper guide 
=====
This document provides a guide to the numerical results in [Our paper](https://arxiv.org/pdf/2412.19803), and how they can be obtained from the code in the present repository. 

=== 

- *Visualizing the dynamics of Tsirelson's automaton:* A spacetime history of the dynamics of Tsirelson's automaton (an example of which appears in Fig. 1 of the paper) can be visualized by running 
- *Numerical analysis of error correction in Tsirelson's automaton*: Logical fidelities and relaxation times in Tsirelson's automata can be computed using . This can be used to generate the data given in Figs. 3 and 4 of the paper.  
- *Lemma 5.5:* This lemma is proved with `scripts/nilpotence_tester.jl`. 
- *Brute-force search for low-weight errors:* An explicit brute-force search for small-weight errors that fail the various gadgets can be carried out using. This is used to arrive at the statements about brute-force searches made in Sec. 5.5. 
- *Numerical analysis of error correction in the toric code automaton:* Relaxation times in the toric code automaton can be computed using . This can be used to generate the data in Figs. 10 and 11 of the paper. 