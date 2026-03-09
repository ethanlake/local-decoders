module LocalDecoders

using Random, Statistics, LinearAlgebra
using JLD2, FileIO
using Combinatorics

include("helper_functions.jl")
include("gate_builder.jl")
include("gate_applier.jl")
include("simulation_parameters.jl")

# gate_builder.jl
export get_gate_size, master_gate_builder, permute_circuit
export tsirelson_gate_builder, twod_tsirelson_gate_builder, tc_gate_builder

# gate_applier.jl
export master_gate_applier, tc_gate_applier!, syndrome_tc_gate_applier!
export tsirelson_gate_applier, twod_tsirelson_gate_applier
export get_domainwall_stats, get_corrs
export maj, correct, energy, xy_maj
export detect_logical_failure, detect_both_logical_failures, average_logicals
export recursive_majority
export H, G

# helper_functions.jl
export generate_subsets, even_vertex_sets
export coarse_grain_damage, errs_to_synds, synds_to_links, get_noise_hist, count_error_weights

# simulation_parameters.jl
export parameter_repository

end
