using Random 
using Statistics
using LinearAlgebra 
using JLD2, FileIO

include("gate_builder.jl")
include("gate_applier.jl")
include("helper_functions.jl")

"""
this code tests nilpotence for the gates of the tc automaton using both exhaustive search (which is unfortunately prohibitively slow at l>1 on my laptop) and the recursive error-renormalization method 
"""

function main()

    ### options ### 
    direct_search = ~false # if true, does a direct search through the spacetime volume of level-l gadgets for l ≤ maxlevel. if false, does the method by which errors are fed-forward from one scale to the next (and over-counted in the process)
    clean_input_search = false # if true, does recursive search over inputs with no 1-frame damage and errors in the leading EC + bulk gadget 
    dirty_input_search = true # if true, does recursive search over inputs with arbitrary 1-frame damage on nonlinear points, and errors in the bulk gadget 
    maxlevel = direct_search ? 4 : 10 # maximum number of levels to try before giving up 
    ###############

    gatetypes = ["i" "mx" "my" "tx" "ty"]
    L0s = Dict("i"=>[2 2],"tx"=>[3 2],"ty"=>[2 3],"mx"=>[4 2],"my"=>[2 4])# sizes (in lattice sites) of the scale-0 gates 
    L1s = Dict("i"=>[4 4],"tx"=>[7 4],"ty"=>[4 7],"mx"=>[10 4],"my"=>[4 10]) # sizes (in lattice sites) of the scale-1 gates (can't use PBC for the scale-1 gates since the output damage will be coarse-grained). the gates on the N and E boundaries will be error-free 
    EC_depth = 24 # depth of R0    
    EC_depths = [24 792 25368] # depths of Rl for small l 
    ### direct search through spacetime volume of the circuit ### 

    if direct_search 
        println("doing a direct search up to level $maxlevel")
        # initialize dictionary of level-0 damage types
        damages = Dict()
        for g in 1:length(gatetypes) 
            gtype = gatetypes[g]
            damages[gtype] = even_vertex_sets([L0s[gtype][1],L0s[gtype][2]]) # level-0 damage sets are all the possible damage configurations that can occur when a gate failure occurs on a clean input 
            # each damage pattern given as a list of syndrome coordinates [[x1 y1] [x2 y2] ...]
        end 

        nilpotence = false 
        level = 1 
        while level ≤ maxlevel && ~nilpotence 
            println("investigating nilpotence at level $level")

            # get the depth of the leading error correction part (needed to avoid including errors in the trailing EC step)
            # EC_gates, EC_depth, _ = tc_gate_builder(level-1,"R",false,1,false,false,false,false,true) 
            this_EC_depth = EC_depths[level] # or just use the lookup table 

            # run through the gates 
            failed = false 
            for (g,gtype) in enumerate(gatetypes)
                if failed break end 
                println("gate: $gtype$level")
                # get system sizes of gadgets (including the boundaries)
                Lx = 3^level + 1; Ly = 3^level + 1 
                if gtype == "tx" Lx = 2 * 3^level + 1 end 
                if gtype == "ty" Ly = 2 * 3^level + 1 end 
                if gtype == "mx" Lx = 3^(level+1) + 1 end 
                if gtype == "my" Ly = 3^(level+1) + 1 end 

                init_synds = zeros(Bool,Lx,Ly); init_errs = zeros(Bool,Lx,Ly,2)
                gates, depth, _ , _ = master_gate_builder("tc",1,gtype;all_boundaries=true)

                println("searching spacetime configuration...")
                println("Lx, Ly, T_{search} = $Lx, $Ly, $(depth - EC_depth + 1)")

                for i in 1:size(gates)[1] # gates[gtype] is a matrix and needs to be looped over like this 
                    println(i)
                    if failed break end 
                    gate = gates[i,:]
                    thisgtype = gatetypes[gate[1]]; o = gate[2]; x = gate[3]; y = gate[4]; t = gate[5]
                    if gate[1] == 1 thisgtype = "i" 
                    elseif gate[1] == 2 && o == 1 thisgtype = "tx" 
                    elseif gate[1] == 2 && o == 2 thisgtype = "ty"
                    elseif gate[1] == 3 && o == 1 thisgtype = "mx"
                    elseif gate[1] == 3 && o == 2 thisgtype = "my"
                    end 

                    if t < depth - EC_depth + 1 # don't include errors in the trailing EC step of the gadget currently under investigation (or on the N and E boundaries). 
                        for error in damages[thisgtype] # all the damages that can occur for this gate type 
                            synds = syndrome_tc_gate_applier!(copy(init_synds),gates,level,[t x y],error)
                            if sum(synds) > 0 # error was not cleaned up 
                                failed = true; break 
                            end  
                        end 
                    end 
                    if failed println("no nilpotence up to level $level..."); break end 
                end 
            end 
            level += 1 
        end 
    end 
    
    ### iterative renormalization of errors from one level to the next ### 
    if clean_input_search || dirty_input_search

        println("maximum level to search: $maxlevel")
        gates = Dict(key=>zeros(Int,1,1) for key in gatetypes) # dictionary storing the gate sets of the scale-1 gates; gates are returned by gate_builder.jl as integer valued matrices 

        depths = Dict(key=>1 for key in gatetypes) # dictionary storing depths of aforementioned 
        for gtype in gatetypes        
            gates[gtype], depths[gtype], _ , _ = master_gate_builder("tc",1,gtype;all_boundaries=true)
        end 
        @assert depths["i"] == depths["tx"] == depths["ty"] # current scheme does not need the depth of M to match these values 

        # initialize dictionary of damage types
        # damages[gn] = all types of damage that can appear as output from gate g at level n 
        damages = Dict()
        for g in 1:length(gatetypes) 
            gtype = gatetypes[g]
            for i in 1:maxlevel
                damages[gtype*"$i"] = []
            end 
            damages[gtype*"0"] = even_vertex_sets([L0s[gtype][1],L0s[gtype][2]]) # level-0 damage sets are all the possible damage configurations that can occur when a gate failure occurs on a clean input 
            # each damage pattern given as a list of syndrome coordinates [[x1 y1] [x2 y2] ...]
        end 

        function recursive_search(dirty_input) 
            """
            if dirty_input = true, searches over errors in the bulk gadget, with arbitrary input syndomes on the 1-frame    
            if false, searches over errors in the bulk gadget + leading EC, with clean 1-frame syndromes 
            """

            input_1frame_synds = Dict(gtype => [[]] for gtype in gatetypes)
            output_1frame_synds = copy(input_1frame_synds)
            
            if dirty_input
                # build a list of inputs and outputs for clean M gadgets with input errors on the nonlinear points of the 1-frame 
                nonlinear_vertices_x = [(4, 1) (7, 1) (4, 4) (7, 4)]
                nonlinear_vertices_y = [(1, 4) (1, 7) (4, 4) (4, 7)]
                input_1frame_synds["mx"] = generate_subsets(nonlinear_vertices_x) 
                input_1frame_synds["my"] = generate_subsets(nonlinear_vertices_y)
                output_1frame_synds = copy(input_1frame_synds)

                # generate the output errors one would have under a clean M gadget: 
                for (s,synds) in enumerate(input_1frame_synds["mx"]) 
                    if (4, 1) ∈ synds && (7, 1) ∈ synds 
                        output_1frame_synds["mx"][s] = filter(x -> x ∉ [(4, 1) (7, 1)], output_1frame_synds["mx"][s])
                    end 
                    if (4, 4) ∈ synds && (7, 4) ∈ synds 
                        output_1frame_synds["mx"][s] = filter(x -> x ∉ [(4, 4) (7, 4)], output_1frame_synds["mx"][s])
                    end 
                end 
                for (s,synds) in enumerate(input_1frame_synds["my"]) 
                    if (1, 4) ∈ synds && (1, 7) ∈ synds 
                        output_1frame_synds["my"][s] = filter(x -> x ∉ [(1, 4) (1, 7)], output_1frame_synds["my"][s])
                    end 
                    if (4, 4) ∈ synds && (4, 7) ∈ synds 
                        output_1frame_synds["my"][s] = filter(x -> x ∉ [(4, 4) (4, 7)], output_1frame_synds["my"][s])
                    end 
                end 
            end 

            for level in 1:maxlevel 
                println("\nsending errors from level $(level-1) to $level...")
                for gtype in gatetypes # level-1 gates  

                    # search only over errors that occur at tmin ≤ t ≤ tmax 
                    tmin = dirty_input ? EC_depth + 1 : 1 
                    tmax = depths[gtype] - EC_depth 

                    Lx = L1s[gtype][1]; Ly = L1s[gtype][2] # sizes of level-1 gadget that is to experience level-0 failures

                    lvl1_gates = gates[gtype]

                    for (s,input_1frame_noise) in enumerate(input_1frame_synds[gtype]) # loop over the different input noise 

                        init_synds = falses(Lx,Ly)
                        for syndloc in input_1frame_noise 
                            init_synds[syndloc[1],syndloc[2]] = true 
                        end 

                        for i in 1:size(lvl1_gates)[1] 
                            gate = lvl1_gates[i,:]
                            thisgtype = gatetypes[gate[1]]; o = gate[2]; x = gate[3]; y = gate[4]; t = gate[5]
                            if t ≥ tmin && t ≤ tmax && x < Lx && y < Ly # errors come in only during the application of the bulk gadget. (note: if errors happen in the leading EC, self-similar damage is generated) 

                                if gate[1] == 1 thisgtype = "i" 
                                elseif gate[1] == 2 && o == 1 thisgtype = "tx" 
                                elseif gate[1] == 2 && o == 2 thisgtype = "ty"
                                elseif gate[1] == 3 && o == 1 thisgtype = "mx"
                                elseif gate[1] == 3 && o == 2 thisgtype = "my"
                                end 

                                for error in damages[thisgtype*"$(level-1)"] # all the damages that can occur for this gate type  
                                    synds = syndrome_tc_gate_applier!(copy(init_synds),lvl1_gates,1,[t x y],error)
                                    # now extract the output damage 
                                    renormalized_synds = coarse_grain_damage(synds) # an array (possibly empty) of xy coordinate tuples 
                                    renormalized_damage = copy(renormalized_synds)
                                    if gtype ∈ ["mx" "my"] && dirty_input # compare to synds that you'd get for the ideal gate 
                                        ideal_output_synds = zeros(Bool,Lx,Ly)
                                        if length(output_1frame_synds[gtype][s]) > 0
                                            for (i,j) in output_1frame_synds[gtype][s] # output_1frame_synds is a list of tuples 
                                                ideal_output_synds[i,j] = true 
                                            end 
                                        end 
                                        ideal_renormalized_synds = coarse_grain_damage(ideal_output_synds)
                                        # find the mod-2 sum of anyon locations in the renormalized syndromes and the error-free renormalized syndromes 
                                        renormalized_damage = setdiff(union(renormalized_synds,ideal_renormalized_synds), intersect(renormalized_synds,ideal_renormalized_synds)) 
                                    end 

                                    # determine if the output damage has already been accounted for 
                                    in_damage = false 
                                    for dam in damages[gtype*"$level"]
                                        if Set(renormalized_damage) == Set(dam) # needed to check if damages differ by a permutation; taking sets gets rid of the ordering that lists have 
                                            in_damage = true 
                                            break 
                                        end 
                                    end 
                                    if ~in_damage
                                        push!(damages[gtype*"$level"],renormalized_damage) 
                                    end 

                                    # check for self-similar errors 
                                    if renormalized_damage == error && length(error) > 0 && gtype == thisgtype 
                                        println("found self-similar error!")
                                        println("input_1frame_noise = $(input_1frame_synds[gtype][s])") 
                                        println("output_1frame_noise = $(output_1frame_synds[gtype][s])") 
                                        println("gate type = $gtype")
                                        println("ren damage = ",renormalized_damage)
                                        println("error = ",error)
                                        println("(t,x,y,o) = $t, $x, $y, $o")
                                        println("renormalized damage = ",renormalized_damage)
                                        return                                         
                                    end     
                                end 
                            end 
                        end 
                    end 
                    damages[gtype*"$(level)_weights"] = count_error_weights(damages[gtype*"$level"])
                end 

                # check to see if the errors in all types of gates have been cleaned up 
                cleaned_up = true 
                for gtype in gatetypes
                    ndamages = length(damages[gtype*"$level"])
                    if ndamages > 1
                        cleaned_up = false 
                    end 
                    println("number of $gtype failures = $(ndamages-1)")
                    println("$gtype damages:")
                    for dam in damages[gtype*"$level"] println(dam) end 
                    println("weight distribution = ",damages[gtype*"$(level)_weights"])
                end 
                if cleaned_up
                    println("terminating at level $level"); break 
                end 
            end 
        end 
    
        if clean_input_search
            println("doing search over leading EC + bulk gadget with clean input...")
            recursive_search(false)
        end 
        if dirty_input_search 
            println("doing search over bulk gadget with dirty input...")
            recursive_search(true)
        end 
    end 
end 

main()
