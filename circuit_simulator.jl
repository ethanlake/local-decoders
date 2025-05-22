using JLD2
using Combinatorics 

include("gate_builder.jl")
include("gate_applier.jl")
include("helper_functions.jl")
include("simulation_parameters.jl")

"""
this code runs monte carlo simulations of various different types of self-correcting circuit dynamics.
can compute: 
    * time dependence of logical fidelity 
    * relaxation times 
    * statistical properties and correlation functions in long-time steady states 

arguments and options given in the first section of main()
"""

function decode(model,state,synds;decoder_gates=[])
    if model ∈ ["rep" "rep_5bit" "twod_rep"]
        return mean(state) < .5 ? 1 : 0 # gobal majority vote is ofc not the same as actually running the decoder --- but it's faster and emperically is equivalent in the TDL. can instead use recursive_majority() to use the exact decoding 
    elseif model == "k2" # [5,2,3] code -- decode just using the code distance of 3^l (always using the all-0s codeword)
        return sum(state) < (3^(log(5,size(state)[1])) -1)/2 ? 1 : 0 
    elseif model == "tc" # toric code -- decode by running ideal decoder and then measuring fixed logicals; harrington's heuristic of straight logicals unsurprisingly doesn't work that well  
        decoded_state, decoded_synds = master_gate_applier(model,state,synds,decoder_gates,1,[0],0,false,false,false,false); @assert ~any(decoded_synds) 
        return detect_logical_failure(decoded_state) ? 0 : 1          
    end            
end 

function noisy_state(p,L,d)
    """
    generates a random state subject to iid errors of strength p.
    p: density of errors
    L: system size
    d: number of dimensions (1 or 2)
    """
    ind(i,L) = i > L ? i - L : i 
    state = d == 1 ? falses(L) : falses(L,L,2) 
    synds = d == 1 ? falses(L) : falses(L,L)
    if d == 1 
        for i in 1:L
            if rand() < p
                state[i] = true
                synds[i] ⊻= true; synds[ind(i+1,L)] ⊻= true
            end 
        end
    elseif d == 2
        for i in 1:L, j in 1:L, o in 1:2
            if rand() < p
                state[i,j,o] = true
                if o == 1 
                    synds[i,j] ⊻= true; synds[ind(i+1,L),j] ⊻= true
                else
                    synds[i,j] ⊻= true; synds[i,ind(j+1,L)] ⊻= true
                end
            end 
        end
    else 
        error("dimension must be 1 or 2")
    end

    return state, synds
end 

function main() 
    ##### arguments and options #####
    
    """
    basic structural options 
    mode: what the code calculates. one of the following: 
        * "Ft": time evolution of logical fidelity (starting in 0^L codeword by default)
        * "trel": relaxation time of logical information 
        * "stats": collects statistics about the long time (t >> trel) steady state. used when computing critical properties. 
        * "hist": just runs a single sample of the dynamics; used for visualization / testing specific noise configurations 
        * "offline": computes the logical error rate under offline decoding
    model ∈ ["rep" "rep_5bit" "k2" "tc" "twod_rep"]
    gate ∈ ["Z" "R" "decoder" "Y" "id" ...]
    test: if true, uses less samples / less periods to speed things up 
    l: level, equals log_n(system size)-1 for EC gates and log_n(system_size) for other gates (in our current (somewhat unfortunate) conventions)
    """
    mode = "offline" #"stats" # "hist" # "stats" #"trel" 
    model = "tc" # "rep"
    gate = "R"
    test = true 
    l = 1 # pool
    n = model ∈ ["rep" "tc" "twod_rep"] ? 3 : 5
    Lx, Ly = get_gate_size(l,model,gate); L = Lx  
    
    """
    noise options 
    squarenoise: if true, applies errors by randomly applying plaquette errors (rather than link errors). 
    measurementnoise: if true, noise comes only from failures in the syndrome measurements
    bias: bias of iid noise (not used for TC)
    gadgetnoise: if true, gadgets fail (all bits in the output of a failed gadget are randomized). if false, "wires" between gadgets fail.
    nps: number of noise values to sample 
    logdist: if true, spaces out noise values evenly on a log scale. if false, spaces out evenly on a linear scale. 
    """
    squarenoise = false 
    measurementnoise = false 
    bias = 0.  
    gadgetnoise = false 
    nps = mode == "trel" ? (l > 1 ? 7 : 9) : (mode ∈ ["Ft" "offline"] ? 10 : (mode == "hist" ? 1 : (mode == "stats" ? 11 : 1)))  
    logdist = mode == "trel" ? true : false 

    """
    options for analyzing steady states of 1d codes  
    calc_corrs: if true, calculates 2-point spin correlation functions at all unique spacetime points in a single floquet period 
    measure_corrs_multiple: only measure correlation functions / domain wall statistics once every measure_corrs_multiple number of steps; want to be larger than a correlation time which should perhaps be dynamically determined using earlier MC code... 
    calc_dw_stats: if true, calculates statistical distribution of domain wall sizes in the long-time steady state 
    rollstate: if true, performs a random cyclic spatial shift on the state in between floquet periods 
    """
    calc_corrs = ~true
    measure_corrs_multiple = 500 
    calc_dw_stats = ~true 
    rollstate = ~true 

    """
    file input / output options 
    adj: distinguishing information about input file (gates)
    out_adj: ... about output file (simulation data)
    """
    adj = "" 
    out_adj = "_test" 

    ####################################### 

    # give the output file some identifying information 
    if model ∉ ["tc"] out_adj *= (bias == 0 ? "_unbiased" : (bias == 1 ? "_biased" : "_medbiased")) end 
    if calc_corrs out_adj *= "_corrs" end 
    if calc_dw_stats out_adj *= "_dwstats" end 
    if model == "tc" 
        if squarenoise out_adj *= "_sqnoise" end 
        if gadgetnoise out_adj *= "_gadnoise" end 
        if measurementnoise out_adj *= "_measmntnoise" end 
    end 
    
    gadgetnoise_dict = []
    if gadgetnoise # determine failure modes and supports of primitive gadgets (for tc automaton)
        gadgetnoise_dict = Dict((i,o) => [] for i in 1:3, o in 1:2) # first index 1,2,3 is gate type I,T,M; second index 1,2 is orientation

        nilpotence_support = true 
        if nilpotence_support
            # "fattened" gates with nilpotence-proof support: 
            gadgetnoise_dict[1,1] = even_vertex_sets([2,2]) # identity 
            gadgetnoise_dict[2,1] = even_vertex_sets([3,2]) # Tx 
            gadgetnoise_dict[2,2] = even_vertex_sets([2,3]) # Ty 
            gadgetnoise_dict[3,1] = even_vertex_sets([4,2]) # Mx 
            gadgetnoise_dict[3,2] = even_vertex_sets([2,4]) # My  
        else 
            # "skinny" gates with support only on those vertices which are either measured or can have their syndrome values changed by the gadget 
            gadgetnoise_dict[1,1] = even_vertex_sets([1,1]) # identity 
            gadgetnoise_dict[2,1] = even_vertex_sets([3,1]) # Tx (T can fail on any of the three vertices on which it either measures or (has the possibility of) creating them)
            gadgetnoise_dict[2,2] = even_vertex_sets([1,3]) # Ty 
            gadgetnoise_dict[3,1] = even_vertex_sets([2,1]) # Mx (M can fail on any of the two vertices on which ...)
            gadgetnoise_dict[3,2] = even_vertex_sets([1,2]) # My  

            # shift the supports so that the errors only occur on the desired vertices in the `middle` of M: 
            for i in 2:length(gadgetnoise_dict[3,1]) for j in 1:length(gadgetnoise_dict[3,1][i]) gadgetnoise_dict[3,1][i][j] = (gadgetnoise_dict[3,1][i][j][1] + 1, gadgetnoise_dict[3,1][i][j][2]) end end # sadly tuples are immutable in julia... 
            for i in 2:length(gadgetnoise_dict[3,2]) for j in 1:length(gadgetnoise_dict[3,2][i]) gadgetnoise_dict[3,2][i][j] = (gadgetnoise_dict[3,2][i][j][1], gadgetnoise_dict[3,2][i][j][2] + 1) end end  
        end 
    end 

    ps,pc,samps,thermal_periods,periods = parameter_repository(model,mode,l,nps,logdist,bias,gadgetnoise,measurementnoise,test)
    fout = "data/$(model)_l$(l)_$(mode)_$(gate)"*(mode == "trel" ? "_$(periods)pers" : "")*"$(out_adj).jld2"
    if mode == "trel"
        fout = "data/$(model)_l$(l)_$(mode)_$(gate)$(out_adj).jld2"
    end 
    if mode == "stats"
        fout = "data/$(model)_l$(l)_$(mode)_$(gate)$(out_adj).jld2"
    end 
    println("will save data as: $fout")
    
    # load / generate gate data 
    gates = []; depth = 24
    # check to see if gate data already exists: 
    fin = "data/$(model)_l$(l)_gates_$(gate)$(adj).jld2" 
    if isfile(fin)
        println("loading gates from file: $fin")
        f = jldopen(fin,"r")
        gates = read(f,"gates"); depth = read(f,"depth")
        println("circuit depth = $depth")
        close(f) 
    else # if it doesn't, make it ourselves 
        println("building gates...")
        gates, depth, Lx, Ly = master_gate_builder(model,l,gate) # to do open bconds for tc need to call all_boundaries = true as an optional argument 
        L = Lx # just for more concise notation in the 1d case 

        # and then save the gate data so that we don't have to re-generate it next time 
        println("writing gate data to file: $fin")
        f = jldopen(fin,"w")
        write(f,"model",model); write(f,"gates",gates); write(f,"gate",gate); write(f,"n",n); write(f,"l",l); write(f,"depth",depth); write(f,"L",Lx); write(f,"Lx",Ly); write(f,"Ly",Ly) 
        close(f)    
    end 
    println("system size = ($Lx,$Ly)")
    ngates = size(gates)[1]

    # all the things that we may have reason to measure 
    datakeys = ["trels" "avg_trels" "Fs" "Ps" "|M|s" "Ms" "floquet_|M|s" "floquet_Ms" "chis" "floquet_chis" "binds" "floquet_binds" "average_history" "dw_stats" "corrs" "tcorrs" "noise_hist" "err_hist" "synd_hist" "floquet_hamming_statistics" "Es" "floquet_Es" "Cs" "floquet_Cs" "dEdps" "dMdps" "halfps" "offline_error_rate"] 
    """
    hamming_statistics: histogram of hamming weight distributions 
    M: magnetization
    chi: susceptibility of magnetization 
    bind: binder cumulant of magnetization (absolute value or signed magnetization, depending on whether noise is biased)
    dMdps: derivative of ⟨|M|⟩ wrt p 
    P: codeword density: ∑_blocks Θ(spins at block is a codeword) / (# of blocks)
    chiP: susceptibility of codeword density 
    bindP: binder cumulant of codeword density 
    Es: total energy (number of violated stabilizers)
    Cs: "specific heat" (susceptibility of total energy)
    dEdps: derivative of ⟨E⟩ wrt p
    """
    
    data = Dict{String, Any}(); for key in datakeys data[key] = zeros(1) end 

    codeword = []; codesynds = []
    if model ∈ ["rep" "rep_5bit" "k2"] codeword = falses(L)  
    elseif model == "twod_rep" codeword = falses(Lx,Ly)
    elseif model == "tc" codeword = falses(Lx,Ly,2); codesynds = falses(Lx,Ly) 
    end 

    ### relaxation times ### 
    if mode == "trel"
        println("calculating trel for $samps samples...")

        data["trels"] = zeros(nps,samps) # full distribution of relaxation times 
        data["avg_trels"] = zeros(nps) # the average thereof (median and mean will agree as long as ⟨trel⟩ sufficiently bigger than a single floquet period)
    
        maxsteps = 10000000 # max number of floquet periods to run the simulation for 

        for (pind,p) in enumerate(ps)
            println("p = $p")

            thistrel = 0
            for samp in 1:samps 
                step = -1 
                state = copy(codeword); synds = copy(codesynds) # what we will be evolving 

                F = 1 # F = 0 is the termination condition for trel 
                while F == 1 && step < maxsteps 
                    
                    # apply one period of the circuit to state 
                    state, synds = master_gate_applier(model,state,synds,gates,1,[p],bias,gadgetnoise_dict,squarenoise,measurementnoise,false)
                    
                    # decode and check for failure (F = 0 if fails)
                    F = decode(model,state,synds,decoder_gates=gates)

                    step += 1 
                    if step%100000 == 0 print(" $(step/maxsteps) ") end # keep track of progress 
                end 

                if step == maxsteps 
                    println("exceeded maximum steps!")
                end 

                thistrel += step * depth / samps # sample-averaged relaxation time (in units of actual circuit depth, not floquet periods)
                data["trels"][pind,samp] = step * depth 
            end 
            println("thistrel = $thistrel")
            data["avg_trels"][pind] = thistrel
        end 


    ### time evolution of logical fidelity ### 
    elseif mode == "Ft" 
        println("calculating F(t) for $samps samples (with $periods periods each)...")

        data["Fs"] = zeros(nps,periods)

        for (pind,p) in enumerate(ps)
            println("p = $p")
            for samp in 1:samps 
                state = copy(codeword); synds = copy(codesynds)
                for per in 1:periods
                    state, synds = master_gate_applier(model,state,synds,gates,1,[p],bias,gadgetnoise_dict,squarenoise,measurementnoise,false)

                    # decode 
                    data["Fs"][pind,per] += decode(model,state,synds,decoder_gates=gates) / samps 
                end  
            end 
            if l < 3 println("⟨F(t)⟩ = $(data["Fs"][pind,:])") end 
        end 
    
    ### offline decoding error rate ### 
    elseif mode == "offline" 
        println("doing offline decoding for $samps samples (with $periods periods each)...")

        data["offline_error_rate"] = zeros(nps)

        for (pind,p) in enumerate(ps)
            println("p = $p")
            for samp in 1:samps 
                state,synds = noisy_state(p,L,model == "tc" ? 2 : 1) # generate a random state subject to iid errors of strength p
        
                state, synds = master_gate_applier(model,state,synds,gates,1,[0],0,gadgetnoise_dict,squarenoise,false,false)

                # decode 
                data["offline_error_rate"][pind,per] += detect_logical_failure(state) ? 1/samps : 0 

            end 
            println("⟨offline error rate⟩ = $(data["offline_error_rate"][pind])")
        end 

    ### measure statistical properties of the noneq steady state ### 
    elseif mode == "stats"
        @assert model ∈ ["rep" "rep_5bit" "k2"] # only supports 1d codes for now 
        println("getting statistics about long time steady states...")

        # define quantities to measure (some are redundant; are so for ease of plotting): 
        data["average_history"] = zeros(nps,depth,L) # the (signed) averaged spins at each spacetime point in the floquet period 
        data["Ms"] = zeros(nps) # steady-state value of magnetization averaged over *all* times 
        data["|M|s"] = zeros(nps) # ... of absolute value of magnetization 
        data["floquet_Ms"] = zeros(nps) # magnetization averaged over only floquet times  
        data["floquet_|M|s"] = zeros(nps) 
        data["chis"] = zeros(nps) # standard deviation of |M| at all time steps  
        data["floquet_chis"] = zeros(nps) # ... at all floquet periods  
        data["dMdps"] = zeros(nps-1)
        data["binds"] = zeros(nps) # binder cumulants of |M| at all time steps  
        data["floquet_binds"] = zeros(nps) # ... at all floquet periods 
        data["Ps"] = zeros(nps) # codeword density 
        data["floquet_Ps"] = zeros(nps) # ... just at floquet steps 
        data["chiPs"] = zeros(nps) # variance of codeword density 
        data["floquet_chiPs"] = zeros(nps) # ... just at floquet steps 
        data["bindPs"] = zeros(nps) # binder cumulants of P at all time steps  
        data["floquet_bindPs"] = zeros(nps) # ... at all floquet periods 

        data["dw_stats"] = zeros(nps,L) # statistical distribution of domain wall sizes (at floquet times only for now...)
        data["floquet_hamming_statistics"] = zeros(nps,periods)
        data["Es"] = zeros(nps) # average steady state energies at all times  
        data["floquet_Es"] = zeros(nps) # average steady state energies at floquet times 
        data["dEdps"] = zeros(nps-1)
        data["Cs"] = zeros(nps) # susceptibility of energy 
        data["floquet_Cs"] = zeros(nps) # susceptibility measured only at floquet time steps 

        if calc_corrs # 
            println("initializing correlation data...")
            # saving all nps * depth * L * L and nps * depth * depth * L of correlation data creates files of several tens of GB when l = 4, n = 3; hence will just save both an averaged and a particular space / time point 
            data["corrs"] = zeros(nps,2,L,L) # if second two indices are i,j, then stores C_{i,j}. equal-time correlation functions at all pairs of spatial points. second dimension: either averaged over all times [1] or just at floquet times [2]
            data["tcorrs"] = zeros(nps,2,depth,depth) # if second two indices are t1, t2, then stores C_{t1,t1+t2-1}; different format from spatial correlations since correlations not periodic over one circuit depth. store correlations in time up to temporal distance of depth (one floquet period). second dimension: either averaged over all spatial locations [1] or just at the first site [2]. sec
            println("tcorrs has size: $(Base.summarysize(data["tcorrs"]/ 1_073_741_824)) GB") 
        end 

        moving_history = calc_corrs ? zeros(2depth,L) : [] # stores moving history of the spins, to be used in calculating temporal correlation functions (equal-time correlations just calculated using a single history). 
        
        state = rand(Bool,L) # this state will be recyled in between changes in p to help with equilibrating. start with the highest value of p, so begin thermalization from a disordered state 
        for (pind,p) in enumerate(ps)
            println("\np = $p")

            # things to keep track of the history (full or floquet) of the global magnetization 
            Ehist = zeros(periods * depth) # history of energy for all time steps (number of violated stabilizers)
            floquet_Ehist = zeros(periods)

            Mhist = zeros(periods * depth) # history of M for ALL time steps (probably better to do a moving average or something for space reasons...)
            absMhist = zeros(periods * depth) # ... of |M|
            floquet_Mhist = zeros(periods) # history of M for all floquet time steps 
            floquet_absMhist = zeros(periods) # ... of |M|
            Phist = zeros(periods * depth) # history of codeword densities 
            floquet_Phist = zeros(periods) # ... at floquet steps 

            # state = copy(codeword) 
            println("thermalizing...")
            for per in 1:thermal_periods 
                state, _ = master_gate_applier(model,state,[],gates,1,[p],bias,model == "tc" ? gadgetnoise_dict : gadgetnoise,squarenoise,measurementnoise,false) # currently only uniform random gadget failures supported for non-TC codes 
            end 

            println("collecting data...")
            num_corr_measurements = 0 # will divide by this at the end to do averaging 
            for per in 1:periods 
                if per%10000 == 0 print(" $(per/periods) ") end # progress report 

                this_history, _ = master_gate_applier(model,state,[],gates,1,[p],bias,model == "tc" ? gadgetnoise_dict : gadgetnoise,squarenoise,measurementnoise,true) 
                
                state .= this_history[end,:]
                
                data["average_history"][pind,:,:] .+= this_history ./ periods  
                
                if rollstate state = roll(state,rand(1:L)) end 
                
                totals_state = sum(state)
                data["floquet_hamming_statistics"][pind,per] = totals_state 

                # get magnetization 
                floquet_absMhist[per] = abs(2*totals_state/L - 1)
                floquet_Mhist[per] = 2*totals_state/L - 1
                totals_hist = (mean(2 .* this_history .- 1,dims=2)) # spatial average of magnetization at each time step 
                Mhist[(per-1)*depth+1:per*depth] .= totals_hist # magnetization *density*!
                absMhist[(per-1)*depth+1:per*depth] .= abs.(totals_hist) # density of |magnetization|

                # codeword projectors 
                Phist[(per-1)*depth+1:per*depth] .= mean_codeword_projectors(this_history,model)
                floquet_Phist[per] = mean_codeword_projectors(state,model)

                # energies 
                Ehist[(per-1)*depth+1:per*depth] .= energy(this_history,model) # energy density 
                floquet_Ehist[per] = energy(state,model)

                if calc_dw_stats && per%measure_corrs_multiple == 0 
                    data["dw_stats"][pind,:] .+= get_domainwall_stats(this_history[end,:]) / depth 
                end 

                # calculate correlation functions if desired 
                if calc_corrs 
                    if per%measure_corrs_multiple == measure_corrs_multiple-1 # prep the time correlation function part 
                        moving_history[1:depth,:] .= this_history[1:end-1,:] 
                    end 
                        
                    if per%measure_corrs_multiple == 0 
                        num_corr_measurements += 1 
                                                    
                        # spatial correlations 
                        for dt in 1:depth 
                            these_corrs = get_corrs(this_history[dt,:])
                            if dt == 1 
                                data["corrs"][pind,2,:,:] .+= these_corrs # just at floquet times 
                            end 
                            data["corrs"][pind,1,:,:] .+= these_corrs / depth # average over all times 
                        end 

                        moving_history[depth+1:end,:] .= this_history[1:end-1,:] # fill in the second half of the temporal correlation data 
                        for t0 in 1:depth, dt in 0:depth-1 
                            data["tcorrs"][pind,2,t0,dt+1] += moving_history[t0,1] * moving_history[t0+dt,1]  # just a single point 
                            for i in 1:L 
                                data["tcorrs"][pind,1,t0,dt+1] += moving_history[t0,i] * moving_history[t0+dt,i] / L # averaged over all space 
                            end 
                        end 
                    end 
                end 
            end # end loop over periods 

            # now collect the statistics of magnetization / codeword density 
            data["Ms"][pind] = mean(Mhist)
            data["|M|s"][pind] = mean(absMhist)
            data["floquet_Ms"][pind] = mean(floquet_Mhist)
            data["floquet_|M|s"][pind] = mean(floquet_absMhist) 
            data["Es"][pind] = mean(Ehist)
            data["Cs"][pind] = (mean(Ehist.^2) - mean(Ehist)^2) * L # normalization so that scales as L^0 away from critical point (by CLT)
            data["floquet_Es"][pind] = mean(floquet_Ehist)
            data["floquet_Cs"][pind] = (mean(floquet_Ehist.^2) - mean(floquet_Ehist)^2) * L

            if bias != 0 # calculate χ and B with M 
                Mhist_mean = mean(Mhist)
                floquet_Mhist_mean = mean(floquet_Mhist)

                data["binds"][pind] = mean((Mhist .- Mhist_mean).^4) / (mean((Mhist .- Mhist_mean).^2)^2)
                data["floquet_binds"][pind] = mean((floquet_Mhist .- floquet_Mhist_mean).^4) / (mean((floquet_Mhist .- floquet_Mhist_mean).^2)^2)
                data["chis"][pind] = (mean(Mhist.^2) - Mhist_mean^2) * L 
                data["floquet_chis"][pind] = (mean(floquet_Mhist.^2) - floquet_Mhist_mean^2) * L 
    
            else # ...calculate with |M|
                absMhist_mean = mean(absMhist)
                floquet_absMhist_mean = mean(floquet_absMhist)

                data["binds"][pind] = mean((absMhist).^4) / (mean((absMhist).^2)^2)
                data["floquet_binds"][pind] = mean((floquet_absMhist).^4) / (mean((floquet_absMhist).^2)^2)
                data["chis"][pind] = (mean(absMhist.^2) - absMhist_mean^2) * L 
                data["floquet_chis"][pind] = (mean(floquet_absMhist.^2) - floquet_absMhist_mean^2) * L 
            end 

            Phist_mean = mean(absMhist)
            floquet_Phist_mean = mean(floquet_absMhist)

            data["bindPs"][pind] = mean((Phist).^4) / (mean((Phist).^2)^2)
            data["floquet_bindPs"][pind] = mean((floquet_Phist).^4) / (mean((floquet_Phist).^2)^2)
            data["chiPs"][pind] = (mean(Phist.^2) - Phist_mean^2) * L 
            data["floquet_chiPs"][pind] = (mean(floquet_Phist.^2) - floquet_Phist_mean^2) * L 


            # take off the disconnected parts of the correlators 
            if calc_corrs 
                data["corrs"] ./= num_corr_measurements
                data["tcorrs"] ./= num_corr_measurements 
                for i in 1:L, j in 1:i 
                    # average over all times: 
                    discon_part = sum(data["average_history"][pind,k,i] * data["average_history"][pind,k,j] for k in 1:depth) / depth # average the disconnected part over times 
                    data["corrs"][pind,1,i,j] -= discon_part 
                    if j < i 
                        data["corrs"][pind,1,j,i] -= discon_part 
                    end 
                    # just at floquet times 
                    discon_part = data["average_history"][pind,1,i] * data["average_history"][pind,1,j] 
                    data["corrs"][pind,2,i,j] -= discon_part 
                    if j < i 
                        data["corrs"][pind,2,j,i] -= discon_part 
                    end 
                end 
                for t0 in 1:depth 
                    for dt in 0:depth-1 
                        data["tcorrs"][pind,1,t0,dt+1] -= sum(data["average_history"][pind,t0,i] * data["average_history"][pind,(t0+dt+depth-1)%depth+1,i] for i in 1:L)/L # average over space 
                        data["tcorrs"][pind,2,t0,dt+1] -= data["average_history"][pind,t0,1] * data["average_history"][pind,(t0+dt+depth-1)%depth+1,1] # just at a single spatial point 
                    end 
                end 
            end  
        end 

        # also compute heat capacity / magnetization by looking at derivative wrt noise: 
        data["halfps"] = [(ps[i]+ps[i-1])/2 for i in 2:nps]
        data["dEdps"] = [(data["Es"][i]-data["Es"][i-1])/(ps[i]-ps[i-1]) for i in 2:nps]
        data["d|M|dps"] = [(data["|M|s"][i]-data["|M|s"][i-1])/(ps[i]-ps[i-1]) for i in 2:nps]
    
    ### run a single spacetime history (for visualization purposes) or search over histories satisfiying some criterion ### 
    elseif mode == "hist" 
        println("running spacetime histories...")
        initial_state = copy(codeword) # clean initial state 
        periods = 1; err_hist = []; synd_hist = []
        if model ∈ ["rep" "rep_5bit" "k2" "twod_rep"] 
            failed = false 
            p = 0.05; bias = .5
            maxtries = 1000000; thistry = 0 
            maxtries = 1 
            if model == "twod_rep" # large square domain of errors -- for visualization purposes 
                initial_state = falses(Lx,Ly)
                Lthird = round(Int,L/3)
                initial_state[Lthird+1:2Lthird,Lthird+1:2Lthird] .= 1 
            end 
            while ~failed && thistry < maxtries 
                history = master_gate_applier(model,copy(initial_state),[],gates,1,[p],bias,false,false,measurementnoise,true)
                if sum(history[end,:,:]) != 0
                    failed = true
                    println("failed!")
                    println("⟨final spins⟩ = ",mean(history[end,:,:]))
                end 
                thistry += 1 
            end 
            if thistry ≥ maxtries 
                println("no logical failures detected")
            end 

        elseif model == "tc"
        
            function failure_test(final_synds) 
                """
                function determining various failure criteria for different gadgets; searches with this not used in proof of fault tolerance
                """
                failure = false 

                justanyons = ~false # fail if any anyons are created 
                diagonallink = false # fail if a diagonal link is created 
                boxable = ~false # fail if error is not boxable in some particular way 
                allowdiagonallinkbox = true # if false, fails if error cannot be made straight. thus a single-boxable diagonal error will fail. if true, failure will occur only if the error is not 1-boxable 

                if justanyons 
                    return sum(final_synds) > 0  
                end 
                if diagonallink 
                    if sum(final_synds) > 0 
                        anyon_locs = findall(x->x==true,final_synds) # test for diagonal stuff 
                        firstrow = false; secondrow = false 
                        firstcol = false; secondcol = false
                        for aloc in anyon_locs 
                            if aloc[1] == 4 firstrow = true end 
                            if aloc[1] == 7 secondrow = true end 
                            if aloc[2] == 4 firstcol = true end 
                            if aloc[2] == 7 secondcol = true end 
                        end 
                        if firstcol && thirdcol && sum(final_synds) == 2 failure = true end 
                    end 
                end 
                if boxable 
                    if ~single_link_test(final_synds,allowdiagonallinkbox)
                        td_synds = copy(final_synds)
                        for j in 0:2 
                            if td_synds[4,1+3*j] 
                                td_synds[1,1+3*j] ⊻= true; td_synds[7,1+3*j] ⊻= true 
                            end 
                        end 
                        if ~single_link_test(td_synds,allowdiagonallinkbox)
                            println("failed! ")
                            println("final synds    = ",final_synds)
                            println("T(final synds) = ",td_synds)
                            failure = true 
                        end 
                    end 
                end 
                return failure 
            end 
            
            ### options for performing search ### 
            corners = ~true # if true, includes diagonal links in the error model 
            exhaustive_search = true; if exhaustive_search periods = 2 end # if exhaustive search, do a noisy EC and then a clean EC 
            save_history = ~true  
            maxtries = 1000000; thistry = 0 
            # maxtries = 1 # if maxtries = 1, does custom error configurations 
            
            transient = true # if true, search is over errors at all spacetime locations; if false, errors occur only in the initial state 
            look_for_corner_error = false # if true, fails if an l-scale corner error occurs; if false, fails if either straight link or corner errors occur 

            error_weight = 2
            error_weights = zeros(Int,3,3) # weights of errors in different 1-cells for a 9x9 grid (used for testing failure properties of level 1 gadgets)

            R0_depth = 24 
            init_errs = falses(Lx,Ly,2); init_synds = falses(Lx,Ly)
            noise_hist = falses(periods*depth+1,Lx,Ly,4)

            failed = false 
            if exhaustive_search  # do an explicit enumeration over all possible error configurations; errors occur in the first half of the gadget (assumed to be EC)
                
                @assert periods == 2

                tsearch_max = depth 
                Lxp = Lx; Lyp = Ly; op = 2 # just to cut down on the search space a bit if only 1 error 
                if error_weight == 1 Lxp = floor(Int,Lx/2); Lyp = floor(Int,Ly/2); op = 1 end 

                for tsearch in 1:tsearch_max 
                    println("progress: $(tsearch/tsearch_max)")

                    valid_noise_indices = [(tsearch, j, k, l) for j in 1:Lxp, k in 1:Lyp, l in 1:(corners ? 4 : op)] # 
                            
                    all_error_combinations = combinations(valid_noise_indices, error_weight) # automatically gets rid of permutaiton-equivalent indices 
                    errors_to_search = length(all_error_combinations)
                    println("beginning the long march for error weight $error_weight")
                    println("number of errors to search: $errors_to_search")
                    progress_report = 0 
                    for true_indices in all_error_combinations
                        # if progress_report % 10 == 0 print(" $(progress_report / errors_to_search)") end 
                        for true_index in true_indices # fix the current noise 
                            noise_hist[true_index...] ⊻= true 
                        end 
                
                        # diagnose failure 
                        errs, synds = tc_gate_applier!(copy(init_synds),copy(init_errs),gates,periods,noise_hist,gadgetnoise_dict,false,false,false)
                        if look_for_corner_error
                            failed = detect_both_logical_failures(errs)  
                        else 
                            failed = detect_logical_failure(errs)  
                        end 

                        if failed 
                            println("failure with weight $error_weight at locations: $true_indices")
                            break 
                        end 

                        for true_index in true_indices # reset the noise 
                            noise_hist[true_index...] ⊻= true  
                        end 
                        progress_report += 1 
                    end
                end 
                if ~failed println("search complete; no failures found. ") end 

            else # random sampling of error configurations or evaluation on a specific error configuration 
                custom_errors = maxtries == 1 

                while ~failed && thistry < maxtries 
                    if thistry % 10 == 0 print(" $(thistry/maxtries) ") end 
                    # last dimension of noise_hist is 1 / 2 / 4 depending on plaquette / straight link / straight and diag link noise 
                    noise_hist = falses(periods*depth+1,Lx,Ly,4)

                    if ~custom_errors # do a random sampling of errors:

                        if ~gadgetnoise
                            ts = Int[1 for i in 1:error_weight]; ts[1] = rand(1:depth)
                            ts[1] = rand(1:depth) 
                            for repind in 1:error_weight
                                if repind > 1 ts[repind] = min(max(1,ts[1] + rand(-5:5)),depth) end # search only over errors that occur close-ish together in time 
                                noise_hist[transient ? ts[repind] : 1,rand(1:Lx),rand(1:Ly),rand(1:(corners ? 4 : 2))] ⊻= true 
                            end 
                        else 
                            noise_hist = [1 for i in 1:error_weight]; noise_hist[1] = rand(1:ngates) # label of the gadgets that fail 
                            for repind in 1:error_weight
                                noise_hist[repind] = min(max(1,noise_hist[1] + rand(-2*Lx^2:2*Lx^2)),ngates) # again require that different errors not be too far apart in time 
                            end 
                        end 


                    else # custom errors 
                        save_history = true 
                        println("custom error occuring!")

                        ### following code block needs to be manually adjusted to yield the desired spacetime noise history ### 

                        noise_inds = []
                        # random error 
                        # error_weight = 40
                        # for i in 1:error_weight 
                        #     push!(noise_inds,CartesianIndex(1,rand(1:Lx),rand(1:Ly),rand(1:4)))
                        # end 

                        # level-1 corner error 
                        # for i in 1:3 
                        #     push!(noise_inds,CartesianIndex(1,3+i,4,1))
                        #     push!(noise_inds,CartesianIndex(1,4,3+i,2))
                        # end 
                        
                        # push!(noise_inds,CartesianIndex(1,1,1,2))
                        push!(noise_inds,CartesianIndex(1,2,2,4))
                        
                        for nind in noise_inds 
                            noise_hist[nind] ⊻= true 
                        end 
                        println("∑noise = ",sum(noise_hist))
                    end 

                    if save_history 
                        # save the whole history (slower)
                        err_hist, synd_hist = tc_gate_applier!(copy(init_synds),copy(init_errs),gates,periods,noise_hist,gadgetnoise_dict,false,false,true)

                        if look_for_corner_error
                            failed = detect_both_logical_failures(err_hist[end,:,:,:]) # look for logical failures along both directions  
                        else 
                            failed = detect_logical_failure(err_hist[end,:,:,:]) # look for logical failure aong either direction  
                        end 
                        # failed = failure_test(synd_hist[end,:,:]) # look for anyons not being cleaned up in a particular way 

                    else 
                        # just keep track of the final configuration (faster)
                        errs, synds = tc_gate_applier!(copy(init_synds),copy(init_errs),gates,periods,noise_hist,gadgetnoise_dict,false,false,false)
                        if look_for_corner_error
                            failed = detect_both_logical_failures(errs)  
                        else 
                            failed = detect_logical_failure(errs)  
                        end 
                    end                     

                    thistry += 1 
                    if failed || maxtries == 1  
                        if ~gadgetnoise
                            println("spacetime coordinates of noise leading to the failure:")
                            for index in findall(x->x,noise_hist)
                                println(index)
                            end 
                        elseif maxtries > 1 
                            println("failed gates = $noise_hist")
                        end 
                        data["noise_hist"] = noise_hist; data["err_hist"] = err_hist; data["synd_hist"] = synd_hist
                    end 
                end 
                # if thistry ≥ maxtries 
                #     println("logcal failure in final state: ",detect_logical_failure(err_hist[end,:,:,:]))
                # end 
            end             
        end 
    else
        println("unsupported mode: $mode")
    end 

    ### save stuff to file ### 

    # information about the simulation awkwardly compiled into a giant dictionary... (would ofc be prettier to collect stuff with an argument parser)
    params = Dict{String, Any}()
    params["mode"] = mode; params["model"] = model; params["gate"] = gate; params["L"] = L; params["l"] = l; params["squarenoise"] = squarenoise; params["bias"] = bias; params["gadgetnoise"] = gadgetnoise; params["measurementnoise"] = measurementnoise; params["nps"] = nps; params["logdist"] = logdist; params["calc_corrs"] = calc_corrs; params["measure_corrs_multiple"] = measure_corrs_multiple; params["calc_dw_stats"] = calc_dw_stats; params["rollstate"] = rollstate; params["ps"] = ps; params["depth"] = depth; params["samps"] = samps; params["pc"] = pc; params["thermal_periods"] = thermal_periods; params["periods"] = periods

    println("saving data to file: $fout")
    f = jldopen(fout,"w")
    for (key,object) in data ∪ params 
        write(f,key,object)
    end 
    close(f)
    
end 

main()

