using JLD2

include("gate_builder.jl")
include("gate_applier.jl")
include("helper_functions.jl")

"""
this code runs monte carlo simulations of various different types of self-correcting circuit dynamics.
can compute: 
    * time dependence of logical fidelity 
    * relaxation times 
    * statistical properties and correlation functions in long-time steady states 

arguments and options given in the first section of main()
"""

function roll(arr,n)
    """
    a julia version of numpy's roll 
    """
    len = length(arr)
    n = mod(n,len) 
    if n == 0
        return arr
    else
        return vcat(arr[end-n+1:end], arr[1:end-n])
    end
end

function parameter_repository(model,mode,l,nps,logdist,bias,gadgetnoise,measurementnoise)
    """ 
    this function is a place to store numerical information about the locations of thresholds in different types of models + conventions on number of samples to use + etc. 

    model: name of code ∈ ["rep" "k2" "tc" "twod_rep"]
    l: log_n of system size 
    nps: desired number of noise steps 
    logdist: returns linearly-spaced noise values if false; log-spaced values if true 
    bias: bias of iid noise (not relevant for TC)
    gadgetnoise: if true, gadget-based noise 
    measurementnoise: if true, noise occurs only in syndrome measurements 

    returns: (vector of noise values, pc, number of samples)
    """
    pmin = 0; pmax = 0 
    ps = []; pc = 0 
    samps = 1 # O(1) if mode == stats, large for Ft, slightly smaller possible for trel 
    periods = 1 # number of applications of one floquet period for the total time evolution (not needed if mode == "trel")
    thermal_periods = 0 # number of periods to thermalize for before measuring properties of the steady state 

    ####### number of periods #######  
    if mode == "Ft" # this number does not really matter -- just sets how close to 1 the value of F at the crossing point is (closer for smaller periods)
        periods = 50 
        if model == "rep_5bit" periods = 5 end 
        if model == "tc" periods = 10 end 
    end 
    if mode == "stats"
        thermal_periods = 10000 # pooperiods
        if l == 0 periods = 5000000 end 
        if l == 1 periods = 4000000 end 
        if l == 2 periods = 9000000 end 
        if l == 3 periods = 800000 end 
        if l == 4 periods = 90000 end 

        if model ∈ ["k2" "rep_5bit"] # both have block size 5 -- use less periods 
            if l == 1 periods = 200000; thermal_periods = 20000 # 500k periods is very easily enough to get good statistics 
            elseif l == 2 periods = 25000; thermal_periods = 14000 
            elseif l == 3 periods = 5000; thermal_periods = 1200 
            elseif l == 4 periods = 600; thermal_periods = 200
            else periods = 200; thermal_periods = 20 end 
        end 
    end 
    if mode == "hist" periods = 1 end 

    ####### noise strengths and samples: 1d codes ####### poops poosamps 

    ### 3-bit rep code ### 
    if model == "rep"
        pc = bias == 0 ? .0144 : .0063 
        if gadgetnoise pc = .008 end 
        # Fidelity: 
        if mode == "Ft"
            if bias == 0 pmin = .011; pmax = .018  
            else pmin = .005; pmax = .008 end 
            if l == 3 pmin = .012; pmax = .017 end 
            if l == 4 pmin = .012; pmax = .0155 end 
            if bias == 1
                if l == 3 pmin = .0055; pmax = .00725 end 
                if l == 4 pmin = .0055; pmax = .0069 end    
            end  

            if gadgetnoise pmin = 5e-3; pmax = 11e-3 end # note: gadgetnoise currently randomizes gadget outputs and hence currently does not support a bias 
            println("l = $l, samps = $samps")
            if l == 1 samps = 100000 
            elseif l == 2 samps = 25000
            elseif l == 3 samps = 7000 #5000
            elseif l == 4 samps = 5000
            else samps = 50
            end 

        elseif mode == "trel"
            if l == 1 pmin = .006; pmax = .02; samps = 1200  
            elseif l == 2 pmin = .009; pmax = .019; samps = 250 
            elseif l == 3 pmin = .011; pmax = .0155; samps = 90 # this takes quite a while... 
            elseif l == 4 pmin = .0118; pmax = .0144; samps = 40 
            elseif l == 5 pmin = .0058; samps = 25 
            else pmin = .005; samps = 15 
            end  

            if gadgetnoise # crossing of relaxation time points also give same p_c as trel 
                if l == 1 pmin = .0005; pmax = .012; samps = 500 
                elseif l == 2 pmin = .003; pmax = .012; samps = 100 
                elseif l == 3 pmin = .0055; pmax = .012; samps = 25
                else pmin = .006; pmax = .01; samps = 50 
                end 
            end 
        elseif mode == "hist" pmin = pc / 5; pmax = pmin 
        elseif mode == "stats" pmin = .7 * pc; pmax = 1.25 * pc 
        end 

    ### 5-bit rep code ### 
    elseif model == "rep_5bit"
        pc = bias == 0 ? .0142 : .0065 # basically the same as the 3-bit repetition code?! 

        if mode == "Ft"
            if l == 1 pmin = .005; pmax = .02; samps = 1500
            elseif l == 2 pmin = .005; pmax = .02; samps = 400
            elseif l == 3 pmin = .0128; pmax = .019
            elseif l == 4 pmin = .0131; pmax = .0158
            end 
            pmin = .005; pmax = .02
        elseif mode == "trel"
            if l == 1 pmin = .005; pmax = .045; samps = 1
            elseif l == 2 pmin = .0105; pmax = .026
            elseif l == 3 pmin = .0128; pmax = .019
            elseif l == 4 pmin = .0131; pmax = .0158
            end 
        end 

    ### k = 2 code ### 
    elseif model == "k2"
        pc = 7.5e-4 

        if mode == "Ft"
            nothing 
        elseif mode == "trel"
            nothing 
        end 
    
    
    ####### noise strengths and samples: 2d codes #######

    ### 2d version of tsirelson (for repetition codes only; just does 1d Tsirelson along two directions) ### 
    elseif model == "twod_rep" 
        pc = bias == 0 ? .01 : .01 
        if mode == "Ft"
            nothing 
        elseif mode == "trel"
            nothing 
        end 

    ### toric code ### 
    elseif model == "tc"
        # fidelity 
        if mode == "Ft"
            pc = .001 
            pmin = 1e-4; pmax = 5e-3; samps = 50  
            pmin = 1e-2; pmax = 1e-1; samps = 100 # for prim R0
        # relaxation time 
        elseif mode == "trel"
            if l == 0 pmin = 1/exp(10); pmax = 1/exp(7.5); samps = 4000 end # trel bigger than around one period by .01ish
            if l == 1 pmin = 1/exp(11.1); pmax = 7e-4; samps = 1000 end # around 7e-4, may flatten out around 1e-4
            if l == 2 pmin = 1/exp(10.2); pmax = 1/exp(9); samps = 100 end 
            if l == 3 pmin = 1/exp(5); pmax = 1/exp(4); samps = 10 end 

            if measurementnoise
                if l == 0 pmin = 1/exp(10.9); pmax = 1/exp(7.25); samps = 1500 end 
                if l == 1 pmin = 1/exp(10.2); pmax = 1/exp(7.75); samps = 800 end 
                if l == 2 pmin = 1/exp(12.5); pmax = 1/exp(10.25); samps = 150 end # trel > T_2 at about 1/exp(10)... 
            end     

            if gadgetnoise # for "fat" gadget supports as defined in proof of nilpotence 
                if l == 0 pmin = 1/exp(13); pmax = 1/exp(8); samps = 1000 end 
                if l == 1 pmin = 1/exp(13); pmax = 1/exp(8); samps = 600 end 
                if l == 2 pmin = 1/exp(13.8); pmax = 1/exp(11.8); samps = 150 end # trel > T_2 at around 1/exp(11)
            end 

        end 
    end 

    ps = []
    if logdist 
        ps = [10^(x) for x in range(log10(pmax),log10(pmin);length=nps)] # ps from largest to smallest so that states can be reused when going to smaller p (helps with equilibration) when measuring steady-state properties 
    else 
        ps = [x for x in range(pmax,pmin;length=nps)] 
    end   
    reverse!(ps) # better to start from smallest p for exploratory purposes (when measuring trel) to understand how fast things will run

    return ps, pc, samps, thermal_periods, periods 
end 

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

function main() 
    ##### arguments and options #####
    
    """
    basic structural options 
    mode: what the code calculates. one of the following: 
        * "Ft": time evolution of logical fidelity (starting in 0^L codeword by default)
        * "trel": relaxation time of logical information 
        * "stats": collects statistics about the long time (t >> trel) steady state. used when computing critical properties. 
        * "hist": just runs a single sample of the dynamics; used for visualization / testing specific noise configurations 
    model ∈ ["rep" "rep_5bit" "k2" "tc" "twod_rep"]
    gate ∈ ["Z" "R" "decoder" "Y" "id" ...]
    l: level, equals log_n(system size)-1 for EC gates and log_n(system_size) for other gates (in our current (somewhat unfortunate) conventions)
    """
    mode = "hist" # 
    mode = "trel" 
    model = "tc"
    gate = "R"
    l = 2 
    n = model ∈ ["rep" "tc" "twod_rep"] ? 3 : 5
    Lx, Ly = get_gate_size(l,model,gate); L = Lx  
    
    """
    noise options 
    squarenoise: if true, applies errors by randomly applying plaquette errors (rather than link errors). 
    measurementnoise: if true, noise comes only from failures in the syndrome measurements
    bias: bias of iid noise (not used for TC)
    gadgetnoise: if true, gadgets fail (all bits in the output of a failed gadget are randomized). if false, "wires" between gadgets fail.
    nps: number of noise values to sample 
    logdist: if true, spaces out noise values evenly on a log scale. if false, even on a linear scale. 
    """
    squarenoise = false 
    measurementnoise = ~false 
    bias = 0.  
    gadgetnoise = false 
    nps = mode == "trel" ? (l > 1 ? 7 : 9) : (mode == "Ft" ? 10 : (mode == "hist" ? 1 : (mode == "stats" ? 7 : 1)))  
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
    adj: distinguishing information about input file 
    out_adj: ... about output file 
    """
    adj = "" 
    out_adj = "" 

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

    ps,pc,samps,thermal_periods,periods = parameter_repository(model,mode,l,nps,logdist,bias,gadgetnoise,measurementnoise)
    fout = "data/$(model)$(l)_$(mode)_$(gate)_"*(mode == "trel" ? "$(periods)pers" : "")*"$(out_adj).h5"
    if mode == "trel"
        fout = "data/$(model)_l$(l)_$(mode)_$(gate)$(out_adj).h5"
    end 
    if mode == "stats"
        fout = "data/$(model)_l$(l)_$(mode)_$(gate)$(out_adj).h5"
    end 
    println("will save data as: $fout")
    
    # load / generate gate data 
    gates = []; depth = 24
    # check to see if gate data already exists: 
    fin = "data/$(model)_l$(l)_gates_$(gate)$(adj).jld2" 
    if isfile(fin)
        println("loading data from file: $fin")
        f = jldopen(fin,"r")
        gates = read(f,"gates"); depth = read(f,"depth")
        println("depth = $depth")
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
    ngates = size(gates)[1]
    
    # all the things that we may have reason to measure 
    datakeys = ["trels" "avg_trels" "Fs" "|M|s" "Ms" "floquet_|M|s" "floqet_Ms" "chis" "floquet_chis" "binds" "floquet_binds" "average_history" "dw_stats" "corrs" "tcorrs" "noise_hist" "err_hist" "synd_hist" "floquet_hamming_statistics"]
    data = Dict{String, Any}()
    for key in datakeys data[key] = zeros(1) end 

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
        data["binds"] = zeros(nps) # binder cumulants of |M| at all time steps  
        data["floquet_binds"] = zeros(nps) # ... at all floquet periods 
        data["dw_stats"] = zeros(nps,L) # statistical distribution of domain wall sizes (at floquet times only for now...)
        data["floquet_hamming_statistics"] = zeros(nps,periods)
        
        data["corrs"] = zeros(nps,depth,L,L) # equal-time correlation functions at all floquet-distinct times 
        data["tcorrs"] = zeros(nps,L,depth,depth) # store correlations in time up to temporal distance of depth (one floquet period)

        moving_history = calc_corrs ? zeros(2depth,L) : [] # stores moving history of the spins, to be used in calculating temporal correlation functions (equal-time correlations just calculated using a single history). 
        
        state = rand(Bool,L) # this state will be recyled in between changes in p to help with equilibrating. start with the highest value of p, so begin thermalization from a disordered state? 
        for (pind,p) in enumerate(ps)
            println("\np = $p")

            # things to keep track of the history (full or floquet) of the global magnetization 
            Mhist = zeros(periods * depth) # history of M for ALL time steps (probably better to do a moving average or something for space reasons...)
            absMhist = zeros(periods * depth) # ... of |M|
            floquet_Mhist = zeros(periods) # history of M for all floquet time steps 
            floquet_absMhist = zeros(periods) # ... of |M|

            # state = copy(codeword) 
            println("thermalizing...")
            for per in 1:thermal_periods 
                state, _ = master_gate_applier(model,state,[],gates,1,[p],bias,gadgetnoise_dict,squarenoise,measurementnoise,false)
            end 

            println("collecting data...")
            num_corr_measurements = 0 # will divide by this at the end to do averaging 
            for per in 1:periods 
                if per%10000 == 0 print(" $(per/periods) ") end # progress report 

                this_history, _ = master_gate_applier(model,state,[],gates,1,[p],bias,gadgetnoise_dict,squarenoise,measurementnoise,true) 
                
                state .= this_history[end,:]
                
                data["average_history"][pind,:,:] .+= this_history[1:end-1,:] / periods  # don't include the temporally-last configuration since it will be included as the starting configuration in the next period 
                
                if rollstate state = roll(state,rand(1:L)) end 
                
                totals_state = sum(state)
                data["floquet_hamming_statistics"][pind,per] = totals_state 

                floquet_absMhist[per] = abs(2*totals_state - 1)
                floquet_Mhist[per] = 2*totals_state - 1
                
                totals_hist = (mean(2 .* this_history[1:end-1,:] .- 1,dims=2))
                Mhist[(per-1)*depth+1:per*depth] .= totals_hist 
                absMhist[(per-1)*depth+1:per*depth] .= abs.(totals_hist)
                
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
                            data["corrs"][pind,dt,:,:] .+= get_corrs(this_history[dt,:])  
                        end 

                        moving_history[depth+1:end,:] .= this_history[1:end-1,:] # fill in the second half of the temporal correlation data 
                        for t0 in 1:depth for dt in 0:depth-1 for i in 1:L # calculate the temporal correlators
                            data["tcorrs"][pind,i,t0,dt+1] += moving_history[t0,i] * moving_history[t0+dt,i]
                        end end end 
                    end 
                end 
            end  

            # now collect the statistics of magnetization 
            data["Ms"][pind] = mean(Mhist)
            data["|M|s"][pind] = mean(absMhist)
            data["floquet_Ms"][pind] = mean(floquet_Mhist)
            data["floquet_|M|s"][pind] = mean(floquet_absMhist) 

            if bias != 0 # calculate χ and B with M 
                Mhist_mean = mean(Mhist)
                floquet_Mhist_mean = mean(floquet_Mhist)

                data["binds"][pind] = mean((Mhist .- Mhist_mean).^4) / (mean((Mhist .- Mhist_mean).^2)^2)
                data["floquet_binds"][pind] = mean((floquet_Mhist .- floquet_Mhist_mean).^4) / (mean((floquet_Mhist .- floquet_Mhist_mean).^2)^2)
                data["chis"][pind] = (mean(Mhist.^2) - Mhist_mean^2) / L 
                data["floquet_chis"][pind] = (mean(floquet_Mhist.^2) - floquet_Mhist_mean^2) / L 
    
            else # ...calculate with |M|
                Mhist_mean = mean(absMhist)
                floquet_Mhist_mean = mean(floquet_absMhist)

                data["binds"][pind] = mean((absMhist).^4) / (mean((absMhist).^2)^2)
                data["floquet_binds"][pind] = mean((floquet_absMhist).^4) / (mean((floquet_absMhist).^2)^2)
                data["chis"][pind] = (mean(absMhist.^2) - absMhist_mean^2) / L 
                data["floquet_chis"][pind] = (mean(floquet_absMhist.^2) - floquet_absMhist_mean^2) / L 
            end 

            # take off the disconnected parts of the correlators 
            if calc_corrs 
                data["corrs"] ./= num_corr_measurements
                data["tcorrs"] ./= num_corr_measurements 
                for t0 in 1:depth 
                    for i in 1:L 
                        for j in 1:L 
                            data["corrs"][pind,t0,i,j] -= data["average_history"][pind,t0,i] * data["average_history"][pind,t0,j] 
                        end 
                    end
                end 
                for t0 in 1:depth 
                    for dt in 0:depth-1 
                        for i in 1:L 
                            data["tcorrs"][pind,i,t0,dt+1] -= data["average_history"][pind,t0,i] * data["average_history"][pind,(t0+dt+depth-1)%depth+1,i]
                        end 
                    end 
                end 
            end  
        end 
    
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

            # R0s_depth = primitive_R0 ? 1 : 24 
            init_errs = falses(Lx,Ly,2); init_synds = falses(Lx,Ly)
            noise_hist = falses(periods*depth+1,Lx,Ly,4)
            
            error_weights = zeros(Int,3,3) # weights of errors in different 1-cells for a 9x9 grid (used for testing failure properties of level 1 gadgets)

            corners = ~true # if true, includes diagonal links in the error model 
            exhaustive_search = ~true
            save_history = ~true  

            if exhaustive_search  # do an explicit enumeration over all possible n-error configurations (slow) 
                failed = false 
                
                # spatial coordinates included during the search 
                xsearchmin = 1; xsearchmax = 7
                ysearchmin = 2; ysearchmax = 7
                tsearchmin = R0s_depth; tsearchmax = depth-R0s_depth
                
                function long_march() # written as a function to more easily facilitate a break out of multiple for loops 
                    for t1 in tsearchmin:tsearchmax 
                        println("t progress: $(t1 / (depth-R0s_depth))")
                        for x1 in xsearchmin:xsearchmax 
                        for y1 in ysearchmin:ysearchmax for a1 in 1:(corners ? 4 : 2)
                        println("a progress = $a1") 
                            
                        noise_hist[t1,x1,y1,a1] = 1 
                        for t2 in t1:tsearchmax for x2 in xsearchmin:xsearchmax for y2 in ysearchmin:ysearchmax for a2 in 1:(corners ? 4 : 2)
                            if failed 
                                println("failed! \n error locations: \n ($t1,$x1,$y1,$a1) \n ($t2,$x2,$y2,$a2)")
                                return 
                            end 
                            noise_hist[t2,x2,y2,a2] = 1 

                            err_hist, synd_hist = twod_majority_gate_applier(copy(init_synds),copy(init_errs),gates,periods,l,L,noise_hist,true,square_noise)#,vert_part,horz_part,just_central_ms)
                        
                            failed = failure_test(synd_hist[end,:,:]) # look for failure to clean up anyons in a certain way 
                            
                            noise_hist[t2,x2,y2,a2] = 0 # reset the noise 
                        end end end end 
                        noise_hist[t1,x1,y1,a1] = 0 
                    end end end end 
                    println("no failures detected...")
                end 

                long_march()

            else # random error configurations 
                failed = false 
                maxtries = 1000000; thistry = 0 
                maxtries = 1 

                custom_errors = maxtries == 1 
                transient = true # if true, search is over errors at all spacetime locations; if false, errors occur only in the initial state 
                while ~failed && thistry < maxtries 
                    if thistry % 10 == 0 print(" $(thistry/maxtries) ") end 
                    # last dimension of noise_hist is 1 / 2 / 4 depending on plaquette / straight link / straight and diag link noise 
                    noise_hist = falses(periods*depth+1,L,L,4)

                    if ~custom_errors
                        # for level 1 gadgets, the following determine the weights of errors in the 9 different 1-cells of a 9x9 grid 
                        # error_weights[1,1] = 0; error_weights[2,1] = 0; error_weights[3,1] = 0 # top    row 
                        # error_weights[1,2] = 1; error_weights[2,2] = 1; error_weights[3,2] = 0 # middle row 
                        # error_weights[1,3] = 0; error_weights[2,3] = 0; error_weights[3,3] = 0 # bottom row 
                        # # fill in the errors accordingly
                        # for i in 0:2 for j in 0:2 
                        #     for d in 1:error_weights[i+1,j+1]
                        #         direction = rand(1:(corners ? 4 : 2))
                        #         dx = rand(1:3); dy = rand(1:3)
                        #         noise_hist[transient ? rand(1:depth-R0s_depth) : 1,3*i+dx,3*j+dy,direction] = 1 # errors should occur before the ideal 0-lvl decoder at the end 
                        #     end 
                        # end end 

                        # alternatively, just do random errors: 
                        error_weight = 1 # pooerrorweight
                        if ~gadgetnoise
                            ts = Int[1 for i in 1:error_weight]; ts[1] = rand(1:depth)
                            ts[1] = rand(1:depth) #; ts[2] = rand(1:depth) # ts[2] = min(max(1,ts[1] + rand(-5:5)),depth)
                            for repind in 1:error_weight
                                if repind > 1 ts[repind] = min(max(1,ts[1] + rand(-5:5)),depth) end 
                                noise_hist[transient ? ts[repind] : 1,rand(1:Lx),rand(1:Ly),rand(1:(corners ? 4 : 2))] ⊻= true 
                            end 
                        else 
                            noise_hist = [1 for i in 1:error_weight]; noise_hist[1] = rand(1:ngates) # label of the gadgets that fail 
                            for repind in 1:error_weight
                                noise_hist[repind] = min(max(1,noise_hist[1] + rand(-2*Lx^2:2*Lx^2)),ngates)
                            end 
                        end 


                    else # custom errors 
                        save_history = true 
                        println("custom error occuring!")
                        # thing that gives asym ty failure 
                        # noise_synds = Any[(0, 0), (1, 0), (2, 0), (1, 1)]
                        # xerr = 2; yerr = 4; t = 23 

                        # noise_hist = get_noise_hist(xerr,yerr,t,noise_synds,periods*depth+1,L,L)

                        # noise_hist = falses(periods*depth+1,L,L,4)

                        noise_inds = []
                        # random error 
                        # for i in 1:40 
                        #     push!(noise_inds,CartesianIndex(1,rand(1:Lx),rand(1:Ly),rand(1:4)))
                        # end 
                        # level-1 corner error 
                        for i in 1:3 
                            push!(noise_inds,CartesianIndex(1,3+i,4,1))
                            push!(noise_inds,CartesianIndex(1,4,3+i,2))
                        end 

                        for nind in noise_inds 
                            noise_hist[nind] ⊻= true 
                        end 
                        # noise_synds = Any[(4,7), (7,4)]
                        # noise_hist[1,:,:,:] .= synds_to_links(noise_synds,L,L)
                        println("∑noise = ",sum(noise_hist))
                    end 

                    if save_history 
                        # save the whole history (slower)
                        err_hist, synd_hist = tc_gate_applier!(copy(init_synds),copy(init_errs),gates,periods,noise_hist,gadgetnoise_dict,false,true)
                        # failed = detect_logical_failure(err_hist[end,:,:,:]) # look for logical failure 
                        failed = detect_both_logical_failures(err_hist[end,:,:,:]) # look for logical failure 
                        # failed = failure_test(synd_hist[end,:,:]) # look for anyons not being cleaned up in a particular way 

                    else 
                        # just keep track of the final configuration (faster)
                        errs, synds = tc_gate_applier!(copy(init_synds),copy(init_errs),gates,periods,noise_hist,gadgetnoise_dict,false,false)
                        # failed = detect_logical_failure(errs) # look for logical failure 
                        failed = detect_both_logical_failures(errs) # look for logical failure 
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
    params["mode"] = mode; params["model"] = model; params["gate"] = gate; params["L"] = L; params["l"] = l; params["squarenoise"] = squarenoise; params["bias"] = bias; params["gadgetnoise"] = gadgetnoise; params["nps"] = nps; params["logdist"] = logdist; params["calc_corrs"] = calc_corrs; params["measure_corrs_multiple"] = measure_corrs_multiple; params["calc_dw_stats"] = calc_dw_stats; params["rollstate"] = rollstate; params["ps"] = ps; params["depth"] = depth; params["samps"] = samps; params["pc"] = pc; params["thermal_periods"] = thermal_periods; params["periods"] = periods

    println("saving data to file: $fout")
    f = jldopen(fout,"w")
    for (key,object) in data 
        write(f,key,object)
    end 
    for (key,object) in params
        write(f,key,object)
    end 
    close(f)
    
end 

main()

