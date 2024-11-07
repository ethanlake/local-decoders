using ArgParse 
using Random 
using Statistics
using LinearAlgebra 

"""
this code contains functions that apply noisy versions of gates to input states, as well as several functions that aid in decoding and analyzing steady-state properties 
"""

# a simple [5,2,2] code: 
const H = [1 0 0 1 1; 0 1 0 0 1; 0 0 1 1 0] # parity check matrix 
const G = [1 1; 0 1; 1 0; 1 0; 0 1] # generator matrix 

### functions used for probing steady state properties ### 

function get_domainwall_stats(vec)
    """
    computes a histogram of domain wall sizes in a given binary string vec 
    """ 
    L = length(vec) 
    stats = zeros(Int,L)

    current_sign = vec[1]
    current_length = 1 
    first_domain = true 
    first_domain_length = 1 
    for i in 2:L 
        if vec[i] != current_sign  
            stats[current_length] += 1 
            current_sign = vec[i]
            if first_domain 
                first_domain_length = current_length 
                first_domain = false 
            end 
            current_length = 1 
        else 
            current_length += 1 
        end 
    end 

    stats[current_length] += 1 # last domain 
    if first_domain && vec[end] != vec[1] # yet another edge case 
        stats[1] += 1 
    end 

    if vec[end] == vec[1]
        if ~first_domain 
            stats[first_domain_length] -= 1 
            stats[current_length] -= 1 
            stats[current_length+first_domain_length] += 1 
        else 
            stats[end] = 1 
        end 
    end 

    return stats 
end 

function get_corrs(vec)
    """
    returns matrix of (non-connected) 2-point correlation functions between spins in a binary string vec 
    """ 
    lv = length(vec)
    cs = zeros(lv,lv)
    for i in 1:lv for j in 1:(i-1) 
        cij = 4 * (vec[i]-.5) * (vec[j]-.5) 
        cs[i,j] = cij; cs[j,i] = cij
    end end 
    return cs 
end 

### functions related to decoding ###

function maj(vec)
    """ 
    majority vote function --- works for any block size 
    """
    return Int(sum(vec) > length(vec)/2)
end 

function correct(vec)
    """
    performs error correction for the [5,2,3] code given above 
    returns: best guess for codeword given a noisy input vector 
    """
    corrected_vec = copy(vec)
    syndrome = transpose((H * vec) .% 2)
    if syndrome == [1 0 0]
        corrected_vec[1] = 1-corrected_vec[1]
    elseif syndrome == [0 1 0]
        corrected_vec[2] = 1-corrected_vec[2]
    elseif syndrome == [0 0 1]
        corrected_vec[3] = 1-corrected_vec[3]
    elseif syndrome == [1 0 1]
        corrected_vec[4] = 1-corrected_vec[4]
    elseif syndrome == [1 1 0]
        corrected_vec[5] = 1-corrected_vec[5]
    # either zero or > 1 errors have occured; do nothing 
    elseif sum(syndrome) > 0 
        corrected_vec = G[:,rand(1:2)] # replace with a random codeword if an error has occured? seems like a bad strat
    end 
    return corrected_vec
end 

function energy(vec,n) 
    """
    measures number of failed checks for the n=5 code 
    """
    eng = 0 
    ind = 0 
    lv = length(vec)
    while n*ind + n < lv 
        eng += sum(transpose(H * vec[n*ind+1:n*ind+n]) .%2) == 0 ? 0 : 1 
        ind += 1 
    end 
    println("eng ind = $ind")
    return eng 
end 

function xy_maj(spins) 
    """
    takes the continuous majority of 3 rotor spins 
    """
    s1 = spins[1]; s2 = spins[2]; s3 = spins[3]
    x = cos(s1) + cos(s2) + cos(s3)
    y = sin(s1) + sin(s2) + sin(s3)
    ϕ = atan(y,x)
    return ϕ
end 

function detect_logical_failure(errs) 
    """ 
    to be used on a state with no anyons on a torus. 
    returns true if there is nontrivial logical error (along either cycle), and false otherwise 
    """
    Lx,Ly,_ = size(errs)
    xparity = false; yparity = false  
    for i in 1:Lx  
        xparity ⊻= errs[i,1,2]
    end 
    for i in 1:Ly 
        yparity ⊻= errs[1,i,1]
    end 
    return xparity || yparity 
end 

function detect_both_logical_failures(errs) 
    """ 
    same as detect_logical_failure, but returns true only if errors occur along *both* cycles of the torus. used to select for corner errors. 
    """
    Lx,Ly,_ = size(errs)
    xparity = false; yparity = false  
    for i in 1:Lx  
        xparity ⊻= errs[i,1,2]
    end 
    for i in 1:Ly 
        yparity ⊻= errs[1,i,1]
    end 
    return xparity && yparity 
end 


function average_logicals(errs,L,dir) 
    """
    computes harrington's approximate decoder thing -- returns pm1 if even/odd holonomy along all loops, and 0 if random holonomy along all loops. dir ∈ {"x" "y" "xy"};
    if "x" ∈ dir, measures Z wilson lines stretching around the y cycle (supported on x-oriented links)
    if "y" ∈ dir, ... around the x cycle (supported on y-oriented links)
    """
    count = 0 
    for r in 1:L 
        if occursin("x",dir) 
            count += sum(errs[r,:,1])%2 
        end 
        if occursin("y",dir)
            count += sum(errs[:,r,2])%2 
        end 
    end 
    return 2*(0.5 - count / ((dir == "xy" ? 2 : 1)*L)) # for each cycle, +1 if no errors have occured and -1 if a logical error happened
end 

function recursive_majority(vec,n) 
    """
    computes recursive majority of a length n^l vector vec 
    """
    l = round(Int,log(n,length(vec)))
    previous_voted_vec = copy(vec)
    voted_vec = zeros(n^(l-1))
    for k in l:-1:1 
        voted_vec = zeros(n^(k-1))
        for m in 1:n^(k-1) 
            voted_vec[m] = maj(previous_voted_vec[n*(m-1)+1 : (n)*(m-1)+n])
        end 
        previous_voted_vec = copy(voted_vec)
    end 
    return voted_vec 
end      

################ 1D codes ##################

function tsirelson_gate_applier(model,init_state,gates,periods,noise_hist,bias,gadgetnoise,record_history) 
    """
    applies gates for 1d tsirelson codes. arguments: 

    model: a string determining the type of decoder to be applied, currently ∈ ["rep" "rep5" "k2"]
    init_state: input state 
    gates: list of gates to be applied 
    periods: number of repetitions of gate sequence 
    noise_hist: either a full spacetime history of noise locations (at which bit flips are performed), or a single Float. If the latter, applies iid noise with strength noise_hist. 
    bias: noise bias in case of iid noise. 
    gadgetnoise: if true and noisehist is a float, applies unbiased iid noisehist-bounded gadget noise (the outputs of *all* failed gadgets are uniformly random)
    record_history: if true, returns the full spin history over the given amount of time. else returns just the spins in the final state 
    """

    n = model ∈ ["rep5" "k2"] ? 5 : 3 
    L = size(init_state)[1]
    iid_noise = size(noise_hist)[1] == 1 # if true, does iid noise with prob p on each site. if false, uses the noise supplied as the noise_hist matrix
    p = noise_hist[1]

    depth = maximum(gates[:,3])
    ngates = size(gates)[1]
    history = record_history ? zeros(Int,periods*(depth)+1,L) : []
    if record_history history[1,:] = init_state end 
    spins = copy(init_state) # the array that will get updated in-place 

    time = 0; prev_time = 0; hist_time = 1 
    num_corr_samps = 0 

    for period in 1:periods 
        for i in 1:ngates 
            gtype = gates[i,1]; gloc = gates[i,2]
            prev_time = time 
            time = gates[i,3]  

            if time != prev_time # if we are just getting to a new timeslice, apply errors (this includes the initial state)
                if ~gadgetnoise # apply wire noise 
                    if iid_noise 
                        for j in 1:L 
                            if rand() < p # whether or not to apply noise 
                                spins[j] = rand() < (1+bias)/2 ? 1 : 0 
                            end 
                        end 
                    else 
                        for j in 1:L 
                            if noise_hist[hist_time,j] # spacetime location of a spin flip 
                                spins[j] = 1-spins[j]
                            end 
                        end 
                    end 
                end 

                # record history if desired 
                if record_history 
                    history[hist_time,:] .= spins
                end  

                hist_time += 1 
            end 

            if gtype == 1 # identity 
                nothing 
            elseif gtype == 2 # swap gate 
                s1 = spins[gloc]; s2 = spins[gloc+1]
                spins[gloc] = s2; spins[gloc+1] = s1 
            elseif gtype == 3 # error correction 
                if model ∈ ["rep" "rep5"]
                    spins[gloc:gloc+n-1] .= maj(spins[gloc:gloc+n-1])
                elseif model ∈ ["k2"] 
                    spins[gloc:gloc+n-1] .= correct(spins[gloc:gloc+n-1])
                end 
            end 
            if gadgetnoise
                if rand() < p # the gate fails -- all spins in its support get randomized 
                    gatesize = (gtype == 1 ? 1 : (gtype == 2 ? 2 : n)) 
                    spins[gloc:gloc+gatesize-1] .= rand(Bool,gatesize)
                end 
            end 
        end 
    end 

    # don't apply errors to the final state, since we are already adding noise at the first timestep
    if record_history 
        return history, [] # second argument are the syndromes; done just to have the same notation as the TC case 
    else 
        return spins, [] 
    end 
end 

function xy_majority_gate_applier(spins,gates,periods,ϵ,record_history) 
    """
    an attempt at producing a threshold (appropriately defined) for rotors. effective error rate however does *not* reduce with system size... kept here for posterity's sake 
    ϵ = noise amplitude: apply a random rotation by 2πϵ at each gate 
    """
    n = 3
    L = size(spins)[1]

    function ind(i)
        return (i+L-1)%L+1
    end 

    depth = maximum(gates[:,3]) 
    ngates = size(gates)[1]
    spin_history = record_history ? zeros(periods*depth+1,L) : zeros(1)

    time = 0; prev_time = 0; hist_time = 1 
    for period in 1:periods 
        for i in 1:ngates 
            gtype = gates[i,1]; x = gates[i,2]
            prev_time = time 
            time = gates[i,3] # starts at 2; 1 is the initial configuration 
            if time != prev_time # when we get to a new time step, or when we are at the last time step
                for xp in 1:L 
                    spins[xp] = (spins[xp] + 2π*ϵ*2*(rand()-.5) + 2π)%(2π)
                end 
                
                # take care of boundary conditions (if appropriate)
                if bcs == "fixed"
                    spins[1] = 0; spins[end] = 0 
                end 

                if record_history # take a snapshot of what has happened 
                    spin_history[hist_time,:] = spins
                end 

                hist_time += 1 
            end 

            # apply gates 
            if gtype == 1 # identity I 
                continue 
            elseif gtype == 2 # T0 gate -- swaps two spins  
                s1 = spins[x]; s2 = spins[x+1]
                spins[x+1] = s1; spins[x] = s2  
            elseif gtype == 3 # primitive majority vote 
                spins[x:x+2] .= xy_maj(spins[x:x+2])
            end 
        end 
    end 
    
    if record_history # take a snapshot of what has happened 
        spin_history[hist_time,:] = spins
    end     

    if record_history
        return spin_history 
    else 
        return spins 
    end 
end 


########### 2D codes ############

function twod_tsirelson_gate_applier(model,init_state,gates,periods,noise_hist,bias,gadgetnoise,record_history)
    """
    applies swaps for the 2d ising model -- really just for the purposes of visualization. gates can be fed in from the 1d tsirelson gate code. 
    arguments same as for the toric code gate applier.
    """

    p = noise_hist[1]
    spins = init_state
    L = size(init_state)[1] 
    n = 3
    function ind(i)
        return (i+L-1)%L+1
    end 
    depth = maximum(gates[:,3])
    spin_history = record_history ? zeros(Bool,2*periods*depth+1,L,L) : zeros(Bool,1)

    time = 0; prev_time = 0; hist_time = 1 

    function apply_noise()
        # beginning a new layer -- apply noise and measure syndromes etc. 
        for xp in 1:L 
            for yp in 1:L 
                if rand() < p
                    spins[xp,yp] = ran
                    if rand() < bias 
                        spins[xp,yp] = (rand() < bias ? true : false)
                    end 
                end 
            end 
        end 
    end 

    for period in 1:periods 
        for t in 1:depth 
            if record_history # take a snapshot of what has happened 
                spin_history[hist_time,:,:] = spins 
            end 
            hist_time += 1 
            apply_noise()

            # apply gates on the horizontal spins 
            zlayer = false 
            thesegates = gates[findall(time->time==t,gates[:,end]),:]
            for gateind in 1:size(thesegates)[1]
                thisgate = thesegates[gateind,:]
                gtype = thisgate[1]; x = thisgate[2]; xp1 = ind(x+1); xp2 = ind(x+2)
                if gtype != 1 # if we are not doing the identity... 
                    if gtype == 2 # swap 
                        for y in 1:L # loop over the columns 
                            s1 = spins[x,y]; s2 = spins[xp1,y]
                            spins[x,y] = s2; spins[xp1,y] = s1
                        end 
                    else # if we are doing a layer of 3x3 majority votes
                        for yblock in 1:round(Int,L/3)
                            y = 3*(yblock-1)+1  
                            spins[x:x+2,y:y+2] .= maj(spins[x:x+2,y:y+2])
                        end 
                        zlayer = true 
                    end 
                end 
            end 

            # apply things in the other direction if the current step consists of doing swaps + identities 
            if ~zlayer
                
                if record_history # take a snapshot of what has happened 
                    spin_history[hist_time,:,:] = spins 
                end 
                hist_time += 1 
                apply_noise()
                
                thesegates = gates[findall(time->time==t,gates[:,end]),:]
                for gateind in 1:size(thesegates)[1]
                    thisgate = thesegates[gateind,:]
                    gtype = thisgate[1]; y = thisgate[2]; yp1 = ind(y+1); yp2 = ind(y+2)
                    @assert gtype != 3 
                    if gtype != 1 
                        for x in 1:L 
                            s1 = spins[x,y]; s2 = spins[x,yp1]
                            spins[x,y] = s2; spins[x,yp1] = s1 
                        end 
                    end 
                end 
            end 
        end 
    end 
    # do the final update of the syndromes + store the final configuration at the last time step 

    if record_history # take a snapshot of what has happened 
        spin_history[hist_time,:,:] = spins 
    end     

    if record_history
        return spin_history, [] #[1:histind-1,:,:], err_history[1:histind-1,:,:,:] 
    else 
        return spins, [] 
    end 
end 

function tc_gate_applier!(synds,errs,gates,periods,noise_hist,gadgetnoise,squarenoise,measurementnoise,record_history)   
    """
    synds, errs: modified in-place 
    gates: list of level-0 gates to be applied 
    periods: number of times to repeat the gate set 
    noise_hist: if length 1, strength of iid noise; else, (time_steps x L x L x 4) array containing the error locations 
    iid noise has the following options: 
        gadgetnoise: if a dictionary (with labels [a,o] with a ∈ {1,2,3} ~ {I,T,M} the gate type and o ∈ {1,2} the orientation), does gadget-based noise and gives a list of the types of failures that can occur for each level-0 gadget. if true, does i) iid gadget noise with strength noise_hist[1] if length(noise_hist) = 1 and noise_hist[1] is a float, and ii) fails a particular sequence of gadgets (labeled by the integer order in which they appear in gates) if noise_hist is an integer-valued list (each gadget failure itself producing a random output)
        squarenoise: if true, noise acts on diagonal links (if false, only on vertical and horizontal links)
        measurementnoise: if true, syndrome measurements give wrong values with probability p (and wires are clean)
    record_history: if true, outputs a spacetime history of the errors being corrected 

    level-0 gate numbering conventions: 
    I: 0; T: 2; M: 3; R: 4
    each gate: [type x y orientation time]
    """

    Lx, Ly = size(synds)
    time = 0; prev_time = 0; hist_time = 1 

    iid_noise = size(noise_hist)[1] == 1 
    p = noise_hist[1]
    gadgetnoise_flag = isa(gadgetnoise,Dict) || gadgetnoise == true # true if gadget noise is being applied 
    iid_gadget_noise = ~isa(noise_hist,Vector{Int}) # if noise_hist a vector of ints, then the gadgets at positions noise_hist[i] fail. else, iid failures. 
    linknoise = ~(gadgetnoise_flag || measurementnoise)

    function indx(i) return (i+Lx-1)%Lx+1 end 
    function indy(i) return (i+Ly-1)%Ly+1 end 

    depth = maximum(gates[:,end]) 
    ngates = size(gates)[1]

    err_history = record_history ? zeros(Bool,periods*depth+1,Lx,Ly,2) : zeros(Bool,1)
    synd_history = record_history ? zeros(Bool,periods*depth+1,Lx,Ly) : zeros(Bool,1)

    function apply_noise()
        # beginning a new layer -- apply noise and measure syndromes etc. 
        if iid_noise 
            if p > 0 
                if ~squarenoise # straight-link-based noise 
                    for xp in 1:Lx for yp in 1:Ly for a in 1:2 
                        if rand() < p
                            errs[xp,yp,a] ⊻= true 
                        end 
                    end end end 
                else # noise that either acts on plaquettes (0box-bounded noise) or acts on *all* links (including diagonal ones)
                    for xp in 1:Lx for yp in 1:Ly 
                        # treat straight and diagonal links equally 
                        # straight links: 
                        if rand() < p errs[xp,yp,1] ⊻= true end 
                        if rand() < p errs[xp,yp,2] ⊻= true end 
                        # diagonal links: 
                        if rand() < p errs[xp,yp,1] ⊻= true; errs[indx(xp+1),yp,2] ⊻= true end # up and to the right diagonal 
                        if rand() < p errs[xp,yp,1] ⊻= true; errs[xp,yp,2] ⊻= true end # up and to the left diagonal (in the same plaquette as the up and to the right one)

                    end end 
                end 
            end 
        else # given an explicit noise history 
            for xp in 1:Lx for yp in 1:Ly 
                # straight links: 
                if noise_hist[hist_time,xp,yp,1] errs[xp,yp,1] ⊻= true end 
                if noise_hist[hist_time,xp,yp,2] errs[xp,yp,2] ⊻= true end 
                # diagonal links: 
                if noise_hist[hist_time,xp,yp,3] errs[xp,yp,1] ⊻= true; errs[indx(xp+1),yp,2] ⊻= true end # up and to the right diagonal 
                if noise_hist[hist_time,xp,yp,4] errs[xp,yp,1] ⊻= true; errs[xp,yp,2] ⊻= true end # up and to the left diagonal 

            end end 
        end 
    end 

    function update_synds()
        for xp in 1:Lx  
            xpm1 = indx(xp-1)
            for yp in 1:Ly 
                ypm1 = indy(yp-1)
                synds[xp,yp] = errs[xp,yp,1] ⊻ errs[xp,yp,2] ⊻ errs[xpm1,yp,1] ⊻ errs[xp,ypm1,2]
            end 
        end 
    end 

    # the following functions applyT and applyM update both link spins and syndromes 
    function applyT(xloc,yloc,dir)
        control_synd = dir == 1 ? synds[indx(xloc+1),yloc] : synds[xloc,indy(yloc+1)]
        if measurementnoise control_synd ⊻= (rand() < p) end # flip syndrome measurement if measurementnoise occurs 
        if control_synd 
            if dir == 1 
                errs[xloc,yloc,1] ⊻= 1; errs[indx(xloc+1),yloc,1] ⊻= 1 
                synds[xloc,yloc] ⊻= 1; synds[indx(xloc+2),yloc] ⊻= 1 
            else  # y-oriented T gate 
                errs[xloc,yloc,2] ⊻= 1; errs[xloc,indy(yloc+1),2] ⊻= 1 
                synds[xloc,yloc] ⊻= 1; synds[xloc,indy(yloc+2)] ⊻= 1 
            end 
        end 
    end 
    
    function applyM(xloc,yloc,dir) 
        control_synd_1 = dir == 1 ? synds[indx(xloc+1),yloc] : synds[xloc,indy(yloc+1)] 
        control_synd_2 = dir == 1 ? synds[indx(xloc+2),yloc] : synds[xloc,indy(yloc+2)]
        if measurementnoise
            control_synd_1 ⊻= (rand() < p); control_synd_2 ⊻= (rand() < p)
        end 
        if control_synd_1 && control_synd_2 
            if dir == 1 # x facing gates 
                errs[indx(xloc+1),yloc,1] ⊻= 1
                synds[indx(xloc+1),yloc] ⊻= 1; synds[indx(xloc+2),yloc] ⊻= 1 # can't just set to false since an anyon pair can be created as the result of two faulty measurements 
            else 
                errs[xloc,indy(yloc+1),2] ⊻= 1
                synds[xloc,indy(yloc+1)] ⊻= 1; synds[xloc,indy(yloc+2)] ⊻= 1 
            end 
        end 
    end 
    
    # define functions that apply blocks of gates that appear in the action of R -- applied over the WHOLE system (only used when primitive R0 = true and Lx = Ly)
    Lo3 = Lx ÷ 3 
    function horz_Ms()
        for i in 1:Lo3 for k in 1:Lo3 #
            x = (i-1)*3 + 1; y = (k-1)*3 + 1
            for j in 0:2 applyM(x,y+j,1) end
        end end  
    end 
    function vert_Ms()
        for i in 1:Lo3 for k in 1:Lo3 #
            x = (i-1)*3 + 1; y = (k-1)*3 + 1
            for j in 0:2 applyM(x+j,y,2) end 
        end end 
    end 
    function top_Ts()
        for i in 1:Lo3 for k in 1:Lo3 #
            x = (i-1)*3 + 1; y = (k-1)*3 + 1
            for j in 0:2 applyT(x+j,y+1,2) end 
        end end 
    end 
    function bottom_Ts()
        for i in 1:Lo3 for k in 1:Lo3 #
            x = (i-1)*3 + 1; y = (k-1)*3 + 1
            for j in 0:2 applyT(x+j,y,2) end
        end end  
    end 
    function right_Ts()
        for i in 1:Lo3 for k in 1:Lo3 #
            x = (i-1)*3 + 1; y = (k-1)*3 + 1
            for j in 0:2 applyT(x+1,y+j,1) end 
        end end 
    end 
    function left_Ts()
        for i in 1:Lo3 for k in 1:Lo3 #
            x = (i-1)*3 + 1; y = (k-1)*3 + 1
            for j in 0:2 applyT(x,y+j,1) end
        end end  
    end 
    
    function top_Tsp() top_Ts(); horz_Ms() end 
    function bottom_Tsp() bottom_Ts(); horz_Ms() end 
    function left_Tsp() left_Ts(); vert_Ms() end 
    function right_Tsp() right_Ts(); vert_Ms() end 

    function applyR() # applies a *full layer of R0s* across the system (updating both the spins and the syndromes)

        function vert_part()
            vert_Ms()
            bottom_Tsp()
            vert_Ms()
            top_Tsp()
        end 

        function horz_part()                    
            horz_Ms()
            left_Tsp()
            horz_Ms()
            right_Tsp()
        end 

        vert_part(); vert_part() # does full EC 
        horz_part(); horz_part() 

    end 

    gadget_failure_links = falses(Lx,Ly,2); gadget_failure_synds = falses(Lx,Ly) # will only be used if applying gadget-based noise 
    for period in 1:periods 
        for i in 1:ngates 
            gtype = gates[i,1]; o = gates[i,2]; x = gates[i,3]; y = gates[i,4]

            prev_time = time 
            time = gates[i,end] 
            if time != prev_time # when we get to a new time step, or when we are at the last time step
                if linknoise # if not doing gadget-based or measurement noise 
                    apply_noise(); update_synds()
                elseif gadgetnoise_flag # if doing gadget-based noise, will keep track of a pattern of failures for a given layer, and then add on that failure at the end of the time slice (the noise occurs *after* the gadgets operate cleanly)
                    errs .⊻= gadget_failure_links; synds .⊻= gadget_failure_synds 
                    gadget_failure_links = falses(Lx,Ly,2); gadget_failure_synds = falses(Lx,Ly) # reset things for the next layer  
                end 

                if record_history # take a snapshot of what has happened 
                    err_history[hist_time,:,:,:] = errs
                    synd_history[hist_time,:,:] = synds 
                end 
                hist_time += 1 
            end 

            # apply gates 
            if gtype == 1 # identity I 
                continue 
            elseif gtype == 2 # T0 gate 
                applyT(x,y,o)
            elseif gtype == 3 # M0 gate -- acts on 3 sites and matches two anyons in the middle (first part of maj vote)
                applyM(x,y,o)
            elseif gtype == 4 # primitive majority vote on the 3x3 square whose bottom left corner is at (x,y). to apply this in one shot, we need to do the WHOLE system in the same time, in parallel. since Rs only happen in a giant layer (they are never isolated), we can call applyR()---which does R over the whole system---only when we are in the bottom-left corner.
                if x == 1 && y == 1 
                    applyR()
                end 
            end 

            if gadgetnoise_flag 
                if (iid_gadget_noise && rand() < p) || i + (period-1)*ngates ∈ noise_hist 
                    errors = rand(gadgetnoise[gtype,o]) # get a random gadget failure 
                    for syndloc in errors 
                        # for each syndrome, create an error path going from the bottom left corner of the gadget to the syndrome location
                        xsynd = indx(x+syndloc[1]); ysynd = indy(y+syndloc[2])
                        gadget_failure_synds[xsynd,ysynd] ⊻= true # update the syndromes 
                        # now update the links 
                        for dx in 0:(syndloc[1]-1) # first (horizontal) segment of error path 
                            gadget_failure_links[indx(x+dx),y,1] ⊻= true 
                        end 
                        for dy in 0:(syndloc[2]-1) # second (vertical) segment of error path 
                            gadget_failure_links[xsynd,indy(y+dy),2] ⊻= true 
                        end 
                    end 
                end 
            end 
        end 
    end 
    # do the final update of the syndromes + store the final configuration at the last time step 
    update_synds()

    if record_history # take a snapshot of what has happened 
        err_history[hist_time,:,:,:] = errs
        synd_history[hist_time,:,:] = synds 
    end     

    if record_history
        return err_history, synd_history 
    else 
        return errs, synds 
    end 
end 

function syndrome_tc_gate_applier!(synds,gates,eloc,error)
    """
    applies gates to a set of syndromes, with a single level-0 gadget error occuring at spacetime location eloc = [t x y] and with syndromes given by error (a list of tuples providing the locations of the errors relative to (x,y))
    is slightly faster since doesn't keep track of the link spins; currently used only when testing nilpotence 
    
    gate numbering conventions as in tc_gate_applier 
    """
    Lx, Ly = size(synds)
    n = 3

    function indx(i) return (i+Lx-1)%Lx+1 end 
    function indy(i) return (i+Ly-1)%Ly+1 end 

    depth = maximum(gates[:,end]) 
    ngates = size(gates)[1]

    function applyT(xloc,yloc,dir)
        if dir == 1 # x-oriented T gate  
            if synds[indx(xloc+1),yloc] # 
                synds[xloc,yloc] ⊻= 1; synds[indx(xloc+2),yloc] ⊻= 1 
            end 
        else  # y-oriented T gate 
            if synds[xloc,indy(yloc+1)] # 
                synds[xloc,yloc] ⊻= 1; synds[xloc,indy(yloc+2)] ⊻= 1 
            end 
        end 
    end 
    
    function applyM(xloc,yloc,dir) 
        if dir == 1 # x facing gates 
            if (synds[indx(xloc+1),yloc] && synds[indx(xloc+2),yloc])
                synds[indx(xloc+1),yloc] = 0; synds[indx(xloc+2),yloc] = 0 # don't forget to update the syndromes! 
            end 
        else # y-oriented gate 
            if (synds[xloc,indy(yloc+1)] && synds[xloc,indy(yloc+2)]) # don't have to worry about PBC since () will always be at least 3 sites from the boundary
                synds[xloc,indy(yloc+1)] = 0; synds[xloc,indy(yloc+2)] = 0 
            end 
            # end 
        end 
    end 
    
    # define functions that apply blocks of gates that appear in the action of R -- applied over the WHOLE system (only used when primitive R0s are used and Lx = Ly)
    Lo3 = Lx ÷ 3 
    function horz_Ms()
        for i in 1:Lo3 for k in 1:Lo3 #
            x = (i-1)*3 + 1; y = (k-1)*3 + 1
            for j in 0:2 applyM(x,y+j,1) end
        end end  
    end 
    function vert_Ms()
        for i in 1:Lo3 for k in 1:Lo3 #
            x = (i-1)*3 + 1; y = (k-1)*3 + 1
            for j in 0:2 applyM(x+j,y,2) end 
        end end 
    end 
    function top_Ts()
        for i in 1:Lo3 for k in 1:Lo3 #
            x = (i-1)*3 + 1; y = (k-1)*3 + 1
            for j in 0:2 applyT(x+j,y+1,2) end 
        end end 
    end 
    function bottom_Ts()
        for i in 1:Lo3 for k in 1:Lo3 #
            x = (i-1)*3 + 1; y = (k-1)*3 + 1
            for j in 0:2 applyT(x+j,y,2) end
        end end  
    end 
    function right_Ts()
        for i in 1:Lo3 for k in 1:Lo3 #
            x = (i-1)*3 + 1; y = (k-1)*3 + 1
            for j in 0:2 applyT(x+1,y+j,1) end 
        end end 
    end 
    function left_Ts()
        for i in 1:Lo3 for k in 1:Lo3 #
            x = (i-1)*3 + 1; y = (k-1)*3 + 1
            for j in 0:2 applyT(x,y+j,1) end
        end end  
    end 
    
    function top_Tsp() top_Ts(); horz_Ms() end 
    function bottom_Tsp() bottom_Ts(); horz_Ms() end 
    function left_Tsp() left_Ts(); vert_Ms() end 
    function right_Tsp() right_Ts(); vert_Ms() end 

    function applyR() # applies a *full layer of R0s* across the system  

        function vert_part()
            vert_Ms()
            bottom_Tsp()
            vert_Ms()
            top_Tsp()
        end 

        function horz_part()                    
            horz_Ms()
            left_Tsp()
            horz_Ms()
            right_Tsp()
        end 

        vert_part(); vert_part()   
        horz_part(); horz_part() 

    end 

    time = 0; prev_time = 0; hist_time = 1 
    for i in 1:ngates 
        gtype = gates[i,1]; o = gates[i,2]; x = gates[i,3]; y = gates[i,4]

        prev_time = time 
        time = gates[i,end] # starts at 2; 1 is the initial configuration 
        if time != prev_time # when we get to a new time step, or when we are at the last time step
            if hist_time == eloc[1]
                ex = eloc[2]; ey = eloc[3]
                for syndloc in error 
                    synds[ex+syndloc[1],ey+syndloc[2]] ⊻= true                 
                end 
            end         
            hist_time += 1 
        end 

        # apply gates 
        if gtype == 1 # identity I 
            continue 
        elseif gtype == 2 # T0 gate 
            applyT(x,y,o)
        elseif gtype == 3 # M0 gate -- acts on 3 sites and matches two anyons in the middle (first part of maj vote)
            applyM(x,y,o)
        elseif gtype == 4 # primitive majority vote on the 3x3 square whose bottom left corner is at (x,y). to apply this in one shot, we need to do the WHOLE system in the same time, in parallel. since Rs only happen in a giant layer (they are never isolated), we can call applyR()---which does R over the whole system---only when we are in the bottom-left corner.
            if x == 1 && y == 1 
                applyR()
            end 
        end 
    end 

    return synds 
end 

function master_gate_applier(model,init_state,init_synds,gates,periods,noise_hist,bias,gadgetnoise,squarenoise,measurementnoise,record_history)
    """
    master function for applying gates -- just helps streamline code in other places 
    
    model: type of code being run, ∈ {"rep" "rep5" "k2" "twod_rep" "tc"}
    init_state: initial spin configuration 
    init_synds: initial syndrome configuration; used only for tc 
    gates: list of gates to be applied 
    periods: number of times to apply the full gate set 
    noise_hist: if an array, a spacetime noise array detailing the points at which spin flips happen. if a number, the strength of iid noise 
    bias: bias of noise in iid case 
    gadgetnoise: if true, does gadget noise model 
    squarenoise: if true, diagonal 0-links fail at the same probability as straight ones 
    measurementnoise: if true, errors (only) occur by way of faulty syndrome measurements 
    record_history: if true, returns full spin history over course of evolution (else just returns final state)

    returns: final_state if not record_history, else (final_state, spin_history)
    """

    # 1d tsirelson codes 
    if model ∈ ["rep" "rep_5bit" "k2"] 
        return tsirelson_gate_applier(model,init_state,gates,periods,noise_hist,bias,gadgetnoise,record_history)  
            
    # 2d tsirelson for majority vote 
    elseif model == "twod_rep"
        return twod_tsirelson_gate_applier(model,init_state,gates,periods,noise_hist,bias,gadgetnoise,record_history)

    # toric code 
    elseif model == "tc"
        return tc_gate_applier!(init_synds,init_state,gates,periods,noise_hist,gadgetnoise,squarenoise,measurementnoise,record_history)    
    end 
end 