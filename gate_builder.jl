using Random 
using Statistics
using LinearAlgebra 
using JLD2, FileIO

"""
this code is used to create circuits for 
1) Tsirelson's 1d self-correcting memory and generalizations thereof (defined for the concatenation of any classical code)
2) a Tsirelson-inspired local decoder for the 2d toric code
when run as main, each circuit is generated as a list of gates and saved as a .jld2 file; see below for conventions on how gates are labeled. 
"""

function get_gate_size(l,model,gate)
    """
    helper function that stores the spatial sizes of various gates
    """
    Lx = 1; Ly = 1 
    n = model ∈ ["k2" "rep_5bit"] ? 5 : 3 
    if model == "tc" 
        if gate ∈ ["R"] Lx = n^(l+1); Ly = n^(l+1) 
        elseif gate ∈ ["I" "IR" "i"] Lx = n^l; Ly = n^l  
        elseif gate ∈ ["MyR" "my"] Lx = n^l; Ly = 3*n^l   
        elseif gate ∈ ["M" "MR" "MxR" "mx"] Lx = 3*n^l; Ly = n^l 
        elseif gate ∈ ["TyR" "ty"] Lx = n^l; Ly = 2*n^l  
        elseif gate ∈ ["T" "TR" "TxR" "tx"] Lx = 2*n^l; Ly = n^l 
        end 
    elseif model == "twod_rep" 
        if gate == "R" Lx = n^(l+1); Ly = n^(l+1) end 
    elseif model ∈ ["rep" "rep_5bit" "k2"] # 1d codes 
        if gate ∈ ["X" "I"] Lx = n^l  
        elseif gate ∈ ["Y" "T"] Lx = 2n^l 
        elseif gate ∈ ["Z" "R" "M"] Lx = n^(l+1) 
        end 
    end   
    if Lx == 1 && Ly == 1 
        println("gate $gate not found...")
    end 
    return Lx, Ly 
end 

function permute_circuit(n,swap) 
    """
    subroutine called for performing swaps of blocks of bits / cyclings of tensor factors 

    n = size of bit block 
    if swap = true
        consider 2n points on a line. builds a circuit which swaps the first n points with the last n
    if swap = false
        consider n^2 points on a line with coordinates (i,j) = i-1 + (j-1)*n. builds a circuit that swaps (i,j) with (j,i).  

    returns: (list of gates in circuit, depth of gates)
    """

    vec = [i for i in 1:2n]
    if ~swap 
        vec = [i for i in 1:n^2]
        for i in 1:n 
            for j in 1:n 
                ind = (i-1)*n + j 
                swapind = (j-1)*n + i 
                vec[ind] = swapind 
            end 
        end 
    else 
        vec[1:n] = [i for i in n+1:2n]
        vec[n+1:end] = [i for i in 1:n]
    end 
    maxind = ~swap ? n^2-1 : 2n-1 # the largest index on which a swap gate can be placed 
    
    parity = 1 # start from 1 wolog if we are doing the tensor cycling 
    if swap # first swap gate should act on the last part of the first block 
        parity = n%2 == 1 ? 0 : 1 
    end 
    depth = 0
    unsorted = true 
    gates = []
    while unsorted 
        added_gate = false 
        if parity == 1 # take care of the identity gate acting on the first site if we start from the second 
            push!(gates,[1 1 depth])
        end 
        for k in (1+parity):2:maxind
            v1 = vec[k]; v2 = vec[k+1]
            if v1 > v2 
                vec[k] = v2; vec[k+1] = v1 
                push!(gates,[k 2 depth]) # swap gate at location k 
                added_gate = true # added a swap gate 
            else # apply two identity gates 
                push!(gates,[k 1 depth])
                push!(gates,[k+1 1 depth])
            end
            if k == maxind-1 # take care of the identity gate acting on the last site 
                push!(gates,[maxind+1 1 depth])
            end 
        end 
    
        parity = 1-parity 
        if ~added_gate
            unsorted = false 
        else 
            depth += 1 
        end 
    end 
    return gates, depth
end 

function tsirelson_gate_builder(model,l,gate) 
    """
    builds a time-ordered list of gates making up renormalized l-scale gates in a Tsirelson automaton. 

    model: choice of code, ∈ {"rep" "rep_5bit" "k2"}
    l: level of gadget 
    gate: label of gate whose circuit is to be created (see below for options)

    returns: time-ordered list of microscopic gates that make up the renormalized gate circuit.

    notation for microscopic gates: 
    X0 (identity) gate at i: [1 i]
    Y0 (swap)     gate at i,i+1: [2 i]
    Z0 (EC)       gate at i,...,i+n: [3 i]
    """

    n = model ∈ ["k2" "rep_5bit"] ? 5 : 3 # block size 
    L, _ = get_gate_size(l,model,gate)

    gates = []

    local_times = ones(Int,L) 

    # the circuits needed for doing the swaps -- neither one of them includes identity gates 
    cycle_circuit, cycle_depth = permute_circuit(n,false) 
    swap_circuit, swap_depth = permute_circuit(n,true)

    nx = 2n-1 # number of rows of Xn-1s that we stack in the definition of Xn 
    
    function decoder(vec,k) # here input vector should be n^k dimensional 
        lv = n^k 
        if k > 1 
            zn(vec,k-1)
            for j in 1:n 
                decoder(vec[1+(j-1)*n^(k-1):j*n^(k-1)],k-1)
            end 
        else 
            zn(vec,0)
        end 
    end 

    function xn(vec,k) 
        """ 
        renormalized identity: acts on vectors of length n^k 
        """
        if k > 1 
            zn(vec,k-1)
            for i in 1:nx # apply the columns of xs 
                for j in 1:n 
                    xn(vec[1+(j-1)*n^(k-1):j*n^(k-1)],k-1)
                end 
            end 
        elseif k == 1 
            zn(vec,0)
            for i in 1:nx # apply the column of xns 
                for j in 1:n 
                    xn(vec[j],0)
                end 
            end 
        elseif k == 0
            push!(gates,[1 vec[1] local_times[vec[1]]])
            local_times[vec[1]] += 1 
        end 
    end 

    function yn(vec,k)
        """
        renormalized swap: acts on vectors of length 2n^k 
        """
        if k > 0 
            # start with EC 
            zn(vec[1:n^k],k-1)
            zn(vec[n^k+1:end],k-1)

            # swaps 
            for gate in swap_circuit
                if gate[3] < swap_depth # the last row of the swap circuit is just trivial identities 
                    if gate[2] == 1 # identity 
                        xn(vec[1+(gate[1]-1) * n^(k-1):gate[1] * n^(k-1)],k-1) 
                    elseif gate[2] == 2 # swap 
                        yn(vec[1+(gate[1]-1)*n^(k-1):(gate[1]+1) * n^(k-1)],k-1)
                    end 
                end 
            end 
        else 
            push!(gates,[2 vec[1] local_times[vec[1]]])
            local_times[vec[1]] += 1; local_times[vec[2]] += 1 
        end 
    end 

    function zn(vec,k) 
        """
        renormalized EC: acts on vectors of length n^(k+1) 
        """
        if k > 0 
            for j in 1:n 
                zn(vec[1+(j-1)*n^k:j*n^k],k-1)
            end 
            
            for gate in cycle_circuit
                if gate[3] < cycle_depth # the last row of the swap circuit is just trivial identities 
                    if gate[2] == 1 # gate[1] runs from 1 to n^2 
                        xn(vec[1+(gate[1]-1) * n^(k-1):gate[1] * n^(k-1)],k-1) 
                    elseif gate[2] == 2 # swap 
                        yn(vec[1+(gate[1]-1) * n^(k-1):(gate[1]+1) * n^(k-1)],k-1)
                    end 
                end 
            end 

        elseif k == 0 
            push!(gates, [3 vec[1] local_times[vec[1]]])
            for j in 1:n 
                local_times[vec[j]] += 1
            end 
        end 
    end     

    vec = [i for i in 1:L] 
    
    if gate == "X"
        xn(vec,l)
    elseif gate == "Y"
        yn(vec,l)
    elseif gate == "Z" 
        zn(vec,l)
    elseif gate == "decoder"
        decoder(vec,l)
    elseif gate == "long_decoder"
        decoder(vec,l)
        xn(vec,l)
    elseif gate == "two_decoders"
        decoder(vec,l); decoder(vec,l)
    end 

    ngates = length(gates)
    depth = maximum(local_times) - 1 

    # turn the gates into an array (rather than a list) sorted in time order --- and check that local_times gives a consistent flat top
    gates_array = zeros(Int,ngates,3) 
    for i in 1:ngates gates_array[i,:] = gates[i][:] end 
    sort_inds = sortperm(gates_array[:,end])
    gates_array = gates_array[sort_inds,:]
    depth = maximum(local_times) - 1
    @assert depth == minimum(local_times)-1 # need all local times to be equal at the end of the gate construction 

    return gates_array, depth, L 
end 

function twod_tsirelson_gate_builder(l)
    """
    tsirelson-inspired 2d majority voting --- mostly just used for visualization purposes. currently only builds the EC gadget. 
    format for gates: 
    [gtype gate_orientation xcoord ycoord local_time]
    """

    n = 3 
    gates = []
    numrs = 1 

    tsirelson_swap_gates, swap_depth = permute_circuit(n,false)
    Lx, Ly = get_gate_size(l,"twod_rep","R")
       
    local_times = ones(Int,Lx,Ly)

    function indx(x) # for dealing with PBC 
        return (x + Lx - 1)%Lx + 1
    end 
    function indy(x)  
        return (x + Ly - 1)%Ly + 1
    end 

    nxx = 5 # if not putting Ms between Ts 

    function idn(xraw,yraw,o,k) # does a stack of 3 identity gates along the chosen orientation 
        x = indx(xraw); y = indy(yraw)
        xf = o == 1 ? 1 : 0; yf = 1 - xf 
        km = k-1
        if k > 0 
            sb = n^(k-1)
            bb = n^k 
            for j in 0:2 
                for repind in 1:numrs 
                    rn(x + xf*j*bb,y + yf*j*bb,km)
                end 
                for i in 0:2 
                    for repind in 1:nxx # depth of the T gate -- 5*2, where 2 comes from the perp. Ms 
                        idn(x + xf*j*bb + yf*i*sb,y + yf*j*bb + xf*i*sb,o,km)
                    end 
                end 
            end 
        else
            push!(gates,[1 o x y local_times[x,y]]) 
            for i in 0:2 # in a stack 
                local_times[indx(x+i*xf),indy(y+i*yf)] += 1 
            end 
        end 
    end 

    function tn(xraw,yraw,o,k) # either a 2x3 or 3x2 block, depending on orientation. if o = 1 then a 3x2 block and the Ts run along the x direction; reversed if o = 2 
        x = indx(xraw); y = indy(yraw)
        xf = o == 1 ? 1 : 0; yf = 1-xf 

        if k > 0 
            km = k-1
            bb = n^k; sb = n^(k-1) # big and small block sizes 

            # first stack of Rs 
            for repind in 1:numrs
                for d in 0:2 
                    rn(x + xf*d*bb,y + yf*d*bb,km) 
                    rn(x + yf*bb + xf*d*bb,y + xf*bb + yf*d*bb,km) 
                end 
            end 

            # second layer: (T in center)
            for j in 0:2 # act on the three big columns 
                xp = x + bb*j*xf # horiz position -- = x if o = 2 
                yp = y + bb*j*yf # vertical position -- = y if o = 1 
                idn(xp,yp,o,km); idn(xp + yf*sb,yp + xf*sb,o,km) 
                tn(xp + 2sb*yf,yp + 2sb*xf,o,km) # T in center 
                idn(xp + yf*4sb,yp + xf*4sb,o,km); idn(xp + yf*5sb,yp + xf*5sb,o,km)  
            end 

            # third layer: (two Ts)
            for j in 0:2 # 
                xp = x+bb*j*xf
                yp = y+bb*j*yf
                idn(xp,yp,o,km); tn(xp + yf*sb,yp + xf*sb,o,km) 
                tn(xp + yf*3sb,yp + xf*3sb,o,km); idn(xp + yf*5sb,yp + xf*5sb,o,km) 
            end 

            # fourth layer: (all Ts)
            for j in 0:2 # 
                xp = x+bb*j*xf
                yp = y+bb*j*yf
                tn(xp,yp,o,km)
                tn(xp + yf*2sb,yp + xf*2sb,o,km)
                tn(xp + yf*4sb,yp + xf*4sb,o,km) 
            end 

            # fifth layer: same as third 
            for j in 0:2 # 
                xp = x+bb*j*xf
                yp = y+bb*j*yf
                idn(xp,yp,o,km); tn(xp + yf*sb,yp + xf*sb,o,km) 
                tn(xp + yf*3sb,yp + xf*3sb,o,km); idn(xp + yf*5sb,yp + xf*5sb,o,km) 
            end 

            # # sixth layer: same as second 
            for j in 0:2 # act on the three big columns 
                xp = x + bb*j*xf # horiz position -- = x if o = 2 
                yp = y + bb*j*yf # vertical position -- = y if o = 1 
                idn(xp,yp,o,km); idn(xp + yf*sb,yp + xf*sb,o,km) 
                tn(xp + 2sb*yf,yp + 2sb*xf,o,km) # T in center 
                idn(xp + yf*4sb,yp + xf*4sb,o,km); idn(xp + yf*5sb,yp + xf*5sb,o,km)  
            end 

        elseif k == 0 # 
            push!(gates,[2 o x y local_times[x,y]])
            for i in 0:2
                local_times[indx(x + xf*i),indy(y + yf*i)] += 1
            end 
            for i in 0:2
                local_times[indx(x + xf*i + yf),indy(y + yf*i + xf)] += 1
            end 
        end 
    end

    # for the current gate this will only ever be called at zeroth level 
    function mn(xraw,yraw,o;flag="123") # "match" gate, acts on 3x3 block; x-oriented gates if o = 1, y-oriented if o = 2 
        x = indx(xraw); y = indy(yraw)
        xf = o == 1 ? 1 : 0; yf = 1-xf 
        op = o == 1 ? 2 : 1 

        f1 = occursin("1",flag); f2 = occursin("2",flag); f3 = occursin("3",flag)
        # flag determines the parts of the gate that we include, determined by looking in the direction of o (we only ever displace M_a gates along the a direction)

        if f1
            push!(gates,[3 o x y local_times[x,y]]) # primitive M0 acts on 3 sites, just to make things notationally consistent with higher layers -- x y labels the leftmost / bottommost point 
        end 
        # update local times 
        for i in 0:2
            if f1
                local_times[ind(x+yf*i),ind(y+xf*i)] += 1
            end 
            if f2
                local_times[ind(x+xf+yf*i),ind(y+yf+xf*i)] += 1 
            end 
            if f3
                local_times[ind(x+2*xf+yf*i),ind(y+2yf+xf*i)] += 1 
            end 
        end 
    end 

    function perp_mns(x,y,o) # applys a staggered depth-3 layer of Ms with orientation o
        sb = 1 # so sb = 1 is lattice spacing at smallest scale 
        if o == 1 
            mn(x,y,1)
            mn(x+sb,y,1,flag="12")
            mn(x-2sb,y,1,flag="3")
            mn(x+2sb,y,1,flag="1")
            mn(x-sb,y,1,flag="23")
        else 
            mn(x,y,2)
            mn(x,y+sb,2,flag="12")
            mn(x,y-2sb,2,flag="3")
            mn(x,y+2sb,2,flag="1")
            mn(x,y-sb,2,flag="23")
        end 
    end 

    function rn(xraw,yraw,k) # "renormalize" gate (majority vote): size of input vector is 3n^k
        x = indx(xraw); y = indy(yraw)
        km = k-1 

        bb = n^k; sb = k > 0 ? n^(k-1) : 0
        
        if k > 0 
            # tile things with R0s 
            for d in 1:numrs
                for i in 0:2 for j in 0:2 
                    rn(x+i*bb,y+j*bb,km)
                end end 
            end 

            function apply_tsirelson_swaps(o)
                xf = (o == 1 ? 1 : 2); yf = 1 - xf 
                lasttime = tsirelson_swap_gates[1][end]
                for thisgate in tsirelson_swap_gates
                    gtype = thisgate[2]; gloc = thisgate[1]; gtime = thisgate[3]
                    if thisgate[3] < swap_depth # the swap gates contain a redundant layer of identity gates at the end which we want to ignore 
                        for d in 0:2 
                            if o == 1 
                                gtype == 1 ? idn(x+(gloc-1)*sb,y+d*bb,2,km) : tn(x+(gloc-1)*sb,y+d*bb,2,km)
                            else 
                                gtype == 1 ? idn(x+d*bb,y+(gloc-1)*sb,1,km) : tn(x+d*bb,y+(gloc-1)*sb,1,km)
                            end 
                        end 
                    end 
                    lasttime = gtime
                end 
            end 
        
            # apply x-facing tsirelson circuit 
            apply_tsirelson_swaps(1)
            # apply y-facing tsirelson circuit 
            apply_tsirelson_swaps(2)
            
        else # apply regular R0 made from constituent gates 

            # horizontal Ms: 
            mn(x,y,1)

            # left Ts 
            tn(x,y,2,k); idn(x+2,y,2,k)
            
            # vertical_Ms()
            perp_mns(x,y,2)

            # horiz Ms 
            mn(x,y,1)
            
            # right Ts 
            idn(x,y,2,k); tn(x+1,y,2,k) 
            
            perp_mns(x,y,2)

            # horizontal Ms: 
            mn(x,y,1)

            # vertical Ms: 
            mn(x,y,2)

            # bottom Ts 
            tn(x,y,1,k); idn(x,y+2,1,k)  
            
            # horizontal_Ms()
            perp_mns(x,y,1)

            # vert Ms 
            mn(x,y,2)

            # top Ts 
            idn(x,y,1,k);  tn(x,y+1,1,k)
            
            # horiz Ms 
            # horizontal_Ms()
            perp_mns(x,y,1)

            # vertical Ms: 
            mn(x,y,2)

        end 

    end     
    
    rn(1,1,l)

    ngates = length(gates)
    depth = maximum(local_times) - 1 

    # turn the gates into an array (rather than a list) sorted in time order --- and check that local_times gives a consistent flat top
    gates_array = zeros(Int,ngates,5) # gtype x y orientation time 
    for i in 1:ngates gates_array[i,:] = gates[i][:] end 
    sort_inds = sortperm(gates_array[:,end])
    gates_array = gates_array[sort_inds,:]
    depth = maximum(local_times) - 1
    @assert depth == minimum(local_times)-1 # need all local times to be equal at the end of the gate construction 
    return gates_array, depth, Lx, Ly 
end 

function tc_gate_builder(l,gate,all_boundaries)
    """
    all_boundaries: if true, includes level-0 gadgets on the N and E boundaries of the gate to be constructed 
    
    format for gates: 
    [gtype gate_orientation xcoord ycoord local_time]
    I: 1
    T: 2 
    M: 3 
    R: 4 (primitive R0)
    
    depths of gates: (simplest to set d(T) = d(X) = d(M) at each layer)
    d(T1) = 5 + nxt 
    d(M1) = 5 + 1 + 2 + 5 + nxm 
    d(X1) = nxx 
    """

    ### various options ### 
    primitive_R0 = false # if true, does R0 as a single depth-1 layer across the whole system 
    numrs = 1 # EC gadgets appear in all places as EC^{numrs}, where EC = Rv^2 Rh^2)
    short = false # if true, replaces EC by Rv Rh 
    no_vert = false # if true, gets rid of Rv (for debugging purposes)
    no_horz = false # same for Rh 
    #######################

    println("building gate: $gate")
    n = 3 # block size 
    gates = []
    tsirelson_swap_gates, swap_depth = permute_circuit(n,false)

    Lx, Ly = get_gate_size(l,"tc",gate)

    local_times = ones(Int,Lx,Ly) # each primitive gate is actually identified with a collection of *vertices* --- thus the local_times array need only keep track of times on vertices, rather than on links 

    # with non-communicating gadgets we should never actually need to worry about the bconds as long as gadgets are property centered 
    function indx(x) 
        return (x + Lx - 1)%Lx + 1
    end 
    function indy(x)
        return (x + Ly - 1)%Ly + 1
    end 

    # need to at least have I and T to be the same depth. Ms always applied in a way which blankets the system, so their depth can be distinct 
    nxx = 5 # depth of T  

    function idn(xraw,yraw,k) # does an identity gate on a single site 
        x = indx(xraw); y = indy(yraw)
        km = k-1
        if k > 0 
            sb = n^(k-1) # small block size 
            bb = n^k # large block size 
            for repind in 1:numrs 
                rn(x,y,km)
            end 
            for j in 0:2 for i in 0:2 
                for repind in 1:nxx # depth of t 
                    idn(x + i*sb,y + j*sb,km) # while the T gate is happening
                end 
            end end 
        else
            push!(gates,[1 0 x y local_times[x,y]]) 
            local_times[x,y] += 1 
        end 
    end 

    function tn(xraw,yraw,o,k) # acts on 2x1 (if o == 1) or 1x2 (if o == 2) block. 
        # T0 updates the local_time of two vertices 
        x = indx(xraw); y = indy(yraw)
        xf = o == 1 ? 1 : 0; yf = 1-xf # determined by the direction the t gate is pointing 
        op = o == 1 ? 2 : 1 # opposite orientation 

        if k > 0 
            km = k-1
            bb = n^k; sb = n^(k-1) # big and small block sizes 

            # first stack of Rs 
            for repind in 1:numrs
                rn(x,y,km); rn(x+xf*bb,y+yf*bb,km)
            end 

            # second layer: (T in center)
            for i in 0:2 
                xp = x + yf*i*sb; yp = y + xf*i*sb 
                idn(xp,yp,km); idn(xp + xf*sb,yp + yf*sb,km) 
                tn(xp + 2sb*xf,yp + 2sb*yf,o,km) # T in center 
                idn(xp + xf*4sb,yp + yf*4sb,km); idn(xp + xf*5sb,yp + yf*5sb,km) 
            end 

            for i in 0:2 
                xp = x + yf*i*sb; yp = y + xf*i*sb 
                # third layer: (two Ts)
                idn(xp,yp,km); tn(xp + xf*sb,yp + yf*sb,o,km) 
                tn(xp + xf*3sb,yp + yf*3sb,o,km); idn(xp + xf*5sb,yp + yf*5sb,km) 
            end 

            for i in 0:2 
                xp = x + yf*i*sb; yp = y + xf*i*sb 
                # fourth layer: (all Ts)
                tn(xp,yp,o,km) 
                tn(xp + xf*2sb,yp + yf*2sb,o,km)
                tn(xp + xf*4sb,yp + yf*4sb,o,km) 
            end 

            for i in 0:2 
                xp = x + yf*i*sb; yp = y + xf*i*sb 
                # fifth layer: same as third 
                idn(xp,yp,km); tn(xp + xf*sb,yp + yf*sb,o,km) 
                tn(xp + xf*3sb,yp + yf*3sb,o,km); idn(xp + xf*5sb,yp + yf*5sb,km) 
            end 

            for i in 0:2 
                xp = x + yf*i*sb; yp = y + xf*i*sb 
                # sixth layer: same as second
                idn(xp,yp,km); idn(xp + xf*sb,yp + yf*sb,km) 
                tn(xp + 2sb*xf,yp + 2sb*yf,o,km) # T in center 
                idn(xp + xf*4sb,yp + yf*4sb,km); idn(xp + xf*5sb,yp + yf*5sb,km) 
            end 

        elseif k == 0 # 
            push!(gates,[2 o x y local_times[x,y]])
            local_times[indx(x),indy(y)] += 1
            local_times[indx(x + xf),indy(y + yf)] += 1
        end 
    end     

    function mn(xraw,yraw,o,k) # "match" gate, acts on 3x1 block in the x (o == 1) or y (o == 2) direction
        x = indx(xraw); y = indy(yraw)
        xf = o == 1 ? 1 : 0; yf = 1-xf 
        op = o == 1 ? 2 : 1 

        if k > 0 
            km = k-1
            bb = n^k; sb = n^(k-1) # big and small block sizes. acts on three big blocks 

            # apply first array of rs 
            for d in 1:numrs
                rn(x,y,km) 
                rn(x+bb*xf,y+bb*yf,km) 
                rn(x+2bb*xf,y+2bb*yf,km) 
            end 

            function apply_tsirelson_swaps() # apply the swap gates that permute the tensor factors of the inputs in the tsirelson way 
                lasttime = tsirelson_swap_gates[1][end]
                for thisgate in tsirelson_swap_gates
                    gtype = thisgate[2]; gloc = thisgate[1]; gtime = thisgate[3]

                    if thisgate[3] < swap_depth # the swap gates contain a redundant layer of identity gates at the end which we want to ignore 
                        for δ in 0:2 
                            if o == 1 # the M_k gate we are applying is oriented in the x direction 
                                gtype == 1 ? idn(x+(gloc-1)*sb,y+δ*sb,km) : tn(x+(gloc-1)*sb,y+δ*sb,1,km)
                            else # it is oriented along y 
                                gtype == 1 ? idn(x+δ*sb,y+(gloc-1)*sb,km) : tn(x+δ*sb,y+(gloc-1)*sb,2,km)
                            end 
                        end 
                    end 
                    lasttime = gtime
                end 
            end 

            # apply first set of swaps 
            apply_tsirelson_swaps()

            if o == 1 # x-oriented M gates (apply x-oriented M_{k-1}s and y-oriented Ts in between the tsirelson swap pattern)
                for k in 0:2 for j in 0:2 mn(x+k*bb,y+j*sb,1,km) end end 
            else # y-oriented gates (applying M_{k-1} in the y direction after swapping)
                for k in 0:2 for j in 0:2 mn(x+j*sb,y+k*bb,2,km) end end 
            end     

            # undo the swap circuit 
            apply_tsirelson_swaps()

        else # level 0 M --- just acts on a 1x3 or 3x1 block (despite the fact that these blocks will always be tiled)
            push!(gates,[3 o x y local_times[x,y]]) # primitive M0 acts on 3 sites, just to make things notationally consistent with higher layers -- x y labels the leftmost / bottommost point 
            # update local times on the 3 vertices in question 
            local_times[indx(x),indy(y)] += 1
            local_times[indx(x+xf),indy(y+yf)] += 1
            local_times[indx(x+2*xf),indy(y+2yf)] += 1 
        end 
    end 

    function rn(xraw,yraw,k) # EC gadget 
        x = indx(xraw); y = indy(yraw)
        
        sb = n^k # so sb = 1 is lattice spacing at smallest scale 

        # if k > 0 && just_R0s # cover the latice with R0s 
        #     for repind in 1:numrs 
        #         for i in 0:2 for j in 0:2 
        #             rn(x+sb*i,y+sb*j,k-1)
        #         end end 
        #     end 
        # else 
            if ~primitive_R0 || k > 0 
                function horz_Ms()
                    for j in 0:2 mn(x,y+j*sb,1,k) end 
                end 
                function vert_Ms()
                    for j in 0:2 mn(x+j*sb,y,2,k) end 
                end 
                function top_Ts()
                    for j in 0:2 tn(x+j*sb,y+sb,2,k); idn(x+j*sb,y,k) end 
                end 
                function bottom_Ts()
                    for j in 0:2 tn(x+j*sb,y,2,k); idn(x+j*sb,y+2sb,k) end 
                end 
                function right_Ts()
                    for j in 0:2 tn(x+sb,y+j*sb,1,k); idn(x,y+j*sb,k) end 
                end 
                function left_Ts()
                    for j in 0:2 tn(x,y+j*sb,1,k); idn(x+2sb,y+j*sb,k) end 
                end 

                function top_Tsp() top_Ts(); horz_Ms() end 
                function bottom_Tsp() bottom_Ts(); horz_Ms() end 
                function left_Tsp() left_Ts(); vert_Ms() end 
                function right_Tsp() right_Ts(); vert_Ms() end 
                
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

                if ~no_vert   
                    vert_part()
                    if ~short vert_part() end 
                end 
                if ~no_horz
                    horz_part()
                    if ~short horz_part() end 
                end 
            
            else # primitive majority (covers a 3x3 square)
                push!(gates,[4 0 x y local_times[x,y]])
                for i in 1:3 for j in 1:3 
                    local_times[indx(x+(i-1)),indy(y+(j-1))] += 1 
                end end 
            end 
        # end 
    end     
    
    # apply the chosen gates: 
    if gate == "I1"
        idn(1,1,l)
    elseif gate ∈ ["IR" "i"]
        idn(1,1,l); rn(1,1,l-1)
    elseif gate ∈ ["TxR" "tx"] # x-oriented T (leading and trailing ECs)
        tn(1,1,1,l); rn(1,1,l-1); rn(1+n^l,1,l-1)
    elseif gate ∈ ["TyR" "ty"] # y-oriented 
        tn(1,1,2,l); rn(1,1,l-1); rn(1,1+n^l,l-1)
    elseif gate ∈ ["MxR" "mx"] # x-oriented M (leading and trailing ECs)
        mn(1,1,1,l); rn(1,1,l-1); rn(1+n^l,1,l-1); rn(1+2*n^l,1,l-1) 
    elseif gate ∈ ["MyR" "my"]
        mn(1,1,2,l); rn(1,1,l-1); rn(1,1+n^l,l-1); rn(1,1+2*n^l,l-1) 
    elseif gate == "R" 
        rn(1,1,l) 
    elseif gate ∈ ["Rsq" "RR"]
        rn(1,1,l); rn(1,1,l)
    elseif gate == "just_R0s"
        nlayers = n^l
        for i in 0:nlayers-1 for j in 0:nlayers-1 
            rn(1+3i,1+3j,0)
        end end 
    elseif gate == "just_R0ssq"
        for repind in 0:1  
            for i in 0:2 for j in 0:2 
                rn(1+3i,1+3j,0)
            end end 
        end 
    else
        println("gate type ($gate) not supported yet")
        return 
    end 

    ngates = length(gates)
    depth = maximum(local_times) - 1 

    if all_boundaries # go through and add on the lvl-0 gadgets on the N and E boundaries if desired 
        boundary_gates = []
        for gate in gates 
            gtype = gate[1]; o = gate[2]; x = gate[3]; y = gate[4]; t = gate[5]
            if y == 1 # copy gates from S boundary to N boundary 
                if gtype == 1 || (gtype == 3 && o == 1) || (gtype == 2 && o == 1)
                    push!(boundary_gates,[gtype o x Ly+1 t])
                end 
            end 
            if x == 1 # copy gates from W boundary to E boundary 
                if gtype == 1 || (gtype == 3 && o == 2) || (gtype == 2 && o == 2)
                    push!(boundary_gates,[gtype o Lx+1 y t])
                end 
            end 
        end 
        for bgate in boundary_gates push!(gates,bgate) end # add boundary gates to big list 
        ngates = length(gates)
    end 

    # turn the gates into an array (rather than a list) sorted in time order --- and check that local_times gives a consistent flat top
    gates_array = zeros(Int,ngates,5) # gtype x y orientation time 
    for i in 1:ngates gates_array[i,:] = gates[i][:] end 
    sort_inds = sortperm(gates_array[:,end])
    gates_array = gates_array[sort_inds,:]
    depth = maximum(local_times) - 1
    # println("size(local times) = $(size(local_times))")
    # println("final local times = ",local_times .- 1)

    @assert depth == minimum(local_times)-1 # need all local times to be equal at the end of the gate construction 

    return gates_array, depth, Lx, Ly  
end 

function master_gate_builder(model,l,gate;all_boundaries=false)
    """ 
    master function for building gates -- just helps streamline code slightly 
    all_boundaries is only passed when building tc gates 
    """
    gates = []; depth = 1; Lx = 1; Ly = 1 
    if model ∈ ["rep" "k2" "rep_5bit"] 
        gates, depth, Lx = tsirelson_gate_builder(model,l,gate); Ly = Lx 
    elseif model == "twod_rep"
        gates, depth, Lx, Ly = twod_tsirelson_gate_builder(l,gate)
    elseif model == "tc"
        gates, depth, Lx, Ly = tc_gate_builder(l,gate,all_boundaries)
    end 

    return gates, depth, Lx, Ly  
end 

function main() 
    model = "tc" # ∈ {"rep" "rep_5bit" "k2" "tc" "twod_rep"}
    gate = "R" # see individual gate builder functions above for options 

    l = 0
    n = model ∈ ["rep_5bit" "k2"] ? 5 : 3 # block size 

    adj = ""
    fout = "data/$(model)$(l)_gates_$(gate)$(adj).jld2" 

    gates = []; depth = 0; L = 0 
    gates, depth, Lx, Ly = master_gate_builder(model,l,gate)

    println("writing gates to file: $fout")
    f = jldopen(fout,"w")
    write(f,"model",model); write(f,"gates",gates); write(f,"gate",gate); write(f,"n",n); write(f,"l",l)
    write(f,"depth",depth); write(f,"Lx",Lx); write(f,"Ly",Ly); write(f,"L",Lx)
    close(f)
end 

if abspath(PROGRAM_FILE) == @__FILE__ 
    main()
end
