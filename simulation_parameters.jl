function parameter_repository(model,mode,l,nps,logdist,bias,gadgetnoise,measurementnoise,test)
    """ 
    this is a place to store numerical information about the locations of thresholds in different types of models + conventions on number of samples to use + etc. 

    model: name of code ∈ ["rep" "k2" "tc" "twod_rep"]
    l: log_n of system size 
    nps: desired number of noise steps 
    logdist: returns linearly-spaced noise values if false; log-spaced values if true 
    bias: bias of iid noise (not relevant for TC)
    gadgetnoise: if true, gadget-based noise 
    measurementnoise: if true, noise occurs only in syndrome measurements 
    test: if true, does a reduced number of samples / measurement time to speed things up 

    returns: (vector of noise values, pc, number of samples)
    """
    pmin = 0; pmax = 0 
    ps = []; pc = 0 
    samps = 1 # O(1) if mode == stats, large for Ft, slightly smaller possible for trel 
    periods = 2 # number of applications of one floquet period for the total time evolution (not needed if mode == "trel")
    thermal_periods = 0 # number of periods to thermalize for before measuring properties of the steady state 

    ####### number of periods #######  
    if mode == "Ft" # this number does not really matter -- just sets how close to 1 the value of F at the crossing point is (closer for smaller periods)
        periods = 50 
        if model == "rep_5bit" periods = 5 end 
        if model == "tc" periods = 10 end 
    end 
    if mode == "stats"
        thermal_periods = 10000 #
        if l == 0 periods = 5000000 # poopers 
        elseif l == 1 periods = 4000000 
        elseif l == 2 periods = 9000000  
        elseif l == 3 periods = 800000 
        elseif l == 4 periods = 90000 
        end 
        
        if model ∈ ["k2" "rep_5bit"] # both have block size 5 -- use less periods 
            if l == 1 periods = 200000; thermal_periods = 20000 # 500k periods is very easily enough to get good statistics 
            elseif l == 2 periods = 25000; thermal_periods = 14000 
            elseif l == 3 periods = 5000; thermal_periods = 1200 
            elseif l == 4 periods = 600; thermal_periods = 200
            else periods = 200; thermal_periods = 20 end 
        end 
    end 
    if mode == "hist" periods = 1 end 

    ####### noise strengths and samples: 1d codes ####### 

    ### 3-bit rep code ### 
    if model == "rep"
        pc = bias == 0 ? .0139 : .0063 # .0063 is for maximal bias 
        if gadgetnoise 
            pc = .008 # when each gadget randomizes its output 
            pc = .0142 # when maj gates give the wrong outputs (not fully ergodic since magnetization always a multiple of 3)
        end 
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

            if gadgetnoise # note: gadgetnoise currently randomizes gadget outputs and hence currently does not support a bias 
                pmin = 5e-3; pmax = 11e-3 # all gadgets fail w/ prob p, and a failed gadget randomizes its output 
                pmin = .01; pmax = .016 # only Z gadgets fail, and randomize their outputs 
            end 
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
        elseif mode == "stats" 
            pmin = .001; pmax = .055  
            # sweep values with old code: 
            if l == 1 pmin = .001; pmax = .175
            elseif l == 2 pmin = .001; pmax = .055 
            elseif l == 3 pmin = .003; pmax = .035
            elseif l == 4 pmin = .01; pmax = .021 
            end 

            # rather wide range about critical point (includes maximum of chi, farish away from p_c for moderate L)
            # if l == 1 pmin = .001; pmax = .2 #
            # elseif l == 2 pmin = .001; pmax = .11 
            # elseif l == 3 pmin = .003; pmax = .065
            # elseif l == 4 pmin = .01; pmax = .05 
            # end 

            # narrow range about critical point -- for fixing value of pc 
            # pmin = .01; pmax = .018 

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

    if test # shorten stuff if want to be fast 
        if periods > 2 
            periods = round(Int,periods/10)
        end 
        if samps > 2 
            samps = round(Int,samps/10) 
        end 
    end 

    return ps, pc, samps, thermal_periods, periods 
end 