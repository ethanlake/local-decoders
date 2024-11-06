function generate_subsets(vertices)
    """
    function for enumerating all subsets of a collection of vertices 
    input: collection of x/y coordinates of vertices: [[x1 y1] [x2 y2] ...]
    """

    n = length(vertices)
    subsets = []
    for i in 0:(2^n - 1) # loop over all binary strings of length |vertices| --- a 1 in the jth position of the ith binary string indicates that vertex j will be inlcuded in the ith set. 
        subset = []
        for j in 1:n
            if (i >> (j - 1)) & 1 == 1 # true if the jth digit of the binary string is 1
                push!(subset, vertices[j])
            end
        end
        push!(subsets, subset)
    end
    return subsets
end

function even_vertex_sets(dims)
    """
    define an even vertex set to be a set containing an even number of vertices of a square grid of n x m vertices, where dims = [n m]. this function returns a list containing all even vertex sets of this square grid. 
    """
    # generate a set of all vertices 
    vertices = [(x, y) for x in 0:dims[1]-1, y in 0:dims[2]-1]
    vertices = vec(vertices)  # flatten the array of tuples
    # println(vertices)

    all_subsets = generate_subsets(vertices)

    # small-brained solution: filter out all those subsets with an odd number of vertices 
    even_subsets = [subset for subset in all_subsets if length(subset) % 2 == 0]

    return even_subsets
end

function coarse_grain_damage(synds) 
    """
    given Lx x Ly array of syndromes, with nontrivial syndromes on the 1-frame, returns coordinates of syndromes in coarse-grained lattice 
    """
    Lx, Ly = size(synds) 
    output_synds = []
    for syndloc in findall(x->x == true,synds)
        @assert syndloc[1] % 3 == 1; @assert syndloc[2] % 3 == 1 
        sx = (syndloc[1]-1) ÷ 3; sy = (syndloc[2]-1) ÷ 3
        push!(output_synds,(sx,sy))
    end 
    # @assert length(output_synds) % 2 == 0 
    return output_synds 
end 

function synds_to_links(synds,Lx,Ly)
    """ 
    given a list of syndrome locations, returns a configuration of link spins on an Lx x Ly lattice with syndromes at the appropriate locations. if |synds| is odd, an addition syndrome will be created at (1,1) (this is not a problem since this is a linear point of all gadgets)
    """ 

    spins = zeros(Bool,Lx,Ly,4)
    for syndloc in synds 
        x = syndloc[1]; y = syndloc[2]
        for xp in 1:(x-1) # first segment of error path 
            spins[xp,1,1] ⊻= true 
        end 
        for yp in 1:(y-1) # second segment of error path 
            spins[x,yp,2] ⊻= true 
        end 
    end 
    return spins 
end 

function get_noise_hist(x,y,t,synds,T,Lx,Ly)
    """
    creates a noise history on a system of spacetime volume T x Lx x Ly with a level-0 gadget failure at spacetime location x,y,t with gadget orientation o and failure determined by error (a list of syndrom locations relative to x,y)
    returns: T x Lx x Ly x 2 array of error locations 
    """

    noise_hist = zeros(Bool,T,Lx,Ly,4)
    for syndloc in synds # syndloc a tuple of coordinates --- smallest values are 0 
        dx = syndloc[1]; dy = syndloc[2]
        # the noise is on links -- so create a pattern of link noise with syndromes at the desired locations. convention is that links are flipped along a path which proceeds from the origin of the gadget to the syndrome location along a 2-segment path which first moves along the x direction, and then along y. 
        xmax = x+dx # x value at which path turns 
        if xmax > x # if the path moves along the x direction by a nonzero amount 
            for xp in x:(xmax-1) # first segment of error path 
                noise_hist[t,xp,y,1] ⊻= true 
            end 
        end 
        ymax = y+dy 
        if ymax > y
            for yp in y:(ymax-1) # second segment of error path 
                noise_hist[t,xmax,yp,2] ⊻= true 
            end 
        end 
    end 

    return noise_hist 
end         

function count_error_weights(errors)
    """
    given a collection of errors, returns a list summarizing the number of errors of different weights (number of anyons)
    """

    error_weights = Dict{Int64, Int64}()
    for error in errors 
        weight = length(error)
        if weight ∉ keys(error_weights)
            error_weights[weight] = 1 
        else 
            error_weights[weight] += 1 
        end 
    end 
    weights_list = []
    for weight in keys(error_weights)
        push!(weights_list,[weight error_weights[weight]])
    end 

    return weights_list 
end 
