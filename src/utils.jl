# Funtion to compute a frequency table for a vector
function tab_table(x::Vector, levels::Vector; prop=true)

    # levels = unique(x)
	x = sort(x) # sorted input aligns with temp (lowest to highest)
	t = zeros(length(levels)) # vector for freqs
	# frequency for each value
	for i in 1:length(levels)
    	t[i] = sum(x .== levels[i])
	end

    if prop
        return t/sum(t)
    else
        return t
    end

end

# Funtion to compute a two way frequency table for two vectors of the same size
function tab_cross_table(x::Vector, y::Vector, levelsx::Vector, levelsy::Vector; prop = "none")
    
    #levelsx = unique(x)
    #levelsy = unique(y)

	t = zeros(length(levelsx),length(levelsy)) # vector for freqs
	# frequency for each value
	for i in 1:length(levelsx), j in 1:length(levelsy)
        t[i,j] = sum((x .== levelsx[i]) .& (y .== levelsy[j]))
    end

    if prop == "none"
	    return t
    elseif prop == "cols"
        return t./sum(t, dims = 1)
    elseif prop == "rows"
        return t./sum(t, dims = 2)
    else 
        return t./sum(t)
    end
end

# Recode a categorical vector
function recode(x,oldcodes,newcodes)

    # Usage  x_recoded = recode(x, old_codes, new_codes)
    N = length(oldcodes)
    out = copy(x)
    
    for i in 1:N
        out[x.==oldcodes[i]] .= newcodes[i]
    end
    
    return out
end

# Trim rows (borrowed from James P. LeSage)
function trimr(x::Vector, n1, n2)

    n = length(x)

    if (n1+n2) >= n
        error("Attempting to trim too much in trimr")
    end

    h1 = n1 + 1   
    h2 = n - n2
    
    return x[h1:h2]
end

function trimr(x::Matrix, n1, n2)

    n, k = size(x)

    if (n1+n2) >= n
        error("Attempting to trim too much in trimr")
    end

    h1 = n1 + 1   
    h2 = n - n2
    
    return x[h1:h2,:]
end

# Creates a matrix (or vector) of lagged values (the first obs are NaN)
# (borrowed from Ambrogio Cesa-Bianchi)
function L(x::Vector; nlag = 0, init = NaN)

    if nlag > 0 # this is lead
        zt = ones(nlag)*init
        return vcat(trimr(x,nlag,0), zt)

    elseif nlag < 0 # this is lag
        zt = ones(abs(nlag))*init;
        return vcat(zt, trimr(x,0,abs(nlag)))
    else
        return x
    end
    
end

function L(x::Matrix; nlag = 0, init = NaN)

    n, k = size(x)

    if nlag > 0 # this is lead
        zt = ones(nlag, k)*init
        return vcat(trimr(x,nlag,0), zt)

    elseif nlag < 0 # this is lag
        zt = ones(abs(nlag), k)*init;
        return vcat(zt, trimr(x,0,abs(nlag)))
    else
        return x
    end
    
end

# Common sample of a matrix
# If a row of x contains a NaN, the row is removed. If dim=2, if a 
# column of x contains a NaN, the column is removed.
function common_sample(x::Matrix; dim = 1)

    if dim == 1
        nonnan = sum(isnan.(x), dims=2) .== 0
        x_common = x[nonnan[:], :]
        index = 1 .- nonnan
    else
        nonnan = sum(isnan.(x), dims=1) .== 0
        x_common = x[:, nonnan[:]]
        index = 1 .- nonnan
    end

    return x_common, index

end 

# Substitute NaN values with value (default -99)
# (borrowed from Ambrogio Cesa-Bianchi)
function nan_2_num(x::Matrix; valrep = -99)

    xnonan = copy(x)
    xnonan[isnan.(x) .== 1] .= valrep
    return xnonan

end

# Substitute values (default -99) with NaN
# (borrowed from Ambrogio Cesa-Bianchi)
function num_2_nan(x::Matrix; valrep = -99)

    xnan = copy(x)
    xnan[x .== valrep] .= NaN
    return xnan

end

# SELIF(x,t) selects the elements of x for which t=1
# x: NxK matrix, t: Nx1 matrix of 0's and 1's
function selif(x::Vector, t::Vector)
    
    nx = length(x)
    nt = length(t)

    if nx != nt
        error("Size must agree")
    end

    ind = t.==1

    xs  = x[ind,:]

    return xs
end

function selif(x::Matrix, t::Vector)
    
    nx, kx = size(x)
    nt = length(t)

    if nx != nt
        error("Size must agree")
    end

    ind = t.==1

    xs  = x[ind,:]

    return xs
end