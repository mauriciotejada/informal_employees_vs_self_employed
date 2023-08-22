"""
    discrete_probability(x,DistFun)

Description: Discretize the distribution
"""
function discrete_probability(x,DistFun)

    N = length(x)
    xmin = x[1]
    xmax = x[N]
    
    TruncDistFun = truncated(DistFun, xmin, xmax)
    GT(x) = cdf(TruncDistFun, x)
    
    ProbG = zeros(N)
    ProbG[1] = GT((x[1]+x[2])/2)
    
    for i = 2:N-1
        mp1 = (x[i-1]+x[i])/2
        mp2 = (x[i]+x[i+1])/2
        ProbG[i] = GT(mp2)-GT(mp1)
    end

    ProbG[N] = 1 - GT((x[N-1]+x[N])/2)
    
    return ProbG
end