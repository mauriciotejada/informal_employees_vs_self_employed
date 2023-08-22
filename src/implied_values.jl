function implied_values(param)

    imv = Array{Float64}(undef,10)

    λs  = param[1]
    λf  = param[2]
    λi  = param[3]
    ηs  = param[4]
    ηf  = param[5]
    ηi  = param[6]
    μs  = param[7]
    σs  = param[8]
    μf  = param[9]
    σf  = param[10]
    μi  = param[11]
    σi  = param[12]
    ρU  = param[13]

    # Hazard rates
    Gs_dist = LogNormal(μs,σs)
    Gf_dist = LogNormal(μf,σf)
    Gi_dist = LogNormal(μi,σi)

    # Hazards
    hus =  λs*(1-cdf(Gs_dist, ρU))
    huf =  λf*(1-cdf(Gf_dist, ρU))
    hui =  λi*(1-cdf(Gi_dist, ρU))
    imv[1] = (hus + huf + hui)
    imv[2] = ηs
    imv[3] = ηf
    imv[4] = ηi

    # Wage offers Distributions
    avewo(mu,sig) = exp(mu + 0.5*sig^2)
    sigwo(mu,sig) = exp(2*mu + sig^2)*(exp(sig^2)-1)

    imv[5] = avewo(μs,σs)
    imv[6] = sigwo(μs,σs)
    imv[7] = avewo(μf,σf)
    imv[8] = sigwo(μf,σf)
    imv[9] = avewo(μi,σi)
    imv[10] = sigwo(μi,σi) 

    return imv
end
