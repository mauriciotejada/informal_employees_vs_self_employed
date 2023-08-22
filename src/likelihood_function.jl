"""
    likefun(x, gp, data; 
                    write_intermediate_res = false)

TBW
"""
function likefun(x, data; 
                    write_intermediate_res = false)
                    

    # Data 
    Ju = Int64.(data[:,1] .== 0.0)
    Js = Int64.(data[:,1] .== 1.0)
    Jf = Int64.(data[:,1] .== 2.0)
    Ji = Int64.(data[:,1] .== 3.0)

    tu  = selif(data[:,2], Ju)
    ts  = selif(data[:,2], Js)
    tf = selif(data[:,2], Jf)
    ti = selif(data[:,2], Ji)

    ws  = selif(data[:,3], Js)
    wf = selif(data[:,3], Jf)
    wi = selif(data[:,3], Ji)

    # Model Definition

    xt = trans(x)

    if write_intermediate_res 
        writedlm("res/est_param_iter.txt", xt)
    end

    # Estimated parameters
    λs = xt[1]
    λf = xt[2]
    λi = xt[3]
    ηs = xt[4]
    ηf = xt[5]
    ηi = xt[6]
    μs = xt[7]
    σs = xt[8]
    μf = xt[9]
    σf = xt[10]
    μi = xt[11]
    σi = xt[12]
    ρU = xt[13]
    σₑ = xt[14]

    # Distributions
    μₑ = -(σₑ^2)/2.0
    V_dist  = LogNormal(μₑ,σₑ)
    Gs_dist = LogNormal(μs,σs)
    Gf_dist = LogNormal(μf,σf)
    Gi_dist = LogNormal(μi,σi)

    # Hazards
    hus =  λs*(1-cdf(Gs_dist, ρU))
    huf =  λf*(1-cdf(Gf_dist, ρU))
    hui =  λi*(1-cdf(Gi_dist, ρU))
    hu = hus + huf + hui
    hs = ηs
    hf = ηf
    hi = ηi

    HAZARDS = [ hus -hs 0.0 0.0;
                huf 0.0 -hf 0.0;
                hui 0.0 0.0 -hi;
                1 1 1 1]

    st = HAZARDS\[0;0;0;1]     
    haz = [hu; hs; hf; hi]

    # Densities for duration
    d_u_pdf(t) = haz[1]*exp.(-haz[1]*t)
    d_s_pdf(t) = haz[2]*exp.(-haz[2]*t)
    d_f_pdf(t) = haz[3]*exp.(-haz[3]*t)
    d_i_pdf(t) = haz[4]*exp.(-haz[4]*t)

    # Densities for wages
    v(e) = pdf.(V_dist, e)

    x_s_pdf(x) = pdf.(truncated(Gs_dist, lower = ρU), x)
    x_f_pdf(x) = pdf.(truncated(Gf_dist, lower = ρU), x)
    x_i_pdf(x) = pdf.(truncated(Gi_dist, lower = ρU), x)

    x_s_pdf_me(x) = quadgk(y -> (1.0./y).*v(x./y).*x_s_pdf(y), ρU, Inf)[1]
    x_f_pdf_me(x) = quadgk(y -> (1.0./y).*v(x./y).*x_f_pdf(y), ρU, Inf)[1]
    x_i_pdf_me(x) = quadgk(y -> (1.0./y).*v(x./y).*x_i_pdf(y), ρU, Inf)[1]

    # Contributions
    Lu = d_u_pdf(tu)*st[1]
    Ls = d_s_pdf(ts).*x_s_pdf_me(ws).*st[2]
    Lf = d_f_pdf(tf).*x_f_pdf_me(wf).*st[3]
    Li = d_i_pdf(ti).*x_i_pdf_me(wi).*st[4]

    # Likelihood
    L = - ( sum(log.(Lu)) + sum(log.(Ls)) + sum(log.(Lf)) + sum(log.(Li)) )
 
    if isnan(L) == 1
        loglike = 999999
    elseif isinf(L) == 1
        loglike = 999999
    else
        loglike = L
    end

    return loglike
    
end


"""
    likefun_r(x, gp, data; 
                    flag_print_iter = false,
                    write_intermediate_res = false)

TBW
"""
function likefun_r(x, data; 
                    write_intermediate_res = false)

    # Data 
    Ju = Int64.(data[:,1] .== 0.0)
    Js = Int64.(data[:,1] .== 1.0)
    Jf = Int64.(data[:,1] .== 2.0)
    Ji = Int64.(data[:,1] .== 3.0)

    tu  = selif(data[:,2], Ju)
    ts  = selif(data[:,2], Js)
    tf = selif(data[:,2], Jf)
    ti = selif(data[:,2], Ji)

    ws  = selif(data[:,3], Js)
    wf = selif(data[:,3], Jf)
    wi = selif(data[:,3], Ji)

    # Model Definition

    xt = trans_r(x)

    if write_intermediate_res 
        writedlm("res/est_param_iter.txt", xt)
    end

    # Estimated parameters
    λs = xt[1]
    λf = xt[2]
    λi = xt[1]
    ηs = xt[3]
    ηf = xt[4]
    ηi = xt[3]
    μs = xt[5]
    σs = xt[6]
    μf = xt[7]
    σf = xt[8]
    μi = xt[5]
    σi = xt[6]
    ρU  = xt[9]
    σₑ = xt[10]

    # Distributions
    μₑ = -(σₑ^2)/2.0
    V_dist  = LogNormal(μₑ,σₑ)
    Gs_dist = LogNormal(μs,σs)
    Gf_dist = LogNormal(μf,σf)
    Gi_dist = LogNormal(μi,σi)

    # Hazards
    hus =  λs*(1-cdf(Gs_dist, ρU))
    huf =  λf*(1-cdf(Gf_dist, ρU))
    hui =  λi*(1-cdf(Gi_dist, ρU))
    hu = hus + huf + hui
    hs = ηs
    hf = ηf
    hi = ηi

    HAZARDS = [ hus -hs 0.0 0.0;
                huf 0.0 -hf 0.0;
                hui 0.0 0.0 -hi;
                1 1 1 1]

    st = HAZARDS\[0;0;0;1]     
    haz = [hu; hs; hf; hi]

    # Densities for duration
    d_u_pdf(t) = haz[1]*exp.(-haz[1]*t)
    d_s_pdf(t) = haz[2]*exp.(-haz[2]*t)
    d_f_pdf(t) = haz[3]*exp.(-haz[3]*t)
    d_i_pdf(t) = haz[4]*exp.(-haz[4]*t)

    # Densities for wages
    v(e) = pdf.(V_dist, e)

    x_s_pdf(x) = pdf.(truncated(Gs_dist, lower = ρU), x)
    x_f_pdf(x) = pdf.(truncated(Gf_dist, lower = ρU), x)
    x_i_pdf(x) = pdf.(truncated(Gi_dist, lower = ρU), x)

    x_s_pdf_me(x) = quadgk(y -> (1.0./y).*v(x./y).*x_s_pdf(y), ρU, Inf)[1]
    x_f_pdf_me(x) = quadgk(y -> (1.0./y).*v(x./y).*x_f_pdf(y), ρU, Inf)[1]
    x_i_pdf_me(x) = quadgk(y -> (1.0./y).*v(x./y).*x_i_pdf(y), ρU, Inf)[1]

    # Contributions
    Lu = d_u_pdf(tu)*st[1]
    Ls = d_s_pdf(ts).*x_s_pdf_me(ws).*st[2]
    Lf = d_f_pdf(tf).*x_f_pdf_me(wf).*st[3]
    Li = d_i_pdf(ti).*x_i_pdf_me(wi).*st[4]

    # Likelihood
    L = - ( sum(log.(Lu)) + sum(log.(Ls)) + sum(log.(Lf)) + sum(log.(Li)) )

    if isnan(L) == 1
        loglike = 999999
    elseif isinf(L) == 1
        loglike = 999999
    else
        loglike = L
    end

    return loglike
    
end