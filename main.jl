include("InformalEmpVsSelfEmp.jl")

# Necessary Packages
using .InformalEmpVsSelfEmp
using BenchmarkTools

using DelimitedFiles
using Statistics
using LinearAlgebra
using Optim
using Roots: fzero

function main(sk_lev; estimate = true)

    # Reading data
    db = readdlm("data/db_geih_2016.csv", ',')

    Isk = Int64.(db[:,4] .== sk_lev)
    db_sk_lev = selif(db, Isk)

    zero_dur = db_sk_lev[:,2].==0
    db_sk_lev[zero_dur,2].= 0.25 # we impute a week for montly duration equal to zero. There are only 695 of these obs.

    Ju = Int64.(db_sk_lev[:,1] .== 0.0)
    Js = Int64.(db_sk_lev[:,1] .== 1.0)
    Jf = Int64.(db_sk_lev[:,1] .== 2.0)
    Ji = Int64.(db_sk_lev[:,1] .== 3.0)

    tu  = selif(db_sk_lev[:,2], Ju)
    ts  = selif(db_sk_lev[:,2], Js)
    tf  = selif(db_sk_lev[:,2], Jf)
    ti  = selif(db_sk_lev[:,2], Ji)

    ws  = selif(db_sk_lev[:,3], Js)
    wf  = selif(db_sk_lev[:,3], Jf)
    wi  = selif(db_sk_lev[:,3], Ji)

    ws_min = minimum(ws)
    wf_min = minimum(wf)
    wi_min = minimum(wi)

    mus, sigs = initial_param_log_normal(mean(ws), std(ws)^2)
    muf, sigf = initial_param_log_normal(mean(wf), std(wf)^2)
    mui, sigi = initial_param_log_normal(mean(wi), std(wi)^2)

    if estimate

        x0 = trans_inv([
            1/mean(tu),
            1/mean(tu),
            1/mean(tu),
            1/mean(ts),
            1/mean(tf),
            1/mean(ti),
            mus, 
            sigs,
            muf, 
            sigf,
            mui, 
            sigi,
            mean([ws_min, wf_min, wi_min]),
            0.5])

        res_est = optimize(Θ ->  likefun(Θ, db_sk_lev), x0, 
                                            NelderMead(),
                                            Optim.Options(show_trace = true, 
                                                            g_tol = 1e-3,
                                                            iterations = 15000))

        writedlm("res/est_param_"*string(Int(sk_lev))*".txt", trans(res_est.minimizer))
        writedlm("res/likfun_"*string(Int(sk_lev))*".txt", -res_est.minimum)
        estparam = trans(res_est.minimizer)

    else

        estparam = readdlm("res/est_param_0.txt")

    end

    # getting b
    gparam = [0.16, 0.1]

    θ = gparam[1]*estparam[13]

    bsol = fzero(x -> getting_b(x, [estparam; θ], gparam), 0.0)

    # Implied values
    impv = implied_values(estparam)
    writedlm("res/imp_values_"*string(Int(sk_lev))*".txt", [impv; bsol; θ])

end

# Unkilled
main(0.0, estimate = true)

