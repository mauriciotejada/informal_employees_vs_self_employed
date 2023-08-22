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
    db_sk_lev[zero_dur,2].= 0.25

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

    mus, sigs = initial_param_log_normal(mean([ws; wi]), std([ws; wi])^2)
    muf, sigf = initial_param_log_normal(mean(wf), std(wf)^2)

    if estimate 
        x0 = trans_inv_r([
            1/mean(tu),
            1/mean(tu),
            1/mean([ts; ti]),
            1/mean(tf),
            mus, 
            sigs,
            muf, 
            sigf,
            mean([ws_min, wf_min, wi_min]),
            0.5])

        res_est = optimize(Θ ->  likefun_r(Θ, db_sk_lev), x0, 
                                            NelderMead(),
                                            Optim.Options(show_trace = true, 
                                                            g_tol = 1e-3,
                                                            iterations = 15000))

        writedlm("res/est_param_r_"*string(Int(sk_lev))*".txt", trans_r(res_est.minimizer))
        writedlm("res/likfun_r_"*string(Int(sk_lev))*".txt", -res_est.minimum)

        estparam = Array{Float64}(undef,14)
        estparam[1] = trans_r(res_est.minimizer)[1]
        estparam[2] = trans_r(res_est.minimizer)[2]
        estparam[3] = trans_r(res_est.minimizer)[1]
        estparam[4] = trans_r(res_est.minimizer)[3]
        estparam[5] = trans_r(res_est.minimizer)[4]
        estparam[6] = trans_r(res_est.minimizer)[3]
        estparam[7:10] = trans_r(res_est.minimizer)[5:8]
        estparam[11] = trans_r(res_est.minimizer)[5]
        estparam[12] = trans_r(res_est.minimizer)[6]
        estparam[13:14] = trans_r(res_est.minimizer)[9:10]
        
    else

        estparam_b = readdlm("res/est_param_r_0.txt")

        estparam = Array{Float64}(undef,14)
        estparam[1] = estparam_b[1]
        estparam[2] = estparam_b[2]
        estparam[3] = estparam_b[1]
        estparam[4] = estparam_b[3]
        estparam[5] = estparam_b[4]
        estparam[6] = estparam_b[3]
        estparam[7:10] = estparam_b[5:8]
        estparam[11] = estparam_b[5]
        estparam[12] = estparam_b[6]
        estparam[13:14] = estparam_b[9:10]

    end


    # getting b
    gparam = [0.16, 0.1]

    θ = gparam[1]*estparam[13]

    bsol = fzero(x -> getting_b(x, [estparam; θ], gparam), 0.0)

    # Implied values
    impv = implied_values(estparam)
    writedlm("res/imp_values_r_"*string(Int(sk_lev))*".txt", [impv; bsol; θ])

end

# Unkilled
main(0.0, estimate = true)
