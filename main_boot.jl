include("InformalEmpVsSelfEmp.jl")

# Necessary Packages
using .InformalEmpVsSelfEmp
using BenchmarkTools

using DelimitedFiles
using Statistics
using LinearAlgebra
using Optim
using Random
using StatsBase

function main(sk_lev)

    # Bootstrap options
    Nboot = 50

    # Random seed
    Random.seed!(170678)

    # Reading data
    db = readdlm("data/db_geih_2016.csv", ',')

    #sk_lev = 0.0
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

    param_boot = Array{Float64}(undef, length(x0), Nboot)

    println("Estimating bootstraped standard errors, be patient!")

    Threads.@threads for i in 1:Nboot

        sam_boot = sample(collect(1:size(db_sk_lev,1)), size(db_sk_lev,1), replace=true)

        db_sk_lev_boot = db_sk_lev[sam_boot, :]

        println("Boostrap iteration: $i")

        res_est = optimize(Θ ->  likefun(Θ, db_sk_lev_boot), x0, 
                                    NelderMead(),
                                    Optim.Options(show_trace = false,  
                                                        g_tol = 1e-3,
                                                        iterations = 15000))        

        param_boot[:,i] = trans(res_est.minimizer)

    end

    stderror = std(param_boot,dims=2)
    writedlm("res/stderror_boot_"*string(Int(sk_lev))*".txt", stderror)
    writedlm("res/est_param_boot_"*string(Int(sk_lev))*".txt", param_boot)

end

# Unkilled
main(0.0)
