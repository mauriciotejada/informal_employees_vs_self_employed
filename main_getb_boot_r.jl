include("InformalEmpVsSelfEmp.jl")

# Necessary Packages
using .InformalEmpVsSelfEmp

using DelimitedFiles
using Statistics
using Roots: fzero

estparam_boot = readdlm("res/est_param_r_boot_0.txt")
Nboot = size(estparam_boot, 2)

estparam = Array{Float64}(undef,14, Nboot)
estparam[1,:] = estparam_boot[1,:]
estparam[2,:] = estparam_boot[2,:]
estparam[3,:] = estparam_boot[1,:]
estparam[4,:] = estparam_boot[3,:]
estparam[5,:] = estparam_boot[4,:]
estparam[6,:] = estparam_boot[3,:]
estparam[7:10,:] = estparam_boot[5:8,:]
estparam[11,:] = estparam_boot[5,:]
estparam[12,:] = estparam_boot[6,:]
estparam[13:14,:] = estparam_boot[9:10,:]

param_boot = Array{Float64}(undef, 2, Nboot)

println("Estimating bootstraped standard errors, be patient!")

Threads.@threads for i in 1:Nboot

    println("Boostrap iteration: $i")

    param = estparam[:,i]

    # getting b
    gparam = [0.16, 0.1]

    θ = gparam[1]*param[13]

    bsol = fzero(x -> getting_b(x, [param; θ], gparam), 0.0)
    

    param_boot[1,i] = bsol
    param_boot[2,i] = θ
    
end

stderror = std(param_boot,dims=2)
writedlm("res/stderror_boot_b_r_0.txt", stderror)
writedlm("res/est_param_boot_b_r_0.txt", param_boot)