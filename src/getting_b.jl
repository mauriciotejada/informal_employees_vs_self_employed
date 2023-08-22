function getting_b(b0, param, gp)

        # Define the Model
        mod = Model(
            # Dynamics
            λs = param[1],
            λf = param[2],
            λi = param[3],
            ηs = param[4],
            ηf = param[5],
            ηi = param[6],
    
            # Earnings distributions
            μs = param[7],
            σs = param[8],
            μf = param[9],
            σf = param[10],
            μi = param[11],
            σi = param[12],
    
            # Other parameters
            b = b0,
            θ = param[15],
            τ = gp[1],
            ρ = gp[2]
        )
    
        # Model solution
        sol = solvemodel(mod, flag_print_iter = false)

        fout = param[13] - gp[2]*sol.U

        return fout

end
