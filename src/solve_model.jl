"""
    solvemodel(mod::Model; 
            Δ = 0.25,
            tol = 1e-6,
            max_iter = 500,
            flag_print_iter = true)

TBW
"""
function solvemodel(mod::Model; 
                    Δ = 0.25,
                    tol = 1e-6,
                    max_iter = 500,
                    flag_print_iter = true)

    # Initial variable setup
    count_iter = 0
    diff = 1    

    # Unpacking parameters
    @unpack λs, λf, λi, ηs, ηf, ηi, μs, σs, μf, σf, μi, σi, b, θ, τ, ρ, Nx, xmin, xmax = mod

    # Grids
    x = collect(range(xmin, xmax, length = Nx))

    # Earningns distributions
    Gs_dist = LogNormal(μs,σs)
    Gf_dist = LogNormal(μf,σf)
    Gi_dist = LogNormal(μi,σi)

    Gs = discrete_probability(x, Gs_dist)
    Gf = discrete_probability(x, Gf_dist)
    Gi = discrete_probability(x, Gi_dist)

    # Initial guess for values
    U₀ = b/ρ
    Es₀ = x/ρ
    Ef₀ = (x*(1-τ) .+ θ)/ρ
    Ei₀ = x/ρ

    while diff > tol

        U₁ = ( b + λs*sum(max.(Es₀.-U₀,0).*Gs) + λf*sum(max.(Ef₀.-U₀,0).*Gf) + 
            λi*sum(max.(Ei₀.-U₀,0).*Gi) ) / ρ

        Es₁ = ( x + ηs*(U₀ .- Es₀) ) / ρ 

        Ef₁ = ( x*(1-τ) .+ θ + ηf*(U₀ .- Ef₀) ) / ρ  

        Ei₁ = ( x + ηi*(U₀ .- Ei₀) ) / ρ 

        # Difference
        diff = maximum([abs.(U₁-U₀); abs.(Es₁[:]-Es₀[:]); abs.(Ef₁[:]-Ef₀[:]); abs.(Ei₁[:]-Ei₀[:])])

        if flag_print_iter
            println(diff)
        end

        # Updating
        U₀  = Δ * U₁ +  (1 - Δ) * U₀
        Es₀ = Δ * Es₁ +  (1 - Δ) * Es₀
        Ef₀ = Δ * Ef₁ + (1 - Δ) * Ef₀
        Ei₀ = Δ * Ei₁ + (1 - Δ) * Ei₀

        count_iter += 1

        if count_iter > max_iter
            if flag_print_iter
                println("The problem did not converged, max number of iterations reached")
            end
            flag_convergence = 0
            return (res_values = NaN, hazard_rates = NaN, states = NaN, U = NaN, 
            Es = NaN, Ef = NaN, Ei = NaN, conv = flag_convergence, iterations = count_iter)
        end

    end

    # Reserations values
    xss = max(ρ*U₀,xmin)
    xfs = max((ρ*U₀ - θ) / (1-τ),xmin)
    xis = max(ρ*U₀,xmin)

    hus = λs*(1-cdf(Gs_dist,xss))
    huf = λf*(1-cdf(Gf_dist,xfs))
    hui = λi*(1-cdf(Gi_dist,xis))
    hu = hus + huf + hui
    hs = ηs
    hf = ηf
    hi = ηi

    HAZ = [ hus -hs 0.0 0.0;
            huf 0.0 -hf 0.0;
            hui 0.0 0.0 -hi;
            1 1 1 1]

    states = HAZ\[0;0;0;1]     
    hazard_rates = [hu; hs; hf; hi]
    res_values = [xss; xfs; xis]
    flag_convergence = 1

    return return (res_values = res_values, hazard_rates = hazard_rates, states = states, U = U₀, 
                   Es = Es₀, Ef = Ef₀, Ei = Ei₀, conv = flag_convergence, iterations = count_iter)

end
