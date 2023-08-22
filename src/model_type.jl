@with_kw struct Model

    # Dynamics
    λs::Float64 = 0.3
    λf::Float64 = 0.2
    λi::Float64 = 0.05
    ηs::Float64 = 0.01
    ηf::Float64 = 0.05
    ηi::Float64 = 0.05

    # Earnings distributions
    μs::Float64 = 0.1
    σs::Float64 = 0.5
    μf::Float64 = 0.5
    σf::Float64 = 0.3
    μi::Float64 = -0.1
    σi::Float64 = 0.9

    # Other parameters
    b::Float64 = -0.1
    θ::Float64 = 0.5
    τ::Float64 = 0.1
    ρ::Float64 = 0.1

    # Grid points
    Nx::Int64 = 500
    xmin::Float64 = 1e-10
    xmax::Float64 = 10
end
