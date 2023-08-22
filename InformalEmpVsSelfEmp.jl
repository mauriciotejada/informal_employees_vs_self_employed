module InformalEmpVsSelfEmp

    # Dependencies
    using Parameters
    using Distributions: Distribution, LogNormal, Exponential, truncated, cdf, pdf, rand
    using Statistics
    using QuadGK
    using DelimitedFiles
    using LinearAlgebra
    using Random
    using LaTeXStrings

    # Export types
    export Model

    # Export functions
    export discrete_probability, solvemodel, getting_b
    
    export likefun, likefun_r

    export tab_table, selif, trans, trans_inv, trans_r, trans_inv_r, initial_param_log_normal

    # Contains

    # Types
    include("src/model_type.jl")

    # Functions
    include("src/probability_discretization.jl")
    include("src/solve_model.jl")
    include("src/likelihood_function.jl")
    include("src/implied_values.jl")
    include("src/getting_b.jl")

    # Utilities
    include("src/utils.jl")
    include("src/trans.jl")
    include("src/initial_param_log_normal.jl")

end