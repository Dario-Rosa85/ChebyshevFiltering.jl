function eigen_cheb(hamiltonian_matrix, log_path, log_file_name, lambda_min::Float64=-0.1, lambda_max::Float64=0.1, search_target_ratio::Float64=3.0, max_degree_KPM::Int64=150, stochastic_dimension_KPM::Int64=30, epsilon_convergence::Float64=10^(-7))
    renormalization_hamiltonian!(hamiltonian_matrix)
    e_max = 1.
    e_min = -1.

    search_vector_numbers, polynomial_degree_optim = optimal_search_degree(hamiltonian_matrix, lambda_min, lambda_max, search_target_ratio, max_degree_KPM, stochastic_dimension_KPM, e_max, e_min, N_0)
    search_vector_numbers = Int64(ceil(search_vector_numbers/Threads.nthreads()) * Threads.nthreads())

    chebyshev_coeff = expansion_coeff(polynomial_degree_optim, e_max, e_min, lambda_max, lambda_min)
    damping_coeff = Lanczos_coeff(polynomial_degree_optim, 2)
    full_coeff = Array{Float64}(undef, length(damping_coeff))
    @inbounds for i in eachindex(full_coeff) 
        full_coeff[i] = chebyshev_coeff[i] * damping_coeff[i]
    end

    converged_target_values, converged_target_vectors = polynomial_filtering(search_vector_numbers, polynomial_degree_optim, full_coeff, hamiltonian_matrix, lambda_max, lambda_min, epsilon_convergence, log_path, log_file_name)

    return converged_target_values, converged_target_vectors

end