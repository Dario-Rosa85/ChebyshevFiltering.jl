function polynomial_filtering(search_vector_numbers, polynomial_degree_optim, full_coeff, hamiltonian_matrix, lambda_max, lambda_min, epsilon_convergence, log_path, log_file_name)
    ########Random vectors inizialization
    search_vectors_list = orthogonalize_QR(randn(ComplexF64, size(hamiltonian_matrix, 2), search_vector_numbers))
    Ritz_matrix = Matrix{ComplexF64}(undef, size(search_vectors_list, 2), size(search_vectors_list, 2))
    provisional_vector = map(x -> zeros(ComplexF64, size(hamiltonian_matrix, 2)), 1:Threads.nthreads())
    u_vectors = map(x -> zeros(ComplexF64, size(hamiltonian_matrix, 2)), 1:search_vector_numbers)
    w_vectors = map(x -> zeros(ComplexF64, size(hamiltonian_matrix, 2)), 1:search_vector_numbers)
    Ritz_vectors = map(x -> zeros(ComplexF64, size(hamiltonian_matrix, 2)), 1:search_vector_numbers) 
    
    convergence_reached = "false"
    number_of_iterations = 0
    while convergence_reached == "false"
        ##########################################################################################################
        ########Filtering step and orthogonalization
        filtering_step!(search_vectors_list, u_vectors, w_vectors, provisional_vector, polynomial_degree_optim, full_coeff, hamiltonian_matrix, lambda_max, lambda_min)
        search_vectors_list = orthogonalize_QR(search_vectors_list)
        ##########################################################################################################
        ##########################################################################################################
        ########Ritz pairs and residuals computation
        Ritz_values, residuals = Rayleigh_Ritz_pairs_residuals!(Ritz_matrix, Ritz_vectors, search_vectors_list, provisional_vector, hamiltonian_matrix)
        ##########################################################################################################
        ##########################################################################################################
        ########Convergence test
        converged_target_values, converged_target_vectors, not_converged_residuals = convergence_test(Ritz_values, Ritz_vectors, residuals, lambda_min, lambda_max, epsilon_convergence)
        if length(not_converged_residuals) == 0
            convergence_reached = "true"
            return converged_target_values, reduce(hcat, converged_target_vectors)
        else
            smallest_not_converged_residual = first(sort(not_converged_residuals))
            if ((epsilon_convergence < smallest_not_converged_residual < epsilon_convergence^(1/2)) || (number_of_iterations < 4) || (length(not_converged_residuals) > 10))
                convergence_reached = "false"
                number_of_iterations += 1
                open(log_path * log_file_name, "a") do io
                    println(io, "convergence not yet reached. Time is ", Dates.now())
                end
            else
                convergence_reached = "true"
                return converged_target_values, converged_target_vectors
            end
        end        
    end
end