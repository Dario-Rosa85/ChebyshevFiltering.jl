function polynomial_filtering(search_vector_numbers, polynomial_degree_optim, full_coeff, hamiltonian_matrix::SparseMatrixCSC{ComplexF64, Int64}, lambda_max, lambda_min, epsilon_convergence, log_path, log_file_name)
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
        filtering_step!(search_vectors_list, u_vectors, w_vectors, provisional_vector, polynomial_degree_optim, full_coeff, hamiltonian_matrix)
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
                if log_path != "none" && log_file_name != "none"
                    open(joinpath(log_path, log_file_name), "a") do io
                        println(io, "Convergence not yet reached. Time is ", Dates.now())
                    end
                end
            else
                convergence_reached = "true"
                if log_path != "none" && log_file_name != "none"
                    open(joinpath(log_path, log_file_name), "a") do io
                        println(io, "Convergence reached. Time is ", Dates.now())
                    end
                end
                return converged_target_values, converged_target_vectors
            end
        end        
    end
end

function filtering_step!(search_vectors_list::Matrix{ComplexF64}, u_vectors::Vector{Vector{ComplexF64}}, w_vectors::Vector{Vector{ComplexF64}}, provisional_vector::Vector{Vector{ComplexF64}}, polynomial_degree_optim, full_coeff, hamiltonian_matrix::SparseMatrixCSC{ComplexF64, Int64})
    length_chunks = collect(Iterators.partition(1:size(search_vectors_list, 2), Int64(floor(size(search_vectors_list, 2) / (Threads.nthreads())))))
    @inbounds Threads.@threads for i in 1:Threads.nthreads() 
        @inbounds for k in length_chunks[i]
            mul!(u_vectors[k], hamiltonian_matrix, search_vectors_list[:, k])
            provisional_vector[Threads.threadid()] = search_vectors_list[:, k]
            mul!(provisional_vector[Threads.threadid()], hamiltonian_matrix, u_vectors[k], 2.0, -1.0)
            w_vectors[k] = provisional_vector[Threads.threadid()]
            search_vectors_list[:, k] .= full_coeff[1] * search_vectors_list[:, k] .+ full_coeff[2] .* u_vectors[k] .+ full_coeff[3] .* w_vectors[k]   
        end
    end
    @inbounds Threads.@threads for i in 1:Threads.nthreads()
        @inbounds for k in length_chunks[i]
            @inbounds for n in 4:polynomial_degree_optim 
                provisional_vector[Threads.threadid()] = u_vectors[k]
                mul!(provisional_vector[Threads.threadid()], hamiltonian_matrix, w_vectors[k], 2.0, -1.0)
                u_vectors[k] = w_vectors[k]
                w_vectors[k] = provisional_vector[Threads.threadid()]
                search_vectors_list[:, k] .+= full_coeff[n] .* w_vectors[k]
            end
        end
    end
    return search_vectors_list        
end

function convergence_test(Ritz_values, Ritz_vectors::Vector{Vector{ComplexF64}}, residuals, lambda_min, lambda_max, epsilon_convergence)
    converged_target_vectors = Vector{ComplexF64}[]
    converged_target_values = Float64[]
    not_converged_residuals = Float64[]
    @inbounds for i in eachindex(Ritz_values) 
        if lambda_min < Ritz_values[i] < lambda_max
            if residuals[i] < epsilon_convergence
                push!(converged_target_vectors, Ritz_vectors[i])
                push!(converged_target_values, Ritz_values[i])
            else
                push!(not_converged_residuals, residuals[i])
            end
        end
    end
    return converged_target_values, converged_target_vectors, not_converged_residuals     
end

function Rayleigh_Ritz_pairs_residuals!(Ritz_matrix::Matrix{ComplexF64}, Ritz_vectors::Vector{Vector{ComplexF64}}, search_vectors_list::Matrix{ComplexF64}, provisional_vector::Vector{Vector{ComplexF64}}, hamiltonian_matrix::SparseMatrixCSC{ComplexF64, Int64})
    Rayleigh_Ritz_matrix_building!(Ritz_matrix, search_vectors_list, provisional_vector, hamiltonian_matrix) 
    Rayleigh_Ritz_pairs = eigen(Ritz_matrix)
    Ritz_values = real(Rayleigh_Ritz_pairs.values)
    residuals = zeros(size(search_vectors_list, 2))
    @inbounds Threads.@threads for i in 1:size(search_vectors_list, 2)
        mul!(Ritz_vectors[i], search_vectors_list, Rayleigh_Ritz_pairs.vectors[:, i])
        residuals[i] = norm(hamiltonian_matrix * Ritz_vectors[i] .- Ritz_values[i] .* Ritz_vectors[i])
    end
    return Ritz_values, residuals    
end

function Rayleigh_Ritz_matrix_building!(Ritz_matrix::Matrix{ComplexF64}, search_vectors_list::Matrix{ComplexF64}, provisional_vector::Vector{Vector{ComplexF64}}, hamiltonian_matrix::SparseMatrixCSC{ComplexF64, Int64})
    @inbounds Threads.@threads for i in axes(search_vectors_list, 2)
        mul!(provisional_vector[Threads.threadid()], hamiltonian_matrix, search_vectors_list[:, i]) 
        @inbounds for j in i:size(search_vectors_list, 2)
            Ritz_matrix[j, i] = dot(search_vectors_list[:, j], provisional_vector[Threads.threadid()])
            if j != i
                Ritz_matrix[i, j] = conj(Ritz_matrix[j, i]) 
            end
        end
    end   
end
