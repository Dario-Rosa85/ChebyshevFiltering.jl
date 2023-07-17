###### new procedure with @spawn

function polynomial_filtering(search_vector_numbers, polynomial_degree_optim, full_coeff, hamiltonian_matrix::SparseMatrixCSC{Float64, Int64}, lambda_max, lambda_min, epsilon_convergence, log_path, log_file_name)
    ########Random vectors inizialization
    search_vectors_list = orthogonalize_QR(randn(Float64, size(hamiltonian_matrix, 2), search_vector_numbers))
    Ritz_matrix = Matrix{Float64}(undef, size(search_vectors_list, 2), size(search_vectors_list, 2))
    # provisional_vector = map(x -> zeros(Float64, size(hamiltonian_matrix, 2)), 1:Threads.nthreads())
    u_vectors = map(x -> zeros(Float64, size(hamiltonian_matrix, 2)), 1:search_vector_numbers)
    w_vectors = map(x -> zeros(Float64, size(hamiltonian_matrix, 2)), 1:search_vector_numbers)
    Ritz_vectors = map(x -> zeros(Float64, size(hamiltonian_matrix, 2)), 1:search_vector_numbers) 
    
    convergence_reached = "false"
    number_of_iterations = 0
    while convergence_reached == "false"
        ##########################################################################################################
        ########Filtering step and orthogonalization
        filtering_step!(search_vectors_list, u_vectors, w_vectors, polynomial_degree_optim, full_coeff, hamiltonian_matrix)
        search_vectors_list = orthogonalize_QR(search_vectors_list)
        ##########################################################################################################
        ##########################################################################################################
        ########Ritz pairs and residuals computation
        Ritz_values, residuals = Rayleigh_Ritz_pairs_residuals!(Ritz_matrix, Ritz_vectors, search_vectors_list, hamiltonian_matrix)
        ##########################################################################################################
        ##########################################################################################################
        ########Convergence test
        converged_target_values, converged_target_vectors, not_converged_residuals = convergence_test(Ritz_values, Ritz_vectors, residuals, lambda_min, lambda_max, epsilon_convergence)
        if length(not_converged_residuals) == 0
            convergence_reached = "true"
            if log_path != "none" && log_file_name != "none"
                open(joinpath(log_path, log_file_name), "a") do io
                    println(io, "Convergence reached. Time is ", Dates.now())
                end
            end
            return converged_target_values, converged_target_vectors
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

function filtering_step!(search_vectors_list::Matrix{Float64}, u_vectors::Vector{Vector{Float64}}, w_vectors::Vector{Vector{Float64}}, polynomial_degree_optim, full_coeff, hamiltonian_matrix::SparseMatrixCSC{Float64, Int64})
    chunks = Iterators.partition(axes(search_vectors_list, 2), Int64(floor(length(axes(search_vectors_list, 2)) / (2 * Threads.nthreads()))))
    Threads.@sync for chunk in chunks
        Threads.@spawn @inbounds for k in chunk
            mul!(u_vectors[k], hamiltonian_matrix, search_vectors_list[:, k])
            provisional_vector = search_vectors_list[:, k]
            mul!(provisional_vector, hamiltonian_matrix, u_vectors[k], 2.0, -1.0)
            w_vectors[k] = provisional_vector
            search_vectors_list[:, k] .= full_coeff[1] * search_vectors_list[:, k] .+ full_coeff[2] .* u_vectors[k] .+ full_coeff[3] .* w_vectors[k]
            Threads.@spawn @inbounds for n in 4:polynomial_degree_optim 
                provisional_vector = u_vectors[k]
                mul!(provisional_vector, hamiltonian_matrix, w_vectors[k], 2.0, -1.0)
                u_vectors[k] = w_vectors[k]
                w_vectors[k] = provisional_vector
                search_vectors_list[:, k] .+= full_coeff[n] .* w_vectors[k]
            end                 
        end
    end 
end


function Rayleigh_Ritz_pairs_residuals!(Ritz_matrix::Matrix{Float64}, Ritz_vectors::Vector{Vector{Float64}}, search_vectors_list::Matrix{Float64}, hamiltonian_matrix::SparseMatrixCSC{Float64, Int64})
    Rayleigh_Ritz_matrix_building!(Ritz_matrix, search_vectors_list, hamiltonian_matrix) 
    Rayleigh_Ritz_pairs = eigen(Ritz_matrix)
    Ritz_values = real(Rayleigh_Ritz_pairs.values)
    residuals = zeros(size(search_vectors_list, 2))
    @inbounds Threads.@threads for i in 1:size(search_vectors_list, 2)
        mul!(Ritz_vectors[i], search_vectors_list, Rayleigh_Ritz_pairs.vectors[:, i])
        residuals[i] = norm(hamiltonian_matrix * Ritz_vectors[i] .- Ritz_values[i] .* Ritz_vectors[i])
    end
    return Ritz_values, residuals    
end

function Rayleigh_Ritz_matrix_building_single_chunk!(Ritz_matrix::Matrix{Float64}, search_vectors_list::Matrix{Float64}, hamiltonian_matrix::SparseMatrixCSC{Float64, Int64}, chunk)
    @inbounds for i in chunk
        temp_vector = hamiltonian_matrix * search_vectors_list[:, i]
        @inbounds for j in i:size(search_vectors_list, 2)
            Ritz_matrix[j, i] = dot(search_vectors_list[:, j], temp_vector)
            if j != i
                Ritz_matrix[i, j] = conj(Ritz_matrix[j, i]) 
            end
        end
    end 
end

function Rayleigh_Ritz_matrix_building!(Ritz_matrix::Matrix{Float64}, search_vectors_list::Matrix{Float64}, hamiltonian_matrix::SparseMatrixCSC{Float64, Int64})
    chunks = Iterators.partition(axes(search_vectors_list,2), Int64(floor(length(axes(search_vectors_list,2)) / (2 * Threads.nthreads()))))
    Threads.@sync for chunk in chunks
        Threads.@spawn Rayleigh_Ritz_matrix_building_single_chunk!(Ritz_matrix, search_vectors_list, hamiltonian_matrix, chunk)
    end
end

function convergence_test(Ritz_values, Ritz_vectors::Vector{Vector{Float64}}, residuals, lambda_min, lambda_max, epsilon_convergence)
    converged_target_vectors = Vector{Float64}[]
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

