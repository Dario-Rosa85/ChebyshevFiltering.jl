function hamiltonian_moments(max_degree, stochastic_dimension, hamiltonian_matrix::SparseMatrixCSC{ComplexF64, Int64})
    hilbert_space_dimension = size(hamiltonian_matrix, 2)
    random_states_matrix = orthogonalize_QR(randn(ComplexF64, hilbert_space_dimension, stochastic_dimension))
    random_states = [random_states_matrix[:, i] for i in axes(random_states_matrix, 2)]
    T_m = deepcopy(random_states)
    T_m_plus = [zeros(ComplexF64, hilbert_space_dimension) for i in 1:length(random_states)]
    moments_computed = zeros(Float64, max_degree + 1)
    moments_computed[1] = 1.0
    auxiliary_vector = map(x -> zeros(ComplexF64, hilbert_space_dimension), 1:Threads.nthreads())
    @inbounds Threads.@threads for r in 1:stochastic_dimension
        mul!(T_m_plus[r],  hamiltonian_matrix, T_m[r])
    end
    @inbounds for n in 1:max_degree
        auxiliary_moments_computed = zeros(Float64, Threads.nthreads())
        @inbounds Threads.@threads for r in 1:stochastic_dimension 
            mul!(T_m[r], hamiltonian_matrix, T_m_plus[r], 2, -1)
            auxiliary_moments_computed[Threads.threadid()] += real(dot(random_states[r], T_m_plus[r]))
            auxiliary_vector[Threads.threadid()] = T_m[r]
            T_m[r] = T_m_plus[r]
            T_m_plus[r] = auxiliary_vector[Threads.threadid()]
        end
        moments_computed[n + 1] = sum(auxiliary_moments_computed) / stochastic_dimension
    end
    return moments_computed
end

function KPM_density_moments_given(moments_computed)
    max_degree = length(moments_computed) - 1
    damping_factors = Jackson_coeff(max_degree)
    prefactors = [moments_computed[1] * damping_factors[1]]
    @inbounds for i in 2:length(moments_computed) 
        push!(prefactors, 2 * moments_computed[i] * damping_factors[i])
    end
    density_KPM(x) = ChebyshevT(prefactors)(x) / (π * sqrt(1 - (x)^2))
    return density_KPM
end


function KPM_density(max_degree, stochastic_dimension, hamiltonian_matrix, normalized::String = "no")
    if normalized == "no"
        normalization_factor = size(hamiltonian_matrix, 2)
    else
        normalization_factor = 1 
    end
    moments_computed = normalization_factor .* hamiltonian_moments(max_degree, stochastic_dimension, hamiltonian_matrix)
    return KPM_density_moments_given(moments_computed)
end