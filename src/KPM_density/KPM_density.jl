function hamiltonian_moments_single_chunk!(T_m, T_m_plus, random_states, hamiltonian_matrix, chunk)
    auxiliary_moments_computed = 0.0
    @inbounds for r in chunk 
        mul!(T_m[r], hamiltonian_matrix, T_m_plus[r], 2, -1)
        auxiliary_moments_computed += real(dot(random_states[r], T_m_plus[r]))
        auxiliary_vector = copy(T_m[r])
        T_m[r] = T_m_plus[r]
        T_m_plus[r] = auxiliary_vector
    end
    return auxiliary_moments_computed
end

function hamiltonian_moments(max_degree, stochastic_dimension, hamiltonian_matrix::SparseMatrixCSC{ComplexF64, Int64})
    hilbert_space_dimension = size(hamiltonian_matrix, 2)
    random_states = eachcol(orthogonalize_QR(randn(ComplexF64, hilbert_space_dimension, stochastic_dimension)))
    T_m = deepcopy(random_states)
    T_m_plus = [zeros(ComplexF64, hilbert_space_dimension) for i in eachindex(random_states)]
    moments_computed = zeros(Float64, max_degree + 1)
    moments_computed[1] = 1.0
    @inbounds Threads.@threads for r in 1:stochastic_dimension
        mul!(T_m_plus[r],  hamiltonian_matrix, T_m[r])
    end
    chunks = Iterators.partition(1:stochastic_dimension, Int64(floor(length(1:stochastic_dimension)/Threads.nthreads())))
    @inbounds for n in 1:max_degree
        tasks = [Threads.@spawn hamiltonian_moments_single_chunk!(T_m, T_m_plus, random_states, hamiltonian_matrix, chunk) for chunk in chunks]
        moments_computed[n + 1] = sum(fetch.(tasks)) / stochastic_dimension
    end
    return moments_computed
end

function hamiltonian_moments(max_degree, stochastic_dimension, hamiltonian_matrix::SparseMatrixCSC{Float64, Int64})
    hilbert_space_dimension = size(hamiltonian_matrix, 2)
    random_states = eachcol(orthogonalize_QR(randn(Float64, hilbert_space_dimension, stochastic_dimension)))
    T_m = deepcopy(random_states)
    T_m_plus = [zeros(Float64, hilbert_space_dimension) for i in eachindex(random_states)]
    moments_computed = zeros(Float64, max_degree + 1)
    moments_computed[1] = 1.0
    @inbounds Threads.@threads for r in 1:stochastic_dimension
        mul!(T_m_plus[r],  hamiltonian_matrix, T_m[r])
    end
    chunks = Iterators.partition(1:stochastic_dimension, Int64(floor(length(1:stochastic_dimension)/Threads.nthreads())))
    @inbounds for n in 1:max_degree
        tasks = [Threads.@spawn hamiltonian_moments_single_chunk!(T_m, T_m_plus, random_states, hamiltonian_matrix, chunk) for chunk in chunks]
        moments_computed[n + 1] = sum(fetch.(tasks)) / stochastic_dimension
    end
    return moments_computed
end

function KPM_density_moments_given(moments_computed)
    max_degree = length(moments_computed) - 1
    damping_factors = Jackson_coeff(max_degree)
    prefactors = 2 .* moments_computed .* damping_factors
    prefactors[1] /= 2
    density_KPM(x) = ChebyshevT(prefactors)(x) / (π * sqrt(1 - (x)^2))
    return density_KPM
end

function KPM_density(hamiltonian_matrix, max_degree, stochastic_dimension)
    moments_computed = size(hamiltonian_matrix, 2) .* hamiltonian_moments(max_degree, stochastic_dimension, hamiltonian_matrix)
    return KPM_density_moments_given(moments_computed)
end
