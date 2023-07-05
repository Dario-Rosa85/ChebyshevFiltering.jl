function renormalization_hamiltonian!(hamiltonian_to_renormalize::SparseMatrixCSC{ComplexF64, Int64})
    matrix_size = size(hamiltonian_to_renormalize, 2)
    sparse_identity = spzeros(ComplexF64, matrix_size, matrix_size)
    @inbounds for i in 1:matrix_size
        sparse_identity[i, i] = 1.0 + im * 0.0
    end
    e_max = real(first(first(eigsolve(hamiltonian_to_renormalize, rand(ComplexF64, matrix_size), 1, :LR))))
    e_min = real(first(first(eigsolve(hamiltonian_to_renormalize, rand(ComplexF64, matrix_size), 1, :SR))))
    rescaling_factor = 2 / (e_max - e_min)
    shifting_factor = (e_max + e_min) / (e_max - e_min)
    hamiltonian_to_renormalize .= rescaling_factor .* hamiltonian_to_renormalize .- shifting_factor .* sparse_identity
end

function renormalization_hamiltonian!(hamiltonian_to_renormalize::SparseMatrixCSC{Float64, Int64})
    matrix_size = size(hamiltonian_to_renormalize, 2)
    sparse_identity = spzeros(Float64, matrix_size, matrix_size)
    @inbounds for i in 1:matrix_size
        sparse_identity[i, i] = 1.0 + im * 0.0
    end
    e_max = real(first(first(eigsolve(hamiltonian_to_renormalize, rand(Float64, matrix_size), 1, :LR))))
    e_min = real(first(first(eigsolve(hamiltonian_to_renormalize, rand(Float64, matrix_size), 1, :SR))))
    rescaling_factor = 2 / (e_max - e_min)
    shifting_factor = (e_max + e_min) / (e_max - e_min)
    hamiltonian_to_renormalize .= rescaling_factor .* hamiltonian_to_renormalize .- shifting_factor .* sparse_identity
end