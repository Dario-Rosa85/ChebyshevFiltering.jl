function orthogonalize_SVD(matrix_of_vectors)
    matrix_rank = rank(matrix_of_vectors)
    U_factor = svd(matrix_of_vectors).U
    U_factor = U_factor[:,1:matrix_rank]
    return U_factor
end

function orthogonalize_QR(matrix_of_vectors)
    qr(matrix_of_vectors).Q * Matrix(I, size(matrix_of_vectors)...)
end