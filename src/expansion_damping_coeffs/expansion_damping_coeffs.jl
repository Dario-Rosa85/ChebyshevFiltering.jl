function Lanczos_coeff(polynomial_order, μ_factor)
    prepend!(map(n -> ((polynomial_order + 1) * sin(π * n / (polynomial_order + 1)) / (π * n))^μ_factor , range(1, polynomial_order)), 1)
end


function Fejer_coeff(polynomial_order)
    map(n -> (polynomial_order - n + 1) / (polynomial_order + 1) , range(0, polynomial_order))
end


function Jackson_coeff(polynomial_order)
    map(n -> polynomial_order^(-1) * ((polynomial_order - n) * cos(π * n / polynomial_order) + sin(π * n / polynomial_order) * cot(π / polynomial_order) ) , range(0, polynomial_order))
end


function expansion_coeff(polynomial_order, e_max, e_min, lambda_max, lambda_min)
    rescaling_factor = 2 / (e_max - e_min)
    shifting_factor = (e_max + e_min) / (e_max - e_min)
    expansion_coeffs = [(1 / π) * (acos(rescaling_factor * lambda_max + shifting_factor) - acos(rescaling_factor * lambda_min + shifting_factor))]
    for n in 1:polynomial_order 
        push!(expansion_coeffs, (2 / (π * n)) * (sin(n * acos(rescaling_factor * lambda_max + shifting_factor)) - sin(n * acos(rescaling_factor * lambda_min + shifting_factor))))
    end
    return expansion_coeffs
end