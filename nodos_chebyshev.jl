using LinearAlgebra
function nodos_cheb(N)  # Nodos de chebyshev
    x = ones(N+1)            
    for k in 0:N
        x[k+1] = cos((pi*k)./N)   
    end
    return x
end
