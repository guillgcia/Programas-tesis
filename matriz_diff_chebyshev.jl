using LinearAlgebra
function matriz_diff_cheb(N,a,b) # Genera una matriz
                                 # de diferenciación de orden (N+1)*(N+1)
    x = ones(N+1)
    for k in 0:N
        x[k+1] = -cos((pi*k)./N) # Nodos de Chebyshev en [-1,1]   
    end
    X = @. (b+a)/2 + (b-a).*(x)./2 # Nodos de Chebyshev en [a,b]
    B = zeros(N+1,N+1)
    B[1,1] = -(2 .*(N.^2)+1)./6
    B[1,N+1] = -((-1)^(N))/2
    B[N+1,1] = -B[1,N+1]
    B[N+1,N+1] = -(B[1,1])
    for j in 1:(N-1)
        B[j+1,j+1] = -(x[j+1])./(2 .*(1 .- (x[j+1]).^2))
        for i in 1:(N-1)
            if (i != j)
                B[j+1,i+1] = -((-1).^((i)+(j)))/(x[i+1]-x[j+1])
            end
            B[1,i+1] = -(2(-1).^(i))/(1+x[i+1])
            B[N+1,i+1] = (2(-1).^(i+N))/(1-x[i+1])
            B[j+1,1] = ((-1).^(j))/(2(1+x[j+1]))
            B[j+1,N+1] = -((-1).^(j+N))/(2(1-x[j+1]))
        end
    end
    B = 2*B/(b-a)   
    return X,B # Da como salida la matriz de diferenciación
               # y los nodos de Chebyshev en un intervalo [a,b]
end
