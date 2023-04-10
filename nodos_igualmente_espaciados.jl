function nodos_espaciados(N,a,b)
    x = zeros(N+1)
    for i in 0:N
        x[i+1] = a*(N-i)/N + (b*i)/N
    end
    return x
end