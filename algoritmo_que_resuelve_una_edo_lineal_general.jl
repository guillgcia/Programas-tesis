using LinearAlgebra
include("matriz_diff_chebyshev.jl")
function Resolucion_PVF(P,Q,R,N,a,b,α,β)
    function vector_diagonal(x)
        n = length(x)
        D = zeros(n,n)
        for i in 1:n
            D[i,i] = x[i]
        end
        return D
    end;
    A = Matrix{Float64}(I, N+1, N+1); V = zeros(N+1)
    V[1] = α; V[N+1] = β
    X,D = matriz_diff_cheb(N,a,b)
    p = vector_diagonal(P.(X)); q = vector_diagonal(Q.(X)); r = R.(X)
    L = D^2 + p*D + q*A
    A[2:N,:] = L[2:N,:]
    V[2:N] = r[2:N]
    Y = A\V
    return X,Y
end;
