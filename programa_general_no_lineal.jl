using LinearAlgebra
function Resolucion_PVF_no_lineal(ϕ,a,b,α,β,N)
    include("fdjac.jl")
    include("levenberg.jl")
    include("matriz_diff_chebyshev.jl")

    x0 = collect(range(α,β,length=N+1)) # primera raíz
    n = length(x0)-1
    x,D = matriz_diff_cheb(n,a,b)

    function sistema_f(u)
        s = (D^2)*u - ϕ.(x,u,D*u)
        s[1] = u[1]-α; s[n+1] = u[n+1]-β
        return s
    end

    u = levenberg(sistema_f,x0)
    return x,u[end]
end