using FFTW, ToeplitzMatrices, AbstractFFTs
function matriz_diff(M)
    e = zeros(M)
    e[2] = 2/3
    e[3] = -1/12
    e[M-1] = -e[3]
    e[M] = -e[2]
    D = Toeplitz(-e,e) #Genera la matriz de diferenciación
    return D
end
