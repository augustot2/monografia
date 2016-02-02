using Jacobi
using PyPlot

#resolve A*X = B
function solver_beta(A,B,nb)
  Abb = A[1:nb,1:nb]
  Abi = A[1:nb,(nb+1):end]
  Aii = A[(nb+1):end,(nb+1):end]
  Aib = A[(nb+1):end,1:nb]
  #########################
  #                       #
  #  |Abb Abi| |xb| = |Bb|#
  #  |Aib Aii| |xi|   |Bi|#
  #########################

  Bi = B[(nb+1):end]
  Bb = B[1:nb]

  Xb = inv(Abb - Abi*inv(Aii)*Aib)* (Bb - Abi*inv(Aii)*Bi)
  Xi = inv(Aii)*Bi - inv(Aii)*transpose(Abi)*Xb

  X = [Xb;Xi]
end

#Jacobi
function ψj(p,E,Q)
    if(p == 1)
     return  (1-E)/2
        elseif(p == 2)
     return  (1+E)/2
    else
     return  (1-E)*(1+E)/4 * jacobi(E, p-3, 1, 1)
    end
end 
#Legendre
function ψle(p,E,Q)
   legendre(E, p-1)
end 
#Lagrange
function ψla(p,E,Q)
    z = zgj(Q)
    lagrange(p,E,z)
end
#Matriz de massa
function Mass_matrix(ψ,Q,M)
    ϕ = zeros(Q,M)
    ξ = zgj(Q)
    w = wgj(ξ,0.,0.)
    for i in 1:M
        for j in 1:Q
            ϕ[j,i] = ψ(i,ξ[j],Q)
        end
    end    
    L = zeros(M,M)
    for i in 1:M
        for j in 1:M
           m= 0.0
            for q in 1:Q
                m = m + ϕ[q,i]*ϕ[q,j]*w[q]
            end
            L[i,j] = m
        end
    end
    return L
end 


#Plota matriz de massa

function plot_matrix(M, posneg=true,eps=1e-5)

    nrow = size(M,1)
    ncol = size(M,2)

    ntot = nrow * ncol

    A = reshape(M, ntot)

    irow = zeros(Int,ntot)
    icol = zeros(Int,ntot)

    er = maxabs(A)*eps
    
    cnt = 1
    for j = 1:ncol
        for i = 1:nrow
            irow[cnt] = i
            icol[cnt] = j
            cnt = cnt + 1
        end
    end

    if posneg
        ipos = A .> er
        ineg = A .< -er
        
        plot(icol[ipos], nrow + 1 - irow[ipos], "ro",label="positive")
        plot(icol[ineg], nrow + 1 - irow[ineg], "bo",label="negative")
        legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)        
    else
        idx = abs(A) .> er
        plot(icol[idx], nrow + 1 - irow[idx], "o", color="black")
    end
    return
end



