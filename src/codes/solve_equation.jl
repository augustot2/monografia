function ψj(p,E,Q)
    if(p == 1)
     return  (1-E)/2
        elseif(p == Q)
     return  (1+E)/2
    else
     return  (1-E)*(1+E)/4 * jacobi(E, p-2, 1, 1)
    end
end
function ψle(p,E,Q)
   legendre(E, p-1)
end
function ψla(p,E,Q)
    z = zgj(Q)
    lagrange(p,E,z)
end

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
