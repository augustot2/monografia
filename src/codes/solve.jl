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
