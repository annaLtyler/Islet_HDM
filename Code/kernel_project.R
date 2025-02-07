# This function "removes a direction from "regresses out" a direction from a kernal matrix.
# The vector 'v' should contain scores that you want to remove from span of the columnx of 'K'

kernel_project <- function(K, v){
  
  # Normalize v
  v_hat = v / sqrt(sum(v^2))
  
  # Preject kernel to orthogonal complement of v
  right_mult = (K %*% v) %*% t(v)
  
  K_new = K - v %*% (t(v) %*% K) - right_mult + v %*% (t(v) %*% right_mult)
  
  return(K_new)
}
