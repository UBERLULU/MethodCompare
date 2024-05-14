var_blup <- function(data, model) {
  # Extract variances parameters from the model
  sig2_u <- as.data.frame(lme4::VarCorr(model))$vcov[1]
  sig2_e <- as.data.frame(lme4::VarCorr(model))$vcov[2]
  
  # Extract model matrices
  Z <- lme4::getME(model, "Z")
  Zt <- lme4::getME(model, "Zt")
  X <- lme4::getME(model, "X")
  D <- sig2_u * diag(length(unique(data$id)))
  R <- sig2_e * diag(nrow(Z))
  invR <- chol2inv(chol(R))
  V <- Z %*% D %*% Zt + R
  
  invV <- chol2inv(chol(V))
  res <- data$y2 - lme4::fixef(model) * X
  
  bi <- D %*% Zt %*% invV %*% res
  
  sumXWX <- 0
  
  for (i in unique(data$id)) {
    Zi <- Z[data$id == i, ]
    
    Vi <- Zi %*% D %*% t(Zi) + R[data$id == i, data$id == i]
    
    invVi <- chol2inv(chol(Vi))
    
    Xi <- as.matrix(X[data$id == i])
    
    sumXWX <- sumXWX + t(Xi) %*% invVi %*% Xi
  }
  
  inv_sumXWX <- chol2inv(chol(sumXWX))
  
  vbi <- D %*% Zt %*% (invV - invV %*% X %*% chol2inv(chol(sumXWX)) %*%
                         t(X) %*% invV) %*% Z %*% D
  
  return(data.frame(id = unique(data$id), v_blup = diag(D - vbi)))
}

