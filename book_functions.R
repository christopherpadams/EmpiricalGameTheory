
# Chapter 4

bs = function(init, fn, y, X, R=100, trace=FALSE,
              optim_trace = 0, write=FALSE, seed=TRUE) {
  # Bootstrap standard errors
  if(seed) {set.seed(123456789)}
  J = length(init)
  X = as.matrix(X)
  N = dim(y)[1]
  if(is.null(N)) {N = length(y)}
  res_mat = matrix(NA, R, J)
  for(r in 1:R) {
    index_r = sample(1:N, N, replace = TRUE)
    if(is.null(dim(y))) {
      a_r = optim(par = init, 
                  fn = fn, 
                  y = y[index_r], 
                  X = X[index_r,],
                  control = list(trace = optim_trace, 
                                 maxit = 1000000))
    } 
    else {
      a_r = optim(par=init, fn=fn, y=y[index_r,], 
                  X=X[index_r,],
                  control = list(trace = optim_trace, maxit = 1000000))
    }
    res_mat[r,] = a_r$par
    if(trace) {
      print(res_mat[r,])
      print(r/R)
    }
    if(write) {
      write.csv(res_mat, "res_mat_temp.csv")
    }
  }
  return(res_mat)
} 


# Chapter 5


f_gmm = function(G, K) {
  G = as.matrix(G)
  N = dim(G)[2]
  if (K == dim(G)[1]) {
    # a check that the matrix G has K rows
    g = rowMeans(G, na.rm = TRUE)
    W = try(solve(G%*%t(G)/N), silent = TRUE)
    # try() lets the function work even if there is an error
    if (is.matrix(W)) {
      # if there is no error, W is a matrix.
      return(t(g)%*%W%*%g)
    }
    else {
      # allow estimation assuming W is identity
      return(t(g)%*%g)
    }
  }
  else {
    return("ERROR: incorrect dimension")
  }
}

bs_gmm = function(init, 
                  fn, 
                  X, 
                  R = 100,  
                  trace = FALSE,
                  trace_optim = 0,
                  write = FALSE) {
  # Bootstrap standard errors
  J = length(init)
  X = as.data.frame(X)
  N = dim(X)[1]
  res_mat = matrix(NA, R, J)
  for(r in 1:R) {
    index_r = sample(1:N, 
                     N, 
                     replace = TRUE)
    a_r = optim(par=init, 
                fn=fn, 
                X=X[index_r,],
                control = list(trace = trace_optim,
                               maxit = 1000000))  
    res_mat[r,] = a_r$par
    if(trace) {
      print(res_mat[r,])
      print(r/R)
    }
    if(write) {
      write.csv(res_mat, "res_mat_temp.csv")
    }
  }
  return(res_mat)
} 

f_loglik_mix = function(X, y, beta_1, beta_2, alpha_1, alpha_2, rho) {
  epsilon = 1e-10
  Lik = f_entry_mix(X, beta_1, beta_2, alpha_1, alpha_2, rho)
  return((y[,1]==0 & y[,2]==0)*log(Lik$p_00 + epsilon) + 
           (y[,1] == 0 & y[,2] == 1)*log(Lik$p_01 + epsilon) + 
           (y[,1]==1 & y[,2]==1)*log(Lik$p_11 + epsilon) + 
           (y[,1]==1 & y[,2]==0)*log(1 - Lik$p_00 - Lik$p_01 - Lik$p_11 + epsilon))
}

f_loglik_mix_int = function(par, X, y) {
  X = as.matrix(cbind(1, X))
  J = dim(X)[2]
  beta_1 = par[1:J]
  beta_2 = par[(J+1):(2*J)]
  alpha_1 = exp(par[2*J+1])
  alpha_2 = exp(par[2*J+2])
  rho = -1 + 2*exp(par[2*J+3])/(1 + exp(par[2*J+3]))
  return(-sum(f_loglik_mix(X, y, beta_1, beta_2, alpha_1, alpha_2, rho)))
}

# Chapter 6

f_loglik_spne = function(X, y, beta_1, beta_2, alpha_1, alpha_2, rho) {
  epsilon = 1e-10
  Lik = f_entry_spne(X, beta_1, beta_2, alpha_1, alpha_2, rho)
  return((y[,1]==0 & y[,2]==0)*log(Lik$p_00 + epsilon) + 
           (y[,1]==0 & y[,2]==1)*log(Lik$p_01 + epsilon) + 
           (y[,1]==1 & y[,2]==1)*log(Lik$p_11 + epsilon) + 
           (y[,1]==1 & y[,2]==0)*log(1 - 
                                       Lik$p_00 - 
                                       Lik$p_01 - 
                                       Lik$p_11 + epsilon))
}

f_loglik_spne_int = function(par, X, y) {
  X = as.matrix(cbind(1, X))
  J = dim(X)[2]
  beta_1 = par[1:J]
  beta_2 = par[(J+1):(2*J)]
  alpha_1 = exp(par[2*J+1])
  alpha_2 = exp(par[2*J+2])
  rho = -1 + 2*exp(par[2*J+3])/(1 + exp(par[2*J+3]))
  return(-mean(f_loglik_spne(X, y, beta_1, beta_2, alpha_1, alpha_2, rho)))
}



# Chapter 9

bs2 = function(init1, init2, fn1, fn2, y, X, R=100, trace=FALSE,
               optim_trace = 0, write=FALSE, seed=TRUE) {
  # Bootstrap standard errors
  if(seed) {set.seed(123456789)}
  J1 = length(init1)
  J2 = length(init2)
  X = as.matrix(X)
  J3 = dim(X)[2] + 1
  N = dim(y)[1]
  if(is.null(N)) {N = length(y)}
  res_mat = matrix(NA, R, J2)
  for(r in 1:R) {
    index_r = sample(1:N, N, replace = TRUE)
    if(is.null(dim(y))) {
      a_r1 = optim(par=init1, fn=fn1, y=y[index_r], 
                   X=X[index_r,],
                   control = list(trace = 0, maxit = 1000000))
      a_r2 = optim(par=init2, fn=fn2, y=y[index_r], 
                   X=X[index_r,],
                   beta_1t = a_r1$par[1:J3],
                   beta_2t = a_r1$par[(J3+1):(2*J3)],
                   rhot = -1 + 
                     2*exp(a_r1$par[2*J3+1])/(1 + exp(a_r1$par[2*J3+1])),
                   control = list(trace = optim_trace, maxit = 1000000))
    } 
    else {
      a_r1 = optim(par=init1, fn=fn1, y=y[index_r,], 
                   X=X[index_r,],
                   control = list(trace = 0, maxit = 1000000))
      a_r2 = optim(par=init2, fn=fn2, y=y[index_r,], 
                   X=X[index_r,],
                   beta_1t = a_r1$par[1:J3],
                   beta_2t = a_r1$par[(J3+1):(2*J3)],
                   rhot = -1 + 
                     2*exp(a_r1$par[2*J3+1])/(1 + exp(a_r1$par[2*J3+1])),
                   control = list(trace = optim_trace, maxit = 1000000))
    }
    res_mat[r,] = a_r2$par
    if(trace) {
      print(res_mat[r,])
      print(r/R)
    }
    if(write) {
      write.csv(res_mat, "res_mat_temp.csv")
    }
  }
  return(res_mat)
} 
