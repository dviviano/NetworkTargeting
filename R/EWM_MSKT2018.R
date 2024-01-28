Est_max_score <- function(
    X,B=1,g_i,params = NA, model_only=FALSE, tolerance_constraint = 10**(-8), scale_Y = F,
    max_treated_units = dim(X)[1], maxtime = 3000,
    min_treated_units = 0, bound_param = NA, numcores = 12, verbose = F)
  {

  n <- dim(X)[1]
  p <- dim(X)[2]
  XX <- cbind(1, X)
  maximum_X <- max(apply(XX, 1, function(h) max(abs(h))))
  XX <- XX/maximum_X

  ## Collapse elements with same covariates
  unique_values <- unique(XX) ## New covariates

  ## Consider the indexes of the unique values
  index_unique <- apply(XX,1, function(y) which(apply(unique_values, 1, function(x) all(y == x))))

  ## Collapse gi accordingly
  g_i <- sapply(c(1:dim(unique_values)[1]), function(x) sum(g_i[which(x == index_unique)]))
  XX <- unique_values
  n <- length(g_i)
  n_indexes <- sapply(c(1:n), function(x) sum(index_unique == x))

  #C <- B*apply(XX, 1, function(x) sum(abs(x)))
  C <- B*max(apply(XX, 1, function(x) sum(abs(x))))
  XX <- XX/C

  #XX <- t(sapply(c(1:n),  function (x)  XX[x,]/C[x]))
  model <- list()

  ## Specify the constraints
  ## Consider the vector of (z_1, ..., z_n,beta_0 , beta_1, ..., beta_p)
  model$A<- rbind(cbind(diag(1, nrow = n), -XX), cbind(diag(1, nrow = n), -XX), c(n_indexes, rep(0, p +1)))
  model$obj<- c(g_i, rep(0, p + 1))
  model$modelsense<-'max'

  ## tolerance_constraint enters here
  model$rhs<- c(rep(1 - tolerance_constraint, dim(model$A)[1]/2), rep(tolerance_constraint, dim(model$A)[1]/2), max_treated_units)
  model$sense<- c(rep('<=', dim(model$A)[1]/2), rep('>', dim(model$A)[1]/2), '<=')

  if(min_treated_units > 0){
    model$A<- rbind(model$A, c(n_indexes, rep(0, p +1)))
    model$rhs<- c(model$rhs, min_treated_units)
    model$sense <- c(model$sense, '>=')
  }

  model$vtype<- c(rep('B', n), rep('C', p+1))

  # Put bounds on the parameter space (If commented, parameter space = real line)
  model$ub<- c(rep(1,n), rep(B,1+p))
  model$lb<- c(rep(0,n), rep(-B,p + 1))
  if(is.na(bound_param) == F){
    model$ub<- c(rep(1,n), bound_param[,1])
    model$lb<- c(rep(0,n), bound_param[,2])
  }

  ## Run gurobi only on one thread by default
  if(is.na(params)){
    params<- list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit =maxtime, BarConvTol = exp(-2),
                  Threads = numcores, Disconnected=0,Heuristics=0, NodeFileStart  = 0.5, OutputFlag = as.numeric(verbose))}

  #params<- list(IntFeasTol = 1e-6, FeasibilityTol = 1e-6, TimeLimit =TimeLimit, BarConvTol = exp(-2))
  if (model_only) return(model)
  result<- gurobi::gurobi(model, params = params)

  ## Compute welfare
  beta_hat <- result$x[(n+1):(n+p+1)]

  ## Since X was rescaled by the same constant we can use X (not XX) here
  pi_est <- apply(cbind(1,X), 1, function(x) ifelse(x%*%beta_hat > 0, 1, 0))

  ## Checking instability of the solution
  agreement <- 1 - mean(apply(XX, 1, function(x) ifelse(x%*%beta_hat > 0, 1, 0)) != result$x[1:dim(XX)[1]])
  return(list(pi_est, beta_hat,  obj_est = result$objval, agreement = agreement, result))
}
