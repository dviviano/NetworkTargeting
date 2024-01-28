
## Linear formulation of MILP for Network Interference.

max_score_NEWM_helper <- function(constraint_matrix, XX, impose_max_treated = NA,b, model,
                                  y_s_position, no_neighb, any_neighb, final_gs_noneighb, my_h, c,
                                  tolerance_constraint1, tolerance_constraint2, B,low_b = -B, params = NA, constant, cost = 0,
                                  constraint_matrix_no_n, maximum_time, additional_param_constraints,
                                  model_return = F, min_treated = 0, numcores = 1, focal_units = rep(1, dim(XX)[1]), verbose = F) {

  n <- dim(XX)[1]
  p <- dim(XX)[2] - 1

  ## Impose constraint on the number of treated units
  if (is.na(impose_max_treated) == F) {
    max_treated <- rep(0, ncol(constraint_matrix))
    max_treated[y_s_position[focal_units[-no_neighb] == 1]] <- 1
    b <- c(b, impose_max_treated)
    ## Number of treated units equal to impose_max_treated
    model$sense<- c(model$sense, '<=')
    constraint_matrix <- rbind(constraint_matrix, max_treated)
  }

  if (min_treated > 0) {
    min_treated_vec <- rep(0, ncol(constraint_matrix))
    min_treated_vec[y_s_position[focal_units[-no_neighb] == 1]] <- 1
    b <- c(b, min_treated)
    model$sense<- c(model$sense, '>=')
    constraint_matrix <- rbind(constraint_matrix, min_treated_vec)
  }

  ## Add back individuals with no_neighbors
  if (length(no_neighb) > 1 & length(any_neighb) > 0) {
    ## Add constraint on maximum number of treated units
    if(is.na(impose_max_treated) == F) constraint_matrix[dim(constraint_matrix)[1] - ifelse(min_treated == 0,0,1), 1:(length(no_neighb))] <- as.numeric(focal_units[no_neighb])
    if(min_treated > 0) constraint_matrix[dim(constraint_matrix)[1], 1:(length(no_neighb))] <- as.numeric(focal_units[no_neighb])
    ## Construct cost function and b
    ## note that final_gs_noneighb already incorporates the costs
    c <- c( final_gs_noneighb, c)
    constraint_matrix <- cbind(matrix(0, nrow = nrow(constraint_matrix), ncol = length(no_neighb)), constraint_matrix)
    constraint_matrix <- rbind(constraint_matrix_no_n, constraint_matrix)
    b <- c(rep(tolerance_constraint1, length(no_neighb)), rep(1, length(no_neighb)), b)
    model$sense<- c(rep('>', length(no_neighb)), rep('<=', length(no_neighb)), model$sense)
  }

  if (verbose) print(dim(constraint_matrix))
  if (verbose) print(length(model$sense))

  model$A <- constraint_matrix
  model$obj <- c
  model$modelsense<-'max'
  model$rhs<- b
  model$vtype<- c(rep('B', length(c) - p - 1), rep('C', p+1))

  # Put bounds on the parameter space (If commented, parameter space = real line)
  if(is.na(additional_param_constraints)){
    model$ub<- c(rep(1,length(c) - p - 1), rep(B,1+p))
    model$lb<- c(rep(0,length(c) - p - 1), rep(low_b,1+p)) }  else {
      model$ub<- c(rep(1,length(c) - p - 1), additional_param_constraints[,1])
      model$lb<- c(rep(0,length(c) - p - 1), additional_param_constraints[,2])
    }

  model$objcon <- constant

  ## Run gurobi only on one thread by default
  if(is.na(params)){
    ## maximum time here
    params <- list(
      OutputFlag = as.numeric(verbose),
      IntFeasTol = 1e-9,
      FeasibilityTol = 1e-9,
      TimeLimit =maximum_time,
      Threads = numcores,
      Heuristics=0,
      Cuts=0)
  }

  if(model_return) return(model)
  result<- gurobi::gurobi(model, params = params)
  return(result)

}

