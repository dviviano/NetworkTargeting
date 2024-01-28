#' Estimate optimal policy
#'
#' This function computes the optimal policy from a network.
#' The policy is the maximum score function.
#' Estimation uses the Mixed integer program solver Gurobi. Installation is available here: https://www.gurobi.com/ - free access for academic users.
#'
#' @param Y vector of outcomes. Y must contain missing values for those units whose neighbors are unknown but who are neighbors of focal points;
#' @param X variables for targeting;
#' @param Z variables used for estimating the conditional mean function and the propensity score;
#' @param Z_int variables included in the regression that interact with individual and neighbors' covariates;
#' @param D binary treatment assignments;
#' @param W_0 adjacency matrix with 0/1 entries;
#' @param W_0_train adjacency matrix for training. Default is that W_0_train equals W_0.
#' @param cost cost of the treatment, default is zero;
#' @param method either "MILP" , "MILP_graph_clustering", "surrogate", or "MILP_sparse". The first return exact solutions but is computationally demanding; the second return approximate solutions obtained by first partitioning the graph into approximately independent components and then computing the solutions in each sub-graph. It both reports each separate solution and the best solution out of sample. 'surrogate' is strongly recommended for graphs with n >= 500, it solves an approximate surrogate objective function. 'MILP_sparse' solves an exact program, but it uses less RAM (using spars matrices).
#' @param focal_units vector of length n, indicating with 1 whether an individual is focal, and zero otherwise. An individual should be focal if all of her neighbors are observed. By default, all individuals are focal.
#' @param model model used for estimation. Any of 'penalized_regression', 'linear_regression', 'random_forest';
#' @param family   family of the outcome, either gaussian or binary;
#' @param maximum_treated maximum number of treated units. NA indicates that no constraint is imposed;
#' @param maximum_time time limit;
#' @param params list with additional parameters (see references)
#' @return final_coef: coefficients of the policy;
#' @return policies: individual treatment assignments
#' @return solutions_clustering: matrix with solutions when MILP_graph_clustering is selected
#' @return partitions: list of indexes with partitions of the graph if MILP_graph_clustering is selected
#' @export
#'


NetworkTargeting <- function(Y, X, Z, Z_int, D, W_0, W_0_train, cost = 0, method = 'MILP',
                             focal_units = rep(1, length(Y)),
                             model = 'penalized_regression', family = 'gaussian', maximum_treated  = NA,
                             maximum_time, params = list()){

  ## handle data missingness

  stopifnot("All values of D must be non-missing"= sum(is.na(D)) == 0)

  stopifnot("All values of X must be non-missing"= sum(apply(X, 1, function(x) sum(is.na(x)) > 0)) == 0)

  missing_Z = apply(Z, 1, function(x) sum(is.na(x)) > 0)
  if(sum(missing_Z) > 0) {
    focal_units[missing_Z] = 0
    print("Some Z are missing! Units with missing Z are set to non-focal units; missing values are imputed with column means")
    for (i in 1:dim(Z)[2]) {
      Z[is.na(Z[,i]), i] <- mean(Z[,i], na.rm = TRUE)
    }
    stopifnot(sum(apply(Z, 1, function(x) sum(is.na(x)) > 0))==0)
  }

  missing_Z_int = apply(Z_int, 1, function(x) sum(is.na(x)) > 0)
  if(sum(missing_Z_int) > 0) {
    focal_units[missing_Z_int] = 0
    print("Some Z_int are missing! Units with missing Z_int are set to non-focal units; missing values are imputed with column means")
    for (i in 1:dim(Z_int)[2]) {
      Z_int[is.na(Z_int[,i]), i] <- mean(Z_int[,i], na.rm = TRUE)
    }
    stopifnot(sum(apply(Z_int, 1, function(x) sum(is.na(x)) > 0))==0)
  }

  missing_W_0 = apply(W_0, 1, function(x) sum(is.na(x)) > 0)
  if(sum(missing_W_0) > 0) {
    focal_units[missing_W_0] = 0
    print("Some W_0 are missing! Units with missing W_0 are set to non-focal units; missing values are imputed with zero")
    for (i in 1:dim(W_0)[2]) {
      W_0[is.na(W_0[,i]), i] <- 0
    }
    stopifnot(sum(apply(W_0, 1, function(x) sum(is.na(x)) > 0))==0)
  }

  missing_W_0_train = apply(W_0_train, 1, function(x) sum(is.na(x)) > 0)
  if(sum(missing_W_0_train) > 0) {
    focal_units[missing_W_0_train] = 0
    print("Some W_0_train are missing! Units with missing W_0_train are set to non-focal units; missing values are imputed with zero")
    for (i in 1:dim(W_0_train)[2]) {
      W_0_train[is.na(W_0_train[,i]), i] <- 0
    }
    stopifnot(sum(apply(W_0_train, 1, function(x) sum(is.na(x)) > 0))==0)
  }

  missing_Y = sapply(Y, is.na)
  if(sum(missing_Y) > 0) focal_units[missing_Y] = 0

  ## setup

  my_params = create_params(params)
  if(method == 'MILP' & length(Y) > 300) warning('Large sample size and computation time can be long. Select method = surrogate to reduce the computational time.')
  if(method %in% c('MILP', 'MILP_graph_clustering')){
  compute_results = optimal_network_linear_rule(Y,X,Z, Z_int, D,W_0,  impose_max_treated = maximum_treated,
                                                method =  method, m1 = my_params$m1,
                                                model = model, cost = cost, cut_high_neighb = my_params$cut_high_neighb,
                                                tolerance_constraint1 = my_params$tolerance_constraint1,
                                                tolerance_constraint2 = my_params$tolerance_constraint1,
                                                B = my_params$B, low_b = my_params$low_b,
                                                params = my_params$params_gurobi,
                                                monotonicity_c = my_params$monotonicity_c,
                                                maxtime_partitions = my_params$maxtime_partitions,
                                                length_small_samples = my_params$length_small_samples,
                                                doubly_robust = my_params$doubly_robust,
                                                penalized= my_params$penalized,
                                                num_strata = my_params$num_strata,
                                                trimming = my_params$trimming,
                                                family = family,
                                                maximum_time = maximum_time,
                                                additional_param_constraints = my_params$additional_param_constraints,
                                                model_return = my_params$model_return,
                                                W_0_train = W_0_train,
                                                names_rf = my_params$names_rf,
                                                print_checking = my_params$print_checking,
                                                focal_units = focal_units,
                                                return_m1 = my_params$return_m1,
                                                bound_outcomes = my_params$bound_outcomes,
                                                trimmed_subpopulation = my_params$trimmed_subpopulation,
                                                min_treated = my_params$min_treated,
                                                Z_int_neighbors = Z_int,
                                                partitions  = my_params$partitions,
                                                numcores = my_params$numcores,
                                                cross_fitting = my_params$cross_fitting,
                                                slack = my_params$slack,
                                                n_folds = my_params$n_folds,
                                                verbose = my_params$verbose)

  results = list(final_coef = compute_results$final_coef, policies = compute_results$policies,
                 solutions_clustering = compute_results$all_solutions,
                 partitions = compute_results$partitions)

  } else {


    compute_results = newm_sparse_or_surrogate(
    Y, Z, D, W_0,
    W_0_train = W_0_train,
    X = X, Z_int = Z_int,
    cost = cost,
    method = method,
    focal_units = focal_units,
    model = model,
    family = family,
    maximum_treated = maximum_treated,
    maximum_time =  maximum_time,
    params = my_params)

    results = list(final_coef = compute_results$result$beta, policies = compute_results$result$assn,
                   milp_mat = compute_results$result$milp_mat,
                   solutions_clustering = 'Inactive argument. Only active for method = MILP_graph_clustering',
                   partitions = 'Inactive argument. Only active for method = MILP_graph_clustering',
                   all_results = compute_results)

  }

  return(results)
}



#' Estimate welfare effects
#'
#' This function computes the welfare effect of a certain policy assignment.
#'
#' @param policies vector of binary policies of length n
#' @param Y vector of outcomes. Y must contain missing values for those units whose neighbors are unknown but who are neighbors of focal points;
#' @param X variables for targeting;
#' @param Z variables used for estimating the conditional mean function and the propensity score;
#' @param Z_int variables included in the regression that interact with individual and neighbors' covariates;
#' @param D binary treatment assignments;
#' @param W_0 adjacency matrix with 0/1 entries;
#' @param W_0_train adjacency matrix for training. Default is that W_0_train equals W_0.
#' @param cost cost of the treatment, default is zero;
#' @param focal_units vector of length n, indicating with 1 whether an individual is focal, and zero otherwise. An individual should be focal if all of her neighbors are observed. By default, all individuals are focal.
#' @param model model used for estimation. Any of 'penalized_regression', 'linear_regression', 'random_forest';
#' @param family   family of the outcome, either gaussian or binary;
#' @param params list with additional parameters (see references)
#' @return predictions: predicted effect on each individual;
#' @return welfare: point estimate of the welfare;
#' @export
#'


NetworkWelfare <- function(policies, Y, X, Z, Z_int, D, W_0, W_0_train, cost = 0, focal_units = rep(1, length(Y)),
                             model = 'penalized_regression', family = 'gaussian', params = list()){

  ## handle data missingness

  stopifnot("All values of D must be non-missing"= sum(is.na(D)) == 0)

  stopifnot("All values of X must be non-missing"= sum(apply(X, 1, function(x) sum(is.na(x)) > 0)) == 0)

  missing_Z = apply(Z, 1, function(x) sum(is.na(x)) > 0)
  if(sum(missing_Z) > 0) {
    focal_units[missing_Z] = 0
    print("Some Z are missing! Units with missing Z are set to non-focal units; missing values are imputed with column means")
    for (i in 1:dim(Z)[2]) {
      Z[is.na(Z[,i]), i] <- mean(Z[,i], na.rm = TRUE)
    }
    stopifnot(sum(apply(Z, 1, function(x) sum(is.na(x)) > 0))==0)
  }

  missing_Z_int = apply(Z_int, 1, function(x) sum(is.na(x)) > 0)
  if(sum(missing_Z_int) > 0) {
    focal_units[missing_Z_int] = 0
    print("Some Z_int are missing! Units with missing Z_int are set to non-focal units; missing values are imputed with column means")
    for (i in 1:dim(Z_int)[2]) {
      Z_int[is.na(Z_int[,i]), i] <- mean(Z_int[,i], na.rm = TRUE)
    }
    stopifnot(sum(apply(Z_int, 1, function(x) sum(is.na(x)) > 0))==0)
  }

  missing_W_0 = apply(W_0, 1, function(x) sum(is.na(x)) > 0)
  if(sum(missing_W_0) > 0) {
    focal_units[missing_W_0] = 0
    print("Some W are missing! Units with missing W are set to non-focal units; missing values are imputed with zero")
    for (i in 1:dim(W_0)[2]) {
      W_0[is.na(W_0[,i]), i] <- 0
    }
    stopifnot(sum(apply(W_0, 1, function(x) sum(is.na(x)) > 0))==0)
  }

  missing_W_0_train = apply(W_0_train, 1, function(x) sum(is.na(x)) > 0)
  if(sum(missing_W_0_train) > 0) {
    focal_units[missing_W_0_train] = 0
    print("Some W_0_train are missing! Units with missing W_0_train are set to non-focal units; missing values are imputed with zero")
    for (i in 1:dim(W_0_train)[2]) {
      W_0_train[is.na(W_0_train[,i]), i] <- 0
    }
    stopifnot(sum(apply(W_0_train, 1, function(x) sum(is.na(x)) > 0))==0)
  }

  missing_Y = sapply(Y, is.na)
  if(sum(missing_Y) > 0) focal_units[missing_Y] = 0

  ## setup

  my_params = create_params(params)

  ## invoke the optimal_network_linear_rule to return the conditional mean function only
  my_params$return_m1 = T
  m1 = optimal_network_linear_rule(Y,X,Z, Z_int, D,W_0,  impose_max_treated = NA,
                                                method =  'MILP', m1 = my_params$m1,
                                                model = model, cost = cost, cut_high_neighb = my_params$cut_high_neighb,
                                   tolerance_constraint1 = my_params$tolerance_constraint1,
                                   tolerance_constraint2 = my_params$tolerance_constraint1,
                                   B = my_params$B, low_b = my_params$low_b,
                                   params = my_params$params_gurobi,
                                   monotonicity_c = my_params$monotonicity_c,
                                   maxtime_partitions = my_params$maxtime_partitions,
                                   length_small_samples = my_params$length_small_samples,
                                   doubly_robust = my_params$doubly_robust,
                                   penalized= my_params$penalized,
                                   num_strata = my_params$num_strata,
                                   trimming = my_params$trimming,
                                   family = family,
                                   maximum_time = maximum_time,
                                   additional_param_constraints = my_params$additional_param_constraints,
                                   model_return = my_params$model_return,
                                   W_0_train = W_0_train,
                                   names_rf = my_params$names_rf,
                                   print_checking = my_params$print_checking,
                                   focal_units = focal_units,
                                   return_m1 = my_params$return_m1,
                                   bound_outcomes = my_params$bound_outcomes,
                                   trimmed_subpopulation = my_params$trimmed_subpopulation,
                                   min_treated = my_params$min_treated,
                                   Z_int_neighbors = Z_int,
                                   partitions  = my_params$partitions,
                                   numcores = my_params$numcores,
                                   cross_fitting = my_params$cross_fitting,
                                   slack = my_params$slack,
                                   n_folds = my_params$n_folds,
                                   verbose = my_params$verbose)
   welfare = suppressWarnings(compute_welfare (list(policies), Y,D,W_0,
                            Z, Z_int, m1 = m1$model,
                            family = family,
                            penalized = my_params$penalized,
                            num_strata = my_params$num_strata,
                            trimming = my_params$trimming,
                            doubly_robust = my_params$doubly_robust,
                            names_rf = m1$names_rf,
                            W_0_train = W_0_train,
                            focal_units = focal_units,
                            cross_fitting = my_params$cross_fitting,
                            bound_outcomes = my_params$bound_outcomes,
                            trimmed_subpopulation = my_params$trimmed_subpopulation,
                            cross_fitting_vector = m1$cross_fitting_vector,
                            cost = cost))
  return(welfare)
}



#' Predict the policy for new observations
#'
#' This function predicts the policy for new observations
#'
#' @param opt solutions to NetworkTargeting
#' @param X variables for targeting;
#' @param params list with additional parameters (see references)
#' @return predictions: predicted effect on each individual;
#' @export


PredictPolicies <- function(opt, X){
  policies <- apply(X, 1, function(x) ifelse(c(1,x)%*%opt$final_coef > 0, 1, 0))
  return(policies)
}

