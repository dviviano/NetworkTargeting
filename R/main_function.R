
## Input:
## Y: Outcome
## Note: NA values in Y for those units whose neighbors are unknown but who can be neighbors of focal points
## X: variables for targeting
## Z: variables for conditional mean and propensity score estimation
## Z_int : variables used for interaction in the conditiona mean function
## D: treatment assignment
## W_0: adjacency matrix
## m1: fitted model for conditional mean(if NA - default estimation is performed)

## The model must take Z, Z_int,num_treated_neighbors/number of neighbors*Z_int,  num_treated_neighbors*Di/num_neighb, num_treated_neighbrs/num_neighb, Di

## cost: cost of treatment (if not included in conditional mean function)
## cut_high_neighb: remove effects on units with high degree from computations for increasing speed
## tolerance constraints: epsilon constraints in strict inequality in gurobi formulation
## B, low_b: upper and lower bound on coefficients
## params: additional parameters to be passted to gurobi (default: 3000 maxtime)
## impose max treated: if numeric maximum number of treated uit
## monotonicity_c: monotonicity constraints on the computation of the maximum score
## print checking: print that constraints are respected from the formulation
## doubly robust: use also propensity score adjustment
## penalized: penalized estimation of proensity score (?)
## num_strata: number of treatment status for propensity score estimation. If NA each number of neighbor correspond to a different treatment status
## trimming: trimming value for propensity score
## family: family of the outcome, either gaussian or binary, denote the support of the outcome

## Method: MILP (exact formulation) or MILP_graph_clustering (approximate)
## maxtime_partitions: maximum time to find partitions
## length_small_samples: size of subsamples when graph clustering algorithm is implemented.
## additional_param_constraints: if not NA, a matrix with nrow = dim(X)[2] + 1 indicating UB (left column) and LB (right column) for each parameter
## W_0_train if not missing it sets the adjacency matrix for the training of the model different from the one used for the policy estimation
## names_rf: column names of the data frame if random forest model is passed to the function (if m1 = NA and randomForest::randomForest is specified the function automatically sets names_rf correctly to the column names)
## focal_units: indexes of units used for the optimization (1 for focal 0 otherwise) - default all
## return_m1: if true the function only returns the estimated conditional mean function
## bound_outcomes: upper and lower bounds on the outcome applied to the estimated conditional mean function hat m
## Z_int_neighbors: interaction terms betwee number of neighbors and individual covariates
## min_treated: minimum number of treated units
## partitions: if selected MILP_graph_clustering, indexes denote a list with each element containing the index of each cluster
## cross_fitting: boolean note: m1 if not NA has to be a list with models for each observation
##slack: numeric indicating the size of the second group equal half of the sample +- slack parameter
## n_folds: number of folds for cross fitting
## share_treated_neighbors: share of treated neighbors (if missing obtained through W_0_train and D)
## pass_num_neighbors: total number of neighbors used in the denominator for the regression. If NA the neighbors from W_0 are used

optimal_network_linear_rule <- function(Y,X,Z, Z_int, D, W_0,  impose_max_treated = NA, method = 'MILP',
                                        m1 = NA, model = 'linear_regression', cost = 0, cut_high_neighb = F,
                                        tolerance_constraint1 = 10**(-6),
                                        tolerance_constraint2 = 10**(-6), B = 1, low_b = - B,params = NA,
                                        monotonicity_c = F,maxtime_partitions = 100, length_small_samples = 50,
                                        doubly_robust = F, penalized= F, num_strata = 4, trimming = 0.01,
                                        family = 'gaussian', maximum_time = 100,
                                        additional_param_constraints = NA, model_return = F,
                                        W_0_train = W_0, names_rf = NA, print_checking = F,
                                        focal_units = rep(1, length(Y)), return_m1 = F,
                                        bound_outcomes = c(-Inf, Inf),
                                        trimmed_subpopulation = F, min_treated = 0,
                                        Z_int_neighbors = Z_int,
                                        partitions  = NA,
                                        numcores = 1, cross_fitting = T, slack = 10, n_folds= 5,
                                        share_treated_neighbors = NA,
                                        pass_num_neighbors = NA, verbose = F){

  cross_fitting_vector = NA

  ## Set the connections between non-focal points equal to zero
  ## (these are not used for estimation purposes and they reduce memory requirements)

  if(sum(focal_units) < length(Y)) {
    W_0[focal_units == 0, focal_units == 0] <- 0
    W_0_train[focal_units == 0, focal_units == 0] <- 0
  }

  p <- dim(X)[2]

  if(method == 'MILP' & dim(W_0)[1] > 500) warning('Running time may take long, try to implement MILP_graph_clustering for moderate size or the surrogate method for large size to reduce running time')

  if (is.na(m1)[1]) {

    if (cross_fitting == F) {

      my_model = compute_conditional_mean_function(Y, Z, Z_int, D, Z_int_neighbors, W_0_train, model, family, cross_fitting = F, cross_fitting_vector = NA,
                                                   focal_units = focal_units, share_treated_neighbors = share_treated_neighbors,
                                                   pass_num_neighbors = pass_num_neighbors)
      m1 = my_model[[1]]
      names_rf = my_model[[2]]
      cross_fitting_vector = NA

      if (return_m1) return(list(model = m1, cross_fitting_vector = cross_fitting_vector, names_rf = names_rf))

    } else {
      my_cut = opposite_min_cut_MILP_1_and_2_degree(W = W_0_train, n_folds, slack = 5, focal_units = focal_units, verbose = verbose)
      my_model = compute_conditional_mean_function(Y, Z, Z_int, D, Z_int_neighbors, W_0_train, model, family, cross_fitting = T, cross_fitting_vector = my_cut,
                                                   focal_units = focal_units, share_treated_neighbors,
                                                   pass_num_neighbors = pass_num_neighbors)
      m1 = my_model[[1]]
      names_rf = my_model[[2]]
      ## Cross fitting vector for the propensity score to be passed to the functions below
      cross_fitting_vector = my_cut$cut1

      if(return_m1) return(list(model = m1, cross_fitting_vector = cross_fitting_vector, names_rf = names_rf))

      if(sum(missing_X) > 0) {
        cross_fitting_vector = cross_fitting_vector[-missing_X]
        m1 = m1[-missing_X]
      }
    }
  } else {
    ## Not missing m1
    if(cross_fitting == F) {
      cross_fitting_vector = NA
    } else {
      my_cut = opposite_min_cut_MILP_1_and_2_degree(W = W_0_train, n_folds, slack = 5, focal_units = focal_units, verbose = verbose)
      ## Cross fitting vector for the propensity score to be passed to the functions below
      cross_fitting_vector = my_cut$cut1
      if(sum(missing_X) > 0) {
        cross_fitting_vector = cross_fitting_vector[-missing_X]
        m1 = m1[-missing_X]
      }
    }
  }

  ## Remove missing values for the subsequent estimation
  # if(sum(missing_X) > 0){
  #   Y = Y[-missing_X]
  #   X = X[-missing_X, ]
  #   Z = Z[-missing_X, ]
  #   Z_int = Z_int[-missing_X, ]
  #   D = D[-missing_X]
  #   W_0 = W_0[-missing_X, -missing_X]
  #   focal_units = focal_units[-missing_X]
  #   Z_int_neighbors = Z_int_neighbors[-missing_X,]
  #   pass_num_neighbors = pass_num_neighbors[-missing_X]
  # }

  if(method == 'MILP') {
    sol <- suppressWarnings(max_score_NEWM(Y,X,Z, Z_int, D,W_0,
                                      m1, cost, cut_high_neighb,
                                      tolerance_constraint1,
                                      tolerance_constraint2, B, low_b,params ,
                                      impose_max_treated, monotonicity_c, doubly_robust = doubly_robust, penalized= penalized, num_strata = num_strata,
                          trimming = trimming, type = family, maximum_time = maximum_time, additional_param_constraints = additional_param_constraints,model_return = model_return, W_0_train = W_0_train,
                          names_rf  = names_rf, print_checking = print_checking, focal_units = focal_units, bound_outcomes = bound_outcomes,
                          min_treated = min_treated,
                          Z_int_neighbors = Z_int_neighbors, numcores = numcores,
                          cross_fitting = cross_fitting, cross_fitting_vector = cross_fitting_vector,
                          trimmed_subpopulation = trimmed_subpopulation,
                          pass_num_neighbors = pass_num_neighbors, verbose = verbose))

    if (model_return) {
      return(sol)
    } else {
      my_sol <- sol[[1]]
      final_coef <- my_sol$x[(length(my_sol$x) - p):length(my_sol$x)]
      objective_function <- sol$objval

      treatment_assignments <- apply(cbind(1, X), 1, function(x) ifelse(x%*%final_coef >= 0, 1, 0))

      ## Diagnostic check
      checks <-sol[[2]]

        ## Sanity check for stable solutions
        discrepancies <- max(checks[2], checks[3], checks[4])
        if(discrepancies > 0.2) warning(paste0('Unstable MILP solution. Increase the tolerance constraint. Discrepancies with solutions of:', discrepancies))

      if (verbose) print(paste0('Number of Treated is:', sum(treatment_assignments)))
      if (verbose) print(paste0('Objective:', my_sol$objval/sum(focal_units)))

      return(list(policies = treatment_assignments,
                  all_gurobi_output = my_sol$x, model = m1, final_coef = final_coef,
                  obj = my_sol$objval/sum(focal_units), miscellaneus = my_sol, checks = checks,
                  names_rf = names_rf,
                  cross_fitting_vector = cross_fitting_vector))
    }

  } else if (method == 'MILP_graph_clustering') {

    ## Estimate the propensity score before the splitting
    reg_prop = estimate_individual_propensity_score(D, Z, penalized, reg = NA, cross_fitting, cross_fitting_vector)

    nn <- dim(W_0)[1]
    num_cut <- floor(nn/length_small_samples)
    ## num_cut denote the number of partitions
    if (is.na(partitions)) {
      partitions <- recursive_min_cut(W_0, maxtime_partitions, params, 1, floor(log(num_cut)/log(2)), c(1:nn), numcores = numcores)
    }
    sol <- suppressWarnings(recursive_max_score_NEWM(Y,X,Z, Z_int, D,W_0,
                                    m1, cost, cut_high_neighb,
                                    tolerance_constraint1,
                                    tolerance_constraint2, B, low_b,params ,
                                    impose_max_treated, monotonicity_c, index_list = partitions,
                                    doubly_robust = doubly_robust, penalized= penalized, num_strata = num_strata,
                                    trimming = trimming, type = family, maximum_time = maximum_time, additional_param_constraints = additional_param_constraints, W_0_train = W_0_train,  names_rf  = names_rf
                                    , print_checking = print_checking, focal_units = focal_units, bound_outcomes = bound_outcomes,
                                    min_treated = min_treated,
                                    Z_int_neighbors = Z_int_neighbors, numcores = numcores, reg = reg_prop,
                                    cross_fitting = cross_fitting, cross_fitting_vector = cross_fitting_vector,
                                    trimmed_subpopulation = trimmed_subpopulation,
                                    pass_num_neighbors = pass_num_neighbors, verbose = verbose))
    ## Return list of coefficients
    sol <- as.data.frame(sol)
    names(sol)[1:(dim(sol)[2] - 2)] <- paste0('Coef_Variable', c(1:(dim(sol)[2] - 2 )))

    names(sol)[dim(sol)[2] - 1] <- 'Objective'

    names(sol)[dim(sol)[2]] <- 'Number of Elements in the cluster'
    rownames(sol) <- paste0('V', c(1:dim(sol)[1]))
    solution_aggregate <- compute_aggregate_solution(sol, X, indexes =  partitions , impose_max_treated = impose_max_treated,
                                                     maximum_time = maximum_time, params, numcores = numcores, verbose = verbose)


    final_coef <- solution_aggregate$all_res[[2]]
    objective_function <- solution_aggregate$objective

    treatment_assignments <- apply(cbind(1, X), 1, function(x) ifelse(x%*%final_coef >= 0, 1, 0))
    if (verbose) print(paste0('Number of Treated is:', sum(treatment_assignments)))
    if (verbose) print(paste0('Percentage of correct classifications is:',  objective_function))


    return(list(policies = treatment_assignments, all_gurobi_output =  solution_aggregate ,
                model = m1, final_coef = final_coef, obj =  objective_function,
                all_solutions = sol, names_rf = names_rf,
                cross_fitting_vector = cross_fitting_vector))


  }
}

## Returns the same as before but also passes the best threshold for the definition of objective function
## To avoid instability complications with strict inequality
## The threshold maximizes in sample welfare

Network_EWM_threshold <- function(Y,X,Z, Z_int, D,W_0,  impose_max_treated = NA, method = 'MILP',
                                  m1 = NA, model = 'linear_regression', cost = 0, cut_high_neighb = F,
                                  tolerance_constraint1 = 10**(-6),
                                  tolerance_constraint2 = 10**(-6), B = 1, low_b = - B,params = NA,
                                  monotonicity_c = F,maxtime_partitions = 100, length_small_samples = 50,
                                  doubly_robust = F, penalized= F, num_strata = 4, trimming = 0.01,
                                  family = 'gaussian', maximum_time = 100,
                                  additional_param_constraints = NA, model_return = F,
                                  W_0_train = W_0, names_rf = NA, print_checking = F,
                                  focal_units = rep(1, length(Y)), return_m1 = F,
                                  bound_outcomes = c(-Inf, Inf),
                                  trimmed_subpopulation = F, min_treated = 0,
                                  Z_int_neighbors = Z_int,
                                  partitions  = NA,
                                  numcores = 1, cross_fitting = T, slack = 10, n_folds= 5,
                                  share_treated_neighbors = NA,
                                  pass_num_neighbors = NA) {

  model_00 = optimal_network_linear_rule(Y,X,Z, Z_int, D,W_0,  impose_max_treated, method,
                              m1, model, cost, cut_high_neighb,
                              tolerance_constraint1,
                              tolerance_constraint2, B, low_b,params,
                              monotonicity_c,maxtime_partitions, length_small_samples,
                              doubly_robust, penalized, num_strata, trimming,
                              family, maximum_time,
                              additional_param_constraints, model_return,
                              W_0_train, names_rf, print_checking,
                              focal_units, return_m1,
                              bound_outcomes,
                              trimmed_subpopulation, min_treated,
                              Z_int_neighbors,
                              partitions,
                              numcores, cross_fitting, slack, n_folds,
                              share_treated_neighbors,
                              pass_num_neighbors )


  my_seq = seq(from = -20, to = 20, length = 200)
  my_seq = sort(my_seq, decreasing= T)
  policies = list()
  acc = 1
  for(j in my_seq){
    policies[[acc]] = sapply(cbind(1, X)%*%model_00$final_coef, function(x) ifelse(x >= j * tolerance_constraint1, 1, 0))
    acc = acc  + 1
  }

 # policy1 =  sapply(cbind(1, X)%*%c(-1, -1, 1, 1, -1), function(x) ifelse(x >= 0, 1, 0))
 # policy2 = sapply(cbind(1, X)%*%c(-1, 1, 1, 1, -1), function(x) ifelse(x >= 0, 1, 0))

  welf1_insample = compute_welfare(my_policies = policies, Y = Y, D = D,
                                   W_0 = W_0,Z = Z, Z_int = Z_int,
                                   Z_int_neighbors = Z_int_neighbors,
                                   m1 = model_00$model, doubly_robust = doubly_robust, num_strata = num_strata,
                                   trimming = trimming, W_0_train = W_0_train, names_rf = model_00$names_rf,
                                   focal_units = focal_units, family = family,
                                   pass_num_neighbors = pass_num_neighbors, cross_fitting = cross_fitting,
                                   bound_outcomes = bound_outcomes, reg = NA, trimmed_subpopulation = trimmed_subpopulation, penalized = penalized,
                                   cross_fitting_vector = model_00$cross_fitting_vector, cost = cost)


  threshold = my_seq * tolerance_constraint1
  ## This deals with parameters of order tolerance_constraint in the presence of instability
  ## Where results may differ by an epsilon difference

  #if(which.max(welf1_insample$welfare) > 1){  model_00$best_threshold = threshold[which.max(welf1_insample$welfare) - 1]
  #if(1 - welf1_insample$welfare[which.max(welf1_insample$welfare) - 1]/welf1_insample$welfare[which.max(welf1_insample$welfare)] > 0.01) model_00$best_threshold = threshold[which.max(welf1_insample$welfare)]
  #} else {
  #  model_00$best_threshold = threshold[which.max(welf1_insample$welfare)]
  #}
  model_00$best_threshold = threshold[which.max(welf1_insample$welfare)]
  model_00$welfare_computations_for_threshold = welf1_insample$welfare
  return(model_00)
}

# X here is a matrix with binary columns
NEWM_for_binary_groups <- function(Y,X,Z, Z_int, D,W_0,  impose_max_treated = NA,
                                   m1 = NA, model = 'linear_regression', cost = 0,
                                   doubly_robust = F, penalized= F, num_strata = 4, trimming = 0.01,
                                   family = 'gaussian', maximum_time = 100,
                                   additional_param_constraints = NA,
                                   W_0_train = W_0,
                                   focal_units = rep(1, length(Y)),
                                   bound_outcomes = c(-Inf, Inf),
                                   trimmed_subpopulation = F, min_treated = 0,
                                   Z_int_neighbors = Z_int,
                                   numcores = 1, cross_fitting = T, slack = 10, n_folds= 5,
                                   share_treated_neighbors = NA,
                                   pass_num_neighbors = NA, ...) {


 n = dim(X)[1]
 unique_values = unique(X)
 p = dim(unique_values)[1]
 policies = list()
 policy_coefs = list()
 acc = 1
 if(is.na(impose_max_treated)) impose_max_treated = n
 ## Generate all combinations of which policies you treat
 my_combinations = list()
 for(j in 1:p){
   my_combinations[[j]] = t(combn(c(1:p), j))
 }



comined_results = foreach(j = 1:length(my_combinations), .combine = append)%dopar%{

  my_list = my_combinations[[j]]
  acc = 1
  for(g in 1:dim(my_list)[1]){

    my_policy = apply(sapply(my_list[g,], function(y) apply(X, 1, function(k) all(k == unique_values[y,]))), 1,
               function(z) sum(z > 0))

    screening = sum(my_policy) <= impose_max_treated & sum(my_policy) >= min_treated
    if(screening){
    policies[[acc]] = my_policy
    policy_coefs[[acc]] = unique_values[my_list[g,],]
    acc = acc + 1
    }


  }
  list(policies, policy_coefs)

 }

all_policies = list()
all_policy_coefs = list()
for(j in seq(from = 1, to = 2 * length(my_combinations), by = 2)){
  all_policies = append(all_policies, comined_results[[j]])
  all_policy_coefs = append(all_policy_coefs, comined_results[[j + 1]])
}
 all_policies[[acc + 1]] = rep(0, n)
 all_policy_coefs[[acc + 1]] = rep(0, dim(unique_values)[2])
if(is.na(m1)[1]){
 model_00 = optimal_network_linear_rule(Y,X,Z, Z_int, D,W_0,  impose_max_treated = NA, method = 'MILP',
                                        m1 = NA, model, cost, cut_high_neighb = F,
                                        tolerance_constraint1 = 10**(-6),
                                        tolerance_constraint2 = 10**(-6), B = 1, low_b = - B,params = NA,
                                        monotonicity_c = F,maxtime_partitions = 100, length_small_samples = 50,
                                        doubly_robust, penalized, num_strata, trimming,
                                        family, maximum_time,
                                        additional_param_constraints, model_return = T,
                                        W_0_train, names_rf = NA, print_checking = F,
                                        focal_units, return_m1 = T,
                                        bound_outcomes,
                                        trimmed_subpopulation, min_treated,
                                        Z_int_neighbors,
                                        partitions  = NA,
                                        numcores, cross_fitting, slack, n_folds,
                                        share_treated_neighbors,
                                        pass_num_neighbors)
} else {
  model_00 = m1
}

 welf1_insample = compute_welfare(my_policies = all_policies, Y = Y, D = D,
                                  W_0 = W_0,Z = Z, Z_int = Z_int,
                                  Z_int_neighbors = Z_int_neighbors,
                                  m1 = model_00$model, doubly_robust = doubly_robust, num_strata = num_strata,
                                  trimming = trimming, W_0_train = W_0_train, names_rf = model_00$names_rf,
                                  focal_units = focal_units, family = family,
                                  pass_num_neighbors = pass_num_neighbors, cross_fitting = cross_fitting,
                                  bound_outcomes = bound_outcomes, reg = NA, trimmed_subpopulation = trimmed_subpopulation, penalized = penalized,
                                  cross_fitting_vector = model_00$cross_fitting_vector, cost = cost,
                                  parallel = T)


 best_policy = which.max(welf1_insample$welfare)
 groups_to_treat = all_policy_coefs[[best_policy]]
 final_policy = all_policies[[best_policy]]
 return(list(groups_to_treat = groups_to_treat, final_policy = final_policy,
             all_results = list(welf1_insample = welf1_insample,
                                model = model_00 ,
                                all_policies_list = all_policies,
                                all_groups_list = all_policy_coefs
                                ) ))
}
