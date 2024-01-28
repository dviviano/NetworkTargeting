

## Input:
## Y: Outcome (NAs: these individuals are the neighbors of the focal units)
## X: variables for targeting
## Z: variables for conditional mean and propensity score estimation
## Z_int : variables used for interaction in the conditiona mean function
## D: treatment assignment
## W_0: adjacency matrix
## m1: fitted model for conditional mean
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
## type: either continuous or binary, denote the support of the outcome
## W_0_train: if specified a different adjacency is used for the training of the model
## names_rf: names of the data frame for prediction (used only for the random forest method)
## focal_units: indexes of units used for the optimization (1 for focal 0 otherwise)
## bound_outcomes: upper and lower bounds on the outcome applied to the estimated conditional mean function hat m
## min_treated: minimum number of treated units
## Z_int_neighbors: interactions with number of treated neighbors
## reg: regression function for propensity score (if NA is estimated)
## cross_fitting: boolean, if True m1 is a list with a model for each unit
## cross_fitting_vector is a vector of 0/1 for each group for cross-fitting
## pass_num_neighbors: # of total neighbors in the denominator of the regression. IF NA use W_0

max_score_NEWM <- function(Y,X,Z, Z_int = NA, D,W_0,
                           m1 = NA, cost = 0, cut_high_neighb = F,
                           tolerance_constraint1 = 10**(-6),
                           tolerance_constraint2 = 10**(-6), B = 1, low_b = - 1,params = NA,
                           impose_max_treated = NA, monotonicity_c = F, print_checking = F, doubly_robust = F,
                           penalized= F, num_strata = 4, trimming = 0.01, type, maximum_time = 1000,
                           additional_param_constraints = NA,model_return = F,
                           W_0_train = W_0, names_rf = NA, focal_units = rep(1, length(Y)),
                           bound_outcomes = c(-Inf, Inf), trimmed_subpopulation = F,
                           min_treated = 0,
                           Z_int_neighbors = Z_int_neighbors,
                           numcores = 1, reg = NA,
                           cross_fitting = F, cross_fitting_vector = NA,
                           pass_num_neighbors = NA, verbose = F){

  final_g0_noneighb <- NA
  final_gs_noneighb <- NA

  if (doubly_robust == F) {

    gs <- compute_gs(W_0, Z,Z_int,m1, type, names_rf = names_rf, focal_units = focal_units,
                     Z_int_neighbors = Z_int_neighbors, cross_fitting = cross_fitting,
                     pass_num_neighbors = pass_num_neighbors, cost = cost)
    ## Apply the bound on the effects if outcome is a bounded random variable with known upper/lower bounds
    final_gs <- gs[[1]]
    final_g0 <- gs[[2]]
    final_gs <- sapply(final_gs, function(x) min(max(x, bound_outcomes[1]), bound_outcomes[2]))
    final_g0 <- sapply(final_g0, function(x) min(max(x, bound_outcomes[1]), bound_outcomes[2]))

  }

  ## Note: missing values passing to DR is not needed since these are computed within the function
  no_neighb <- which(apply(W_0, 1, sum) == 0 & apply(W_0, 2, sum) == 0)
  any_neighb <- which(apply(W_0, 1, sum) > 0 & apply(W_0, 2, sum) > 0)

  if (doubly_robust & length(any_neighb) > 0) {

    gs <- compute_DR_effects(Y,D,W_0, Z, Z_int, also_reg_adj = T, m1, type, penalized,
                             num_strata, trimming, W_0_train = W_0_train, names_rf = names_rf, focal_units = focal_units,
                             trimmed_subpopulation = trimmed_subpopulation,
                             Z_int_neighbors = Z_int_neighbors, reg = reg,
                             cross_fitting = cross_fitting, cross_fitting_vector = cross_fitting_vector,
                             pass_num_neighbors = pass_num_neighbors, cost = cost)

    final_gs <- gs[[1]]
    final_g0 <- gs[[2]]

  }

  ## Consider no_neighb on both columns and rows for asymmetric adjacency matrices

  if (length(no_neighb) > 1) {

    nas_Y <- which(is.na(Y[no_neighb]))
    if(length(nas_Y) == 0) nas_Y = NA

    if(doubly_robust == F) {

      if(cross_fitting) my_m1 = m1[no_neighb]

      gs <- compute_gs(diag(1, length(no_neighb)), Z[no_neighb,],
                       Z_int[no_neighb,],m1, type, no_atoms=F,
                       names_rf = names_rf, focal_units = focal_units[no_neighb],
                       Z_int_neighbors = Z_int_neighbors[no_neighb,], cross_fitting = cross_fitting, cost = cost)

      final_gs_noneighb  <- gs[[1]]
      final_g0_noneighb <- gs[[2]]
      final_gs_noneighb <- sapply(final_gs_noneighb, function(x) min(max(x, bound_outcomes[1]), bound_outcomes[2]))
      final_g0_noneighb <- sapply(final_g0_noneighb, function(x) min(max(x, bound_outcomes[1]), bound_outcomes[2]))
    }

    if(doubly_robust) {

      gs <- compute_DR_effects(Y,D,W_0, Z, Z_int, also_reg_adj = T, m1, type, penalized, num_strata, trimming, no_atoms = F,
                                                W_0_train = W_0_train, names_rf = names_rf, focal_units = focal_units,
                                                bound_outcomes = bound_outcomes, trimmed_subpopulation = trimmed_subpopulation,
                                                Z_int_neighbors = Z_int_neighbors, reg = reg,
                               cross_fitting = cross_fitting, cross_fitting_vector = cross_fitting_vector, cost = cost)

      final_gs_noneighb  <- gs[[1]]
      final_g0_noneighb <- gs[[2]]

      if (length(any_neighb) == 0) {
        print("Not a single person with a neighbor!")
        final_gs_noneighb  <- NA
        final_g0_noneighb <- NA
        final_gs  <- gs[[1]]
        final_g0 <- gs[[2]]
      }

    }

  }

  ## cost is already in the compute gs functions
  constraints_final <- compute_constraints_function(Y,X,Z, Z_int,D,W_0,
                                                    final_gs, final_g0, final_gs_noneighb,
                                                    final_g0_noneighb, m1, cost = 0, cut_high_neighb, B, tolerance_constraint1,
                                                    tolerance_constraint2, monotonicity_c = monotonicity_c, direction = '>=',
                                                    Z_int_neighbors = Z_int_neighbors)
  n <- length(Y)
  p <- dim(X)[2]

  my_result <- max_score_NEWM_helper(
    constraint_matrix = constraints_final$constraint_matrix, XX = constraints_final$XX, impose_max_treated,
    b = constraints_final$b, model = constraints_final$model, y_s_position = constraints_final$y_s_position,
    no_neighb = constraints_final$no_neighb, any_neighb = constraints_final$any_neighb,
    final_gs_noneighb = constraints_final$final_gs_noneigh,
    my_h = constraints_final$my_h,
    c = constraints_final$c, tolerance_constraint1,
    tolerance_constraint2, B = B, low_b = low_b, params = params,
    constant = - sum(final_g0, na.rm = T) + sum(final_g0_noneighb, na.rm = T), cost = cost,
    constraint_matrix_no_n = constraints_final$constraint_matrix_no_n, maximum_time = maximum_time,
    additional_param_constraints = additional_param_constraints, model_return = model_return,
    min_treated= min_treated, numcores = numcores, focal_units = focal_units, verbose = verbose)

  if (model_return) {

    return(my_result)

  } else {

    checks <- model_checking(result = my_result,XX = constraints_final$XX, Z,
                             W_0, y_s_position = constraints_final$y_s_position,
                             positions = constraints_final$positions,
                             beta_hat = my_result$x[c((length(my_result$x) - p):length(my_result$x))],
                             final_gs, final_g0,cc = constraints_final$c, to_print = print_checking,
                             focal_units = focal_units, final_gs_noneighb = constraints_final$final_gs_noneigh,
                             final_g0_noneighb = final_g0_noneighb)

    obj1 <- my_result$objval

    if (monotonicity_c) {
      ## Do with opposite direction
      constraints_final <- compute_constraints_function(Y,X,Z, Z_int,D,W_0,
                                                        final_gs, final_g0, final_gs_noneighb,
                                                        final_g0_noneighb, m1, cost = 0, cut_high_neighb, B, tolerance_constraint1,
                                                        tolerance_constraint2, monotonicity_c = monotonicity_c, direction = '<=',
                                                        Z_int_neighbors = Z_int_neighbors)

      my_result2 <- max_score_NEWM_helper(constraints_final$constraint_matrix, constraints_final$XX, impose_max_treated,
                                          constraints_final$b, constraints_final$model, constraints_final$y_s_position,
                                          constraints_final$no_neighb, constraints_final$final_gs_noneighb, constraints_final$my_h,
                                          constraints_final$c, tolerance_constraint1,
                                          tolerance_constraint2, B, low_b, params,  cost = cost,
                                          constraint_matrix_no_n = constraints_final$constraint_matrix_no_n, maximum_time = maximum_time,
                                          additional_param_constraints = additional_param_constraints,
                                          min_treated = min_treated, numcores = numcores, focal_units = focal_units, verbose = verbose)


      checks <- tryCatch(model_checking(my_result2,constraints_final$XX, Z, W_0,
                                        constraints_final$y_s_position, constraints_final$positions,
                                        my_result$x[c((length(results) - p):length(results))],
                                        final_gs, final_g0,constraints_final$c,
                                        to_print = print_checking, focal_units = focal_units), error = function(e) {NA})
      obj2 <- my_result2$objval
      if (obj2 >= obj1) {
        my_result <- my_result2
      }
    }

    return(list(my_result, checks))
  }

}

## Run the above function recursevely
## index list: list of indexes over which run the NEWM maximum score

recursive_max_score_NEWM <- function(Y,X,Z, Z_int, D,W_0,
                                     m1, cost = 0, cut_high_neighb = F,
                                     tolerance_constraint1 = 10**(-6),
                                     tolerance_constraint2 = 10**(-6), B = 1, low_b = - B,params = NA,
                                     impose_max_treated = NA, monotonicity_c = F,
                                     print_checking = F, index_list,
                                     doubly_robust = F, penalized= F, num_strata = 4,
                                     trimming = 0.01, type,
                                     maximum_time = 1000,
                                     additional_param_constraints = NA, W_0_train = W_0, names_rf = NA,
                                     focal_units = rep(1, length(Y)),
                                     bound_outcomes = c(-Inf, Inf), trimmed_subpopulation = F,
                                     min_treated = min_treated,
                                     Z_int_neighbors = Z_int,
                                     numcores = 1, reg = NA,  cross_fitting = F, cross_fitting_vector = NA,
                                     pass_num_neighbors = NA, verbose = F){


  if(is.list(index_list) == F){
    p <- dim(X)[2]
    if(cross_fitting) m1 = m1[index_list]
    if(length(reg) > 1) reg = reg[index_list]
    sol <- max_score_NEWM(
      Y[index_list],X[index_list,],Z[index_list,], Z_int[index_list,], D[index_list],W_0[index_list, index_list],
                          m1, cost, cut_high_neighb,
                          tolerance_constraint1,
                          tolerance_constraint2, B, low_b,params ,
                          ## Rescale appropriately the number of treated units
                          floor(impose_max_treated/(dim(W_0)[1]))*length(index_list), monotonicity_c,
                          doubly_robust = doubly_robust, penalized= penalized, num_strata = num_strata,
                          trimming = trimming, type = type, maximum_time = maximum_time,
                          additional_param_constraints = additional_param_constraints, W_0_train = W_0_train[index_list, index_list],
                          names_rf = names_rf,
                         print_checking = print_checking, focal_units = focal_units[index_list],
                         bound_outcomes = bound_outcomes, trimmed_subpopulation = trimmed_subpopulation,
                         min_treated = min_treated,
                         Z_int_neighbors = Z_int_neighbors[index_list,], numcores = numcores, reg = reg,
                         cross_fitting = cross_fitting, cross_fitting_vector = cross_fitting_vector[index_list],
      pass_num_neighbors = pass_num_neighbors[index_list], verbose = verbose)
    sol <- sol[[1]]
    final_coef <- sol$x[(length(sol$x) - p):length(sol$x)]
    objective_function <- sol$objval
    return(c(final_coef, objective_function, length(index_list)))
  } else {

    k = 1
    for(j in index_list){

      rec1 <- recursive_max_score_NEWM(Y,X,Z, Z_int, D,W_0,
                                       m1, cost, cut_high_neighb,
                                       tolerance_constraint1,
                                       tolerance_constraint2, B, low_b,params,
                                       impose_max_treated, monotonicity_c, print_checking, index_list = j,
                                       doubly_robust = doubly_robust, penalized= penalized, num_strata = num_strata,
                                       trimming = trimming, type = type, maximum_time = maximum_time,
                                       additional_param_constraints= additional_param_constraints,
                                       W_0_train = W_0_train, names_rf = names_rf, print_checking = print_checking,
                                       focal_units = focal_units, bound_outcomes = bound_outcomes,
                                       min_treated = min_treated,
                                       Z_int_neighbors = Z_int_neighbors, numcores = numcores, reg = reg,
                                       cross_fitting = cross_fitting, cross_fitting_vector = cross_fitting_vector,
                                       pass_num_neighbors = pass_num_neighbors, verbose = verbose)

      if(k == 1) {results = rec1
      } else {
        results = rbind(results, rec1)
      }
      k = k + 1
      }
    return(results)
  }

}


