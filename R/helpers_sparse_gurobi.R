## This function follows similarly to optimal_network_linear_rule function (under 'main_function.R')
## But uses a sparse matrix of constraints and GUROBI optimization, allowing for better RAM usage
## but slower computational time

newm_linear_policy_sparse <- function(
    Y, X, Z, Z_int, D, W_0, W_0_train,
    cost = 0, eps = 1e-4,
    family = 'gaussian',
    model_cond_mean = "linear_regression",
    ## penalized indicates whether to penalize for prop score estimation
    penalized = T,
    doubly_robust=F,
    cross_fitting = T,
    n_folds = 5,
    trimming = 0.01,
    timelimit=100,
    verbose=F,
    solver="gurobi",
    prop_score = NA,
    cond_mean_func = NA,
    Z_int_neighbors = Z_int,
    B = 1,
    low_b = -1,
    max_treated = NA,
    focal_units = rep(1, length(Y)),
    numcores = 2)
{

  cross_fitting_vector = NA

  missing_X = apply(X, 1, function(x) sum(is.na(x)) > 0)
  if (sum(missing_X) > 0) focal_units[missing_X] = 0

  missing_Z = apply(Z, 1, function(x) sum(is.na(x)) > 0)
  if(sum(missing_Z) > 0) focal_units[missing_Z] = 0

  missing_Y = sapply(Y, is.na)
  if(sum(missing_Y) > 0) focal_units[missing_Y] = 0

  ## Set the connections between non-focal points equal to zero
  ## (these are not used for estimation purposes and they reduce memory requirements)

  if(sum(focal_units) < length(Y)) {
    W_0[focal_units == 0, focal_units == 0] <- 0
    W_0_train[focal_units == 0, focal_units == 0] <- 0
  }

  if (is.na(cond_mean_func)[1]) {

    if (cross_fitting == F) {

      my_model = compute_conditional_mean_function(Y, Z, Z_int, D, Z_int_neighbors, W_0_train, model_cond_mean, family, cross_fitting = F, cross_fitting_vector = NA,
                                                   focal_units = focal_units)
      m1 = my_model[[1]]
      names_rf = my_model[[2]]
      cross_fitting_vector = NA


    } else {
      my_cut = opposite_min_cut_MILP_1_and_2_degree(W = W_0_train, n_folds, slack = 5, focal_units = focal_units)
      my_model = compute_conditional_mean_function(Y, Z, Z_int, D, Z_int_neighbors, W_0_train, model_cond_mean, family, cross_fitting = T, cross_fitting_vector = my_cut,
                                                   focal_units = focal_units)
      m1 = my_model[[1]]
      names_rf = my_model[[2]]
      ## Cross fitting vector for the propensity score to be passed to the functions below
      cross_fitting_vector = my_cut$cut1

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
      my_cut = opposite_min_cut_MILP_1_and_2_degree(W = W_0_train, n_folds, slack = 5, focal_units = focal_units)
      ## Cross fitting vector for the propensity score to be passed to the functions below
      cross_fitting_vector = my_cut$cut1
      m1 = cond_mean_func
      if(sum(missing_X) > 0) {
        cross_fitting_vector = cross_fitting_vector[-missing_X]
        m1 = cond_mean_func[-missing_X]
      }
    }
  }



  ## Remove missing values for the subsequent estimation
  if(sum(missing_X) > 0){

    Y = Y[-missing_X]
    X = X[-missing_X, ]
    Z = Z[-missing_X, ]
    Z_int = Z_int[-missing_X, ]
    D = D[-missing_X]
    W_0 = W_0[-missing_X, -missing_X]
    focal_units = focal_units[-missing_X]
    Z_int_neighbors = Z_int_neighbors[-missing_X,]
    pass_num_neighbors = pass_num_neighbors[-missing_X]

  }


  ## optionally compute propensity score for doubly robust estimation

  if (doubly_robust & all(is.na(prop_score))) {
    ps <- estimate_individual_propensity_score(
      Z=Z, D=D,
      cross_fitting = cross_fitting, cross_fitting_vector = cross_fitting_vector,
      penalized = penalized)
  } else if (doubly_robust) {
    ps <- prop_score
  } else {
    ps <- NA
  }

  ## assemble other objects needed to set up optimization problem

  X <- cbind(1, X) ## add intercept
  A <- W_0
  neighb <- apply(A, 1, sum)
  n <- length(Y)
  max_h <- max(neighb)

  Xfunc <- function(i,t){
    X_sub <- X[i, t]
    X_vec <- as.vector(X_sub)
    return(X_vec)
  }

  Afunc <- function(i,t){
    A_sub <- A[i, t]
    A_vec <- as.vector(A_sub)
    return(A_vec)
  }

  q1 <- matrix(NA, nrow = n, ncol = max_h+1)
  q0 <- matrix(NA, nrow = n, ncol = max_h+1)

  for (i in 1:n) {
    for (h in 0:neighb[i]) {
      q1[i, h+1] <- compute_individual_gs(
        d=1, h=h, i=i, m1=m1, W_0=W_0_train, Z=Z, Z_int=Z_int, Y=Y, D=D, ps=ps,
        trimming=trimming, family=family, model=model_cond_mean)
      q0[i, h+1] <- compute_individual_gs(
        d=0, h=h, i=i, m1=m1, W_0=W_0_train, Z=Z, Z_int=Z_int, Y=Y, D=D, ps=ps,
        trimming=trimming, family=family, model=model_cond_mean)
    }
  }

  if(sum(focal_units) > 0){
    q1[focal_units == 0,] = cost ## Use cost since this will cancel out from - cost below
    q0[focal_units == 0,] = 0
  }

  Q1func <- function(i,h){
    q1_sub <- q1[i, h+1]
    q1_vec <- as.vector(q1_sub)
    return(q1_vec)
  }

  Q0func <- function(i,h){
    q0_sub <- q0[i, h+1]
    q0_vec <- as.vector(q0_sub)
    return(q0_vec)
  }

  ## set up optimization problem

  newm_model <- ompr::MILPModel() |>
    # beta
    ompr::add_variable(beta[j,l], type = "continuous", lb = low_b, ub = B, j = 1:dim(X)[2], l=0:0) |>
    # y
    ompr::add_variable(y[i,l], i = 1:n, l=0:0, type = "binary") |>
    ompr::add_constraint(ompr::sum_over(ompr::colwise(Xfunc(i,j))*beta[j,l], j=1:dim(X)[2])/(ompr::sum_over(X[i,j], j=1:dim(X)[2]) + eps) + eps <= y[i], j=1:dim(X)[2], i=1:n, l=0:0) |>
    ompr::add_constraint(ompr::sum_over(ompr::colwise(Xfunc(i,j))*beta[j,l], j=1:dim(X)[2])/(ompr::sum_over(X[i,j], j=1:dim(X)[2]) + eps) + 1   >= y[i], j=1:dim(X)[2], i=1:n, l=0:0) |>
    # x1
    ompr::add_variable(x1[i,h], i = 1:n, h = 0:max_h, type = "binary", h <= neighb[i]) |>
    ompr::add_constraint((ompr::sum_over(ompr::colwise(Afunc(i,k))*y[k,l], k=1:n) - h) / (neighb[i] + 1) + eps <= x1[i,h], i = 1:n, h = 0:max_h, l=0:0, h <= neighb[i]) |>
    ompr::add_constraint((ompr::sum_over(ompr::colwise(Afunc(i,k))*y[k,l], k=1:n) - h) / (neighb[i] + 1) + 1   >= x1[i,h], i = 1:n, h = 0:max_h, l=0:0, h <= neighb[i]) |>
    # x2
    ompr::add_variable(x2[i,h], i = 1:n, h = 0:max_h, type = "binary", h <= neighb[i]) |>
    ompr::add_constraint((h - ompr::sum_over(ompr::colwise(Afunc(i,k))*y[k,l], k=1:n)) / (neighb[i] + 1) + eps <= x2[i,h], i = 1:n, h = 0:max_h, l=0:0, h <= neighb[i]) |>
    ompr::add_constraint((h - ompr::sum_over(ompr::colwise(Afunc(i,k))*y[k,l], k=1:n)) / (neighb[i] + 1) + 1   >= x2[i,h], i = 1:n, h = 0:max_h, l=0:0, h <= neighb[i]) |>
    # z
    ompr::add_variable(z[i,h],  i = 1:n, h = 0:max_h, type = "binary", h <= neighb[i]) |>
    ompr::add_constraint((y[i,l]+x1[i,h]+x2[i,h])/3 - 1 + eps <= z[i,h], i = 1:n, l=0:0, h = 0:max_h, h <= neighb[i]) |>
    ompr::add_constraint((y[i,l]+x1[i,h]+x2[i,h])/3           >= z[i,h], i = 1:n, l=0:0, h = 0:max_h, h <= neighb[i]) |>
    # objective
    ompr::set_objective(ompr::sum_over(
      ((ompr::colwise(Q1func(i,h)) - cost) - ompr::colwise(Q0func(i,h)))*z[i,h] + ompr::colwise(Q0func(i,h))*(x1[i,h]+x2[i,h]-1),
      i = 1:n, h = 0:max_h, h <= neighb[i]), "max")

  if (!is.na(max_treated)) {
    newm_model <- newm_model |> ompr::add_constraint(ompr::sum_over(y[i,l]) <= max_treated, i=1:n, l=0:0)
  }

  ## optimize!

  milp_mat <- newm_model |> ompr::extract_constraints()

  newm_result <- newm_model |>
    ompr::solve_model(
      ompr.roi::with_ROI(
        solver = solver,
        OutputFlag = verbose,
        TimeLimit = timelimit,
        IntFeasTol = 1e-9,
        FeasibilityTol = 1e-9,
        Threads = numcores,
        Presolve = 0,
        Heuristics=0,
        Cuts=0))

  ## result

  assn <- newm_result |> ompr::get_solution(y[i,l]) |> dplyr::select(value) |> unlist()
  beta_star <- newm_result |> ompr::get_solution(beta[i,j]) |> dplyr::select(value) |> unlist()
  newm_res <- list("beta" = beta_star, "assn" = assn)

  return(list("result" = newm_res, "cond_mean_func" = m1, "prop_score" = ps, "milp_mat" = milp_mat))

}

