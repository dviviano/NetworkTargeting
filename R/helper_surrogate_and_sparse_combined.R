`%dopar%` <- foreach::`%dopar%`
# Estimate optimal policy -- surrogate and sparse matrix code
#


newm_sparse_or_surrogate <- function(
  Y, Z, D, W_0,
  W_0_train = W_0,
  X = Z, Z_int = Z,
  cost = 0,
  method = "MILP_sparse",
  focal_units = rep(1, length(Y)),
  model = "penalized_regression",
  family = "gaussian",
  maximum_treated = NA,
  maximum_time = 100,
  params = list(), ...
)
{





  if (method == "MILP_sparse") {

    result <- newm_linear_policy_sparse(
      Y = Y, Z = Z, X = X, Z_int = Z_int, D = D, W_0 = W_0, W_0_train = W_0_train,
      cost = cost, eps = params$tolerance_constraint1,
      family = family,
      model_cond_mean = model,
      focal_units = focal_units,
      penalized = params$penalized,
      doubly_robust = params$doubly_robust,
      cross_fitting = params$cross_fitting,
      n_folds = params$n_folds,
      trimming = params$trimming,
      timelimit = maximum_time,
      verbose = params$verbose,
      solver = params$solver,
      prop_score = NA,
      cond_mean_func = params$m1,
      B = params$B, low_b = params$low_b,
      max_treated = maximum_treated,
      numcores = params$numcores
    )
    return(result)

  } else if (method == "surrogate") {

    if (params$numcores > 1) {
      cl <-parallel::makeCluster(params$numcores)
      doParallel::registerDoParallel(cl)
    }
    rseeds <- sample.int(10000, params$numcores, replace = FALSE)
    results <- foreach::foreach(rseed = rseeds, .packages = "torch") %dopar% {
      newm_linear_policy_surrogate_autograd(
        rseed = rseed, Y = Y, Z = Z, X = X, Z_int = Z_int,
        D = D, W_0 = W_0, W_0_train = W_0_train,
        cost = cost, eps = params$tolerance_constraint1,
        family = family,
        model_cond_mean = model,
        focal_units = focal_units,
        penalized = params$penalized,
        doubly_robust = params$doubly_robust,
        cross_fitting = params$cross_fitting,
        n_folds = params$n_folds,
        trimming = params$trimming,
        timelimit = maximum_time,
        verbose = params$verbose,
        solver = params$solver,
        prop_score = NA,
        cond_mean_func = params$m1,
        B = params$B, low_b = params$low_b,
        max_treated = maximum_treated,
        lambda = params$lambda,
        sigma = params$sigma,
        lr = params$lr,
        epochs = params$epochs,
        numcores = params$numcores,
        cross_validation = params$cross_validation,
        min_size_stop_cv = params$min_size_stop_cv
      )
    }
    losses <- sapply(results, \(x) x$result$loss)
    smallest_loss_idx <- which(losses == min(losses))
    result <- results[[smallest_loss_idx]]
    if (params$numcores > 1) parallel::stopCluster(cl)
    return(result)
  } else if (method == "surrogate_graph_clustering") {

    sg_members <- subgraph_membership(
      W_0_train,
      maxtime = params$maxtime_partitions, numcores = params$numcores,
      k = 1, num_cut = floor(sum(focal_units)/params$length_small_samples),
      minimum_size_stopping_condition = params$length_small_samples
    )

    my_models = newm_linear_policy_surrogate_autograd(rseed = rseed, Y = Y, Z = Z, X = X, Z_int = Z_int, D = D, W_0 = W_0, W_0_train = W_0_train,
                                                      cost = cost, eps = params$tolerance_constraint1,
                                                      family = family,
                                                      model_cond_mean = model,
                                                      focal_units = focal_units,
                                                      penalized = params$penalized,
                                                      doubly_robust = params$doubly_robust,
                                                      cross_fitting = params$cross_fitting,
                                                      n_folds = params$n_folds,
                                                      trimming = params$trimming,
                                                      timelimit = maximum_time,
                                                      verbose = params$verbose,
                                                      solver = params$solver,
                                                      prop_score = NA,
                                                      cond_mean_func = params$m1,
                                                      B = params$B, low_b = params$low_b,
                                                      max_treated = maximum_treated,
                                                      lambda = params$lambda,
                                                      sigma = params$sigma,
                                                      lr = params$lr,
                                                      epochs = params$epochs,
                                                      numcores = params$numcores,
                                                      cross_validation = params$cross_validation,
                                                      min_size_stop_cv = params$min_size_stop_cv,
                                                      return_cond_mean = T)
    m1 = my_models[[1]]
    ps = my_models[[2]]
    cross_fitting_vector = my_models[[3]]
    best_welfare <- 0

    for (sg_id in sort(unique(sg_members))) {
      sg_result <- newm_linear_policy_surrogate_autograd(
        Y = Y[sg_members == sg_id], Z = Z[sg_members == sg_id,],
        X = X[sg_members == sg_id,], Z_int = Z_int[sg_members == sg_id,],
        D = D[sg_members == sg_id], W_0 = W_0[sg_members == sg_id, sg_members == sg_id],
        W_0_train = W_0_train[sg_members == sg_id, sg_members == sg_id],
        cost = cost, eps = params$tolerance_constraint1,
        family = family,
        model_cond_mean = model,
        focal_units = focal_units,
        penalized = params$penalized,
        doubly_robust = params$doubly_robust,
        cross_fitting = params$cross_fitting,
        n_folds = params$n_folds,
        trimming = params$trimming,
        timelimit = maximum_time,
        verbose = params$verbose,
        solver = params$solver,
        prop_score = ps[sg_members == sg_id],
        cond_mean_func = m1[sg_members == sg_id],
        B = params$B, low_b = params$low_b,
        max_treated = maximum_treated,
        lambda = params$lambda,
        sigma = params$sigma,
        lr = params$lr,
        epochs = params$epochs,
        numcores = params$numcores
      )
      suppressWarnings(sg_policy_betas <- (sg_result$result |>
                                             ompr::get_solution(beta[i]) |>
                                             dplyr::select(value))$value)
      sg_policy <- apply(X %*% sg_policy_betas, 1, \(x) ifelse(x > 0, 1, 0))


      sg_welfare <- NetworkWelfare(sg_policy, Y, X, Z, Z_int, D, W_0, W_0_train, cost, focal_units,
                                               model, family, params)


      ## Check compute welfare function above from here
      if (sg_welfare > best_welfare) {
        best_welfare <- sg_welfare
        best_sg_policy_betas <- sg_policy_betas
        best_sg_policy <- sg_policy
      } else if (sg_id == 1) {
        best_welfare <- sg_welfare
        best_sg_policy_betas <- sg_policy_betas
        best_sg_policy <- sg_policy
      }
    }
    return(list("result" = list("policy" = best_sg_policy,
                                "betas" = best_sg_policy_betas,
                                "welfare" = best_welfare),
                "cond_mean_func" = m1, "prop_score" = ps))
  }

}
