newm_linear_policy_surrogate_autograd <- function(
    rseed,
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
    lambda = 1,
    sigma = 1,
    lr = 1e-2,
    epochs = 100,
    numcores = 2,
    return_cond_mean = F,
    cross_validation = list(),
    min_size_stop_cv = 20)
{

  stopifnot(solver %in% c('glpk', 'gurobi'))

  ## compute conditional mean function, used to impute score

  cross_fitting_vector = NA

  ## Set the connections between non-focal points equal to zero
  ## (these are not used for estimation purposes and they reduce memory requirements)

  if(sum(focal_units) < length(Y)) {
    W_0[focal_units == 0, focal_units == 0] <- 0
    W_0_train[focal_units == 0, focal_units == 0] <- 0
  }

  p <- dim(X)[2]


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

  ## optionally compute propensity score for doubly robust estimation

  if (doubly_robust & all(is.na(prop_score))) {
    ps <- estimate_individual_propensity_score(
      Z=Z,  D=D,
      cross_fitting = cross_fitting, cross_fitting_vector = cross_fitting_vector,
      penalized = penalized)
  } else if (doubly_robust) {
    ps <- prop_score
  } else {
    ps <- NA
  }

  if(return_cond_mean) return(list(m1 = m1, ps = ps, cross_fitting_vector = cross_fitting_vector))


  ## assemble other objects needed to set up optimization problem

  A <- W_0
  neighb <- apply(A, 1, sum)
  n <- length(Y)
  X <- cbind(1, X) ## Add intercept
  Xt <- torch::torch_tensor(X)
  At <- torch::torch_tensor(A)
  max_h <- max(neighb)
  max_ht <- torch::torch_arange(0, max_h)

  q1 <- matrix(NA, nrow = n, ncol = max_h+1)
  q0 <- matrix(NA, nrow = n, ncol = max_h+1)

  for (i in 1:n) {
    for (h in 0:max_h) {
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

  q1t <- torch::torch_tensor(q1)
  q0t <- torch::torch_tensor(q0)

  ## set up optimization problem
  # local_torch_manual_seed(99)

  policy_func <- torch::nn_module(
    clasname = "policy_func",
    initialize = function(in_features, out_features, sigma) {
      self$w <- torch::nn_parameter(torch::torch_randn(in_features, out_features))
      self$sigma <- sigma
    },
    forward = function(x) {
      torch::distr_normal(loc=0, scale=self$sigma)$cdf(torch::torch_mm(x, self$w))
    }
  )

  loss_fn <- function(pred, At, q1t, q0t, neighb, lambda=lambda) {

    n <- At$shape[1]

    neighb_sums <- torch::torch_mm(At, pred)
    h_subs <- max_ht$`repeat`(c(n,1))
    distt <- (neighb_sums - h_subs)$pow(2)
    B <- torch::torch_exp(-lambda*distt)

    losst <- B * (torch::torch_multiply((q1t - cost - q0t), pred) + q0t)
    lostvec <- losst$reshape(-1)

    neighbseq <- sapply(neighb, \(x) 0:x)
    ihpairs <- dplyr::bind_rows(mapply(\(i, h) expand.grid(i, h), i=1:n, h=neighbseq, SIMPLIFY=F))
    ihind <- mapply(\(i, h) ((i-1)*(max_h+1)) + (h+1), i=ihpairs[,1], h=ihpairs[,2])

    f_loss <- -1*lostvec[ihind]$sum()

    return(f_loss)

  }

  ### run optimization

  cross_validation_defaults <- list('sigma'=c(sigma),'lambda'=c(lambda),'lr'=c(lr), 'epochs'=c(epochs))
  cv_list <- modifyList(cross_validation_defaults, cross_validation)
  cv_grid <- as.matrix(expand.grid(cv_list))
  best_hp <- NA

  if (nrow(cv_grid) > 1) {
    sg_members <- subgraph_membership(
      W_0_train, maxtime = timelimit, numcores = numcores,
      k = 1, num_cut = floor(sum(focal_units)/min_size_stop_cv),
      minimum_size_stopping_condition = min_size_stop_cv,
      focal_units = focal_units
    )
    num_folds <- length(unique(sg_members))
    cv_losses <- rep(NA, nrow(cv_grid))
    print(paste("Cross-validating on", num_folds, "subgraph-based folds..."))
    for (i in 1:nrow(cv_grid)) {
      curr_grid <- cv_grid[i,]
      print(list(curr_grid))
      cv_inner_loss <- rep(NA, num_folds)
      for (K in 1:num_folds) {
        train_units <- which(sg_members!=K)
        test_units <- which(sg_members==K)
        model <- policy_func(dim(X)[2], 1, curr_grid[['sigma']])
        optimizer <- torch::optim_sgd(model$parameters, lr = curr_grid[['lr']])
        Xt_train <- torch::torch_tensor(X[train_units,])
        Xt_test <- torch::torch_tensor(X[test_units,])
        At_train <- torch::torch_tensor(A[train_units,train_units])
        At_test <- torch::torch_tensor(A[test_units,test_units])
        q1t_train <- torch::torch_tensor(q1[train_units,])
        q0t_train <- torch::torch_tensor(q0[train_units,])
        q1t_test <- torch::torch_tensor(q1[test_units,])
        q0t_test <- torch::torch_tensor(q0[test_units,])
        neighb_train <- apply(A[train_units,train_units], 1, sum)
        neighb_test <- apply(A[test_units,test_units], 1, sum)
        for (epoch in 1:curr_grid[['epochs']]) {
          optimizer$zero_grad()
          pred <- model(Xt_train)
          loss <- loss_fn(
            pred, At_train, q1t_train, q0t_train,
            neighb_train, lambda=curr_grid[['lambda']])
          loss$backward()
          optimizer$step()
        }
        optimizer$zero_grad()
        pred_test <- model(Xt_test)
        cv_inner_loss[K] <- loss_fn(
          pred_test, At_test, q1t_test, q0t_test,
          neighb_test, lambda=curr_grid[['lambda']])$item()
      }
      cv_losses[i] <- mean(cv_inner_loss)
    }
    best_hp <- cv_grid[which(cv_losses==min(cv_losses)),]
  }

  if (!all(is.na(best_hp))) {
    print("Best CV hyperparameters:")
    print(best_hp)
    sigma <- best_hp[['sigma']]
    lambda <- best_hp[['lambda']]
    lr <- best_hp[['lr']]
    epochs <- best_hp[['epochs']]
  }

  model <- policy_func(dim(X)[2], 1, sigma)
  optimizer <- torch::optim_sgd(model$parameters, lr = lr)
  for (epoch in 1:epochs) {
    optimizer$zero_grad()
    pred <- model(Xt)
    loss <- loss_fn(pred, At, q1t, q0t, neighb, lambda=lambda)
    loss$backward()
    optimizer$step()
  }
  final_loss <- loss$item()
  beta_star <- torch::as_array(model$w)
  assn <- as.numeric(X %*% beta_star >= 0)
  probablistic_policies = sapply(X %*% beta_star, function(x) pnorm(x/sigma))
  newm_result <- list("beta" = beta_star, "assn" = assn, "loss" = final_loss,
                      'probablistic_policies' =  probablistic_policies)

  return(list("result" = newm_result, "cond_mean_func" = m1, "prop_score" = ps))

}
