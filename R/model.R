
## Compute conditional mean function

compute_conditional_mean_function = function(Y, Z, Z_int, D, Z_int_neighbors, W_0_train, model, family, cross_fitting = F, cross_fitting_vector = NA,
                                             focal_units = rep(1, length(Y)), share_treated_neighbors = NA,
                                             pass_num_neighbors = NA){

  names_rf = NA
  num_neighb <- apply(W_0_train,1,sum)
  num_neighb <- sapply(num_neighb, function(x) max(x,1)) ## Threshold at 1
  if(all(is.na(pass_num_neighbors))) {
    pass_num_neighbors = num_neighb
  } else {
    pass_num_neighbors = sapply(pass_num_neighbors, function(x) max(x, 1))
  }
  num_treated_neighb <- as.vector(W_0_train%*%D)/pass_num_neighbors
  if(all(is.na(share_treated_neighbors)) == F)  num_treated_neighb = share_treated_neighbors
  ## prepare the model
  matrix_to_pass <- cbind(Z, D*Z_int, num_treated_neighb*Z_int_neighbors, num_treated_neighb*D, num_treated_neighb*D*Z_int_neighbors, num_treated_neighb, D)


  Y_pass <- Y
  Y_pass = Y_pass[focal_units == 1]
  matrix_to_pass = matrix_to_pass[focal_units == 1,]
  if(cross_fitting == F){
    if(model == 'penalized_regression') m1 <- glmnet::cv.glmnet(y = Y_pass, x = as.matrix(matrix_to_pass), family = family)
    if(model == 'linear_regression') m1 <- glmnet::glmnet(y = Y_pass, x = as.matrix(matrix_to_pass), lambda = exp(-12), family = family)
    if(model == 'random_forest' & family == 'binomial') {
      to_pass_rf <- as.data.frame(matrix_to_pass)
      names_rf <- names(to_pass_rf)
      m1 <- randomForest::randomForest(y = as.factor(Y_pass), x = to_pass_rf, lambda = exp(-6))
    }
  } else {

    # my_cut = opposite_min_cut_MILP_1_and_2_degree(W_0_train)
    m1_accumulator = list()
    my_cut2 = cross_fitting_vector$cut2
    acc_num = 1
    for(j in sort(unique(my_cut2))){

      if(model == 'penalized_regression') m1_accumulator[[acc_num]] <- glmnet::cv.glmnet(y = Y_pass[my_cut2[focal_units == 1] != j], x = as.matrix(matrix_to_pass)[my_cut2[focal_units == 1] != j,], family = family)
      if(model == 'linear_regression') m1_accumulator[[acc_num]] <- glmnet::glmnet(y = Y_pass[my_cut2[focal_units == 1] != j], x = as.matrix(matrix_to_pass)[my_cut2[focal_units == 1] != j,], lambda = exp(-12), family = family)
      if(model == 'random_forest' & family == 'binomial') {
        to_pass_rf <- as.data.frame(matrix_to_pass[my_cut2[focal_units == 1] != j,])
        names_rf <- names(to_pass_rf)
        m1_accumulator[[acc_num]] <- randomForest::randomForest(y = as.factor(Y_pass)[my_cut2[focal_units == 1] != j], x = to_pass_rf, lambda = exp(-6))
      }
      acc_num = acc_num + 1
    }
    m1 = list()
    for(i in 1:dim(W_0_train)[1]){
      ## Switch the conditional mean functions that you use for estimation
      m1[[i]] = m1_accumulator[[which(sort(unique(my_cut2)) == my_cut2[i])]]
    }
  }
  return(list(m1, names_rf, matrix_to_pass))

}
