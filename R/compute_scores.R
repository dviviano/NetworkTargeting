## The function computes the effects
## input: W_0: adjacency matrix
##        Z: covariate matrix
##        Z_int: matrix which contains interaction terms with treatment assignment
##        m1: fitted model (class that can predict with newx)
##        type : binomial or continuous (type of outcomes)
##        names_rf: names of the columns of the data (applied only for the random forest method)
##        focal_units: which units are focal (1 for focal 0 otherwise)
##
##  The effect are computed without propensity score adjustment
## With a function m1 which is appropriately trained using Z, D*Z_int, num_treated_neighb/num_neighb*Z_int, num_neighbors, num_neighbors*D_i, num_treated_neighb*D_i, D_i
compute_gs <- function(
    W_0, Z,Z_int, m1, type = 'binomial', no_atoms = T, names_rf = NA, focal_units = rep(1, dim(W_0)[1]),
    Z_int_neighbors = Z_int, cross_fitting = F, pass_num_neighbors = NA, cost = 0, ...) {

  my_model = m1
  non_missing_outcomes <- focal_units

  ## Set up accumulators
  final_gs <- NA ## Effect under treatment
  final_g0 <- NA ## Effect under control

  ## Number of neighbors
  neighbors <- apply(W_0, 1, sum)
  neighbors_column = apply(W_0, 2, sum)
  if(all(is.na(pass_num_neighbors))) pass_num_neighbors = neighbors

  ## Run in loop for all observations
  for (i in 1:dim(W_0)[1]) {

    if (cross_fitting) m1 = my_model[[i]]
    non_missing_outcome_i <- non_missing_outcomes[i]
    nn <- neighbors[i]

    ## THis is the denominator in the regression if passed by the user
    nn_passed = pass_num_neighbors[i]
    nn_column = neighbors_column[i]

    ## Here I am eliminating those units without neighbors by using nn > 0 | nn_column > 0
    if (nn > 0 | nn_column > 0) {                                                                 ##################

      ## Compute effect for any possible number of treated neighbors
      res <- sapply(
        c(0:nn),
        function(z) {
          ## model
          #data_XX1 <- as.matrix(c(Z[i,], Z_int[i,],z/nn*Z_int[i,],  neighbors[i], neighbors[i],z/nn, 1), nrow = 1)
          ### Here is where the model enters
          data_XX1 <- as.matrix(c(Z[i,], Z_int[i,],z/max(nn_passed, 1)*Z_int_neighbors[i,],
                                  z/max(nn_passed, 1), z/max(nn_passed, 1)*Z_int_neighbors[i,], z/max(nn_passed, 1), 1), nrow = 1)
          if (type == 'binomial' && all(is.na(names_rf))) {
            to_return = non_missing_outcome_i * predict(m1, newx = t(data_XX1), type = 'response')  - non_missing_outcome_i * cost
            sapply(to_return, function(x) ifelse(is.na(x) & non_missing_outcome_i == 0, 0, x))
          } else if ( all(is.na(names_rf) == F)) {
              data_predict <- as.data.frame(t(data_XX1))
              colnames(data_predict) <- names(data_predict) <- names_rf
              if(sum(is.na(data_predict)) > 0) {
                data_predict[is.na(data_predict)] = 0
                warning('Missing values in covariates. Replaced with 0s, make sure you select focal_unit = F for such individuals.')
              }
              to_return = non_missing_outcome_i *  as.numeric(as.character(predict(m1, newdata = as.data.frame(data_predict), type = 'prob')[2]))  - non_missing_outcome_i * cost ## predict probability of being one
              sapply(to_return, function(x) ifelse(is.na(x) & non_missing_outcome_i == 0, 0, x))
          } else {
            ## Return the effect under treatment with glmnet (linear model)
            to_return = non_missing_outcome_i * as.numeric(coef(m1))%*%c(1, as.numeric(data_XX1)) - non_missing_outcome_i * cost
            sapply(to_return, function(x) ifelse(is.na(x) & non_missing_outcome_i == 0, 0, x))
          }
      })

      ## Return the effect under control
      res2 <-  sapply(c(0:nn), function(z) {
        #data_XX0 <- as.matrix(c(Z[i,], 0*Z_int[i,],z/nn*Z_int[i,], neighbors[i],0*neighbors[i],z/nn, 0), nrow = 1)
        ## Here is where the model enters
        data_XX0 <- as.matrix(c(Z[i,], 0*Z_int[i,],
                                z/max(nn_passed, 1)*Z_int_neighbors[i,], 0,0*Z_int_neighbors[i,], z/max(nn_passed, 1), 0), nrow = 1)

        if (type == 'binomial' && all(is.na(names_rf))) {
          to_return = non_missing_outcome_i * predict(m1, newx = t(data_XX0), type = 'response')
          sapply(to_return, function(x) ifelse(is.na(x) & non_missing_outcome_i == 0, 0, x))
        } else if  (all(is.na(names_rf)== F)){
          data_predict <- as.data.frame(t(data_XX0))
          names(data_predict) <- colnames(data_predict) <- names_rf
          if(sum(is.na(data_predict)) > 0) data_predict[is.na(data_predict)] = 0
          to_return = non_missing_outcome_i *as.numeric(as.character(predict(m1, newdata = as.data.frame(data_predict), type = 'prob')[2])) ## predict probability of being one
          sapply(to_return, function(x) ifelse(is.na(x) & non_missing_outcome_i == 0, 0, x))
        } else {
          to_return = non_missing_outcome_i * as.numeric(coef(m1))%*%c(1, as.numeric(data_XX0))
          sapply(to_return, function(x) ifelse(is.na(x) & non_missing_outcome_i == 0, 0, x))
        }
      })

      final_gs <- c(final_gs , res)
      final_g0 <- c(final_g0, res2)

    }

  }

  ## Remove first entry
  if (no_atoms) {

    final_gs <- final_gs[-1]
    final_g0 <- final_g0[-1]
    ## Return the effects

    return(list(final_gs - final_g0, final_g0, final_gs))

  } else {
    ## Introduce an error for individuals with no neighbors for simplicity
    if(cross_fitting) m1 = my_model[[1]]
    ## Case where only individuals with no neighbors are included
    data_XX1 <- cbind(Z, Z_int,0*Z_int_neighbors,  0, 0*Z_int_neighbors , 0, 1)
    data_XX0 <- cbind(Z, 0*Z_int,0*Z_int_neighbors,  0, 0*Z_int_neighbors, 0, 0)
    if (type == 'binomial' && all(is.na(names_rf))) {

      final_gs <- non_missing_outcomes * predict(m1, newx = as.matrix(data_XX1), type = 'response') - non_missing_outcomes * cost
      final_g0 <- non_missing_outcomes * predict(m1, newx = as.matrix(data_XX0), type = 'response')
      final_gs = apply(cbind(final_gs, non_missing_outcomes), 1,  function(x) ifelse(is.na(x[1]) & x[2] == 0, 0, x[1]))
      final_g0 = apply(cbind(final_g0, non_missing_outcomes), 1,  function(x) ifelse(is.na(x[1]) & x[2] == 0, 0, x[1]))

    } else if ( is.na(names_rf) == F) {

      data_XX1 <- as.data.frame(data_XX1); names(data_XX1) <- names_rf
      data_XX0 <- as.data.frame(data_XX0); names(data_XX0) <- names_rf
      if(sum(is.na(data_XX1)) > 0) {
        data_XX1[is.na(data_XX1)] = 0
        data_XX0[is.na(data_XX0)] = 0
      }
      final_gs <- non_missing_outcomes * as.numeric(as.character(predict(m1, newdata = data_XX1, type = 'prob')[,2]))  - non_missing_outcomes * cost
      final_g0 <- non_missing_outcomes * as.numeric(as.character(predict(m1, newdata = data_XX0, type = 'prob')[,2]))
      final_gs = apply(cbind(final_gs, non_missing_outcomes), 1,  function(x) ifelse(is.na(x[1]) & x[2] == 0, 0, x[1]))
      final_g0 = apply(cbind(final_g0, non_missing_outcomes), 1,  function(x) ifelse(is.na(x[1]) & x[2] == 0, 0, x[1]))

    } else {

      final_gs <- non_missing_outcomes * predict(m1, newx = as.matrix(data_XX1))  - non_missing_outcomes * cost
      final_g0 <- non_missing_outcomes * predict(m1, newx = as.matrix(data_XX0))
      final_gs = apply(cbind(final_gs, non_missing_outcomes), 1,  function(x) ifelse(is.na(x[1]) & x[2] == 0, 0, x[1]))
      final_g0 = apply(cbind(final_g0, non_missing_outcomes), 1,  function(x) ifelse(is.na(x[1]) & x[2] == 0, 0, x[1]))

    }

    return(list(final_gs - final_g0, final_g0, final_gs))
  }
}



##### The following function compute individual level gs for given
#### number of treated neighbors h and individual index i
## The function assumes you pass m1 and the propensity score
## This function is assumed for sparse and surrogate implementation

compute_individual_gs <- function(
    d, h, i, m1, W_0, Z, Z_int,
    Z_int_neighbors = Z_int, Y, D,
    ps = NA, trimming = 0.00,
    model = 'linear_regression',
    family = 'binomial')
{

  neighbors <- apply(W_0, 1, sum)

  if ("list" %in% class(m1)) m1 <- m1[[i]]

  feature_matrix <- as.matrix(c(
    Z[i,], d*Z_int[i,],
    h/max(neighbors[i], 1)*Z_int_neighbors[i,],
    h/max(neighbors[i], 1)*d,
    h/max(neighbors[i], 1)*d*Z_int_neighbors[i,],
    h/max(neighbors[i], 1),
    d), nrow = 1)

  if (!all(is.na(ps))) {
    indicator <- (as.vector(W_0 %*% D)[i] == h) & (D[i] == d)
    ec <- ps[i] * compute_ps_second_term(i, h, ps, W_0)
    ec <- max(ec, trimming)
  }

  if (family == 'binomial' && model != 'random_forest') {

    mc <- predict(m1, newx = t(feature_matrix), type = 'response')
    if (all(is.na(ps))) {
      to_return <- mc
    } else {
      to_return <- (indicator / ec)*(Y[i] - mc) + mc
    }

  } else if (model == 'random_forest') {

    data_predict <- as.data.frame(t(feature_matrix))

    if(sum(is.na(data_predict)) > 0) {
      data_predict[is.na(data_predict)] <- 0
      warning('Missing values in covariates. Replaced with 0s, make sure you select focal_unit = F for such individuals.')
    }

    mc <- as.numeric(as.character(predict(m1, newdata = as.data.frame(data_predict), type = 'prob')[2]))
    if (all(is.na(ps))) {
      to_return <- mc
    } else {
      to_return <- (indicator / ec)*(Y[i] - mc) + mc
    }

  } else {

    mc <- as.numeric(coef(m1))%*%c(1, as.numeric(feature_matrix))
    if (all(is.na(ps))) {
      to_return <- mc
    } else {
      to_return <- (indicator / ec)*(Y[i] - mc) + mc
    }

  }


  return(as.vector(to_return)  )

}



#####################
### Estimate effects via propensity score
######################

## Helpers

## Estimate individual propensity score
## if reg != NA then the model is not estimated
## reg contains the propensity score predictions
estimate_individual_propensity_score <- function(D, Z, penalized = T, reg = NA, cross_fitting = F, cross_fitting_vector = NA){


  if(all(is.na(reg)) == F){
    return(reg)
  }


  which_is_na = apply(Z, 1, function(x) sum(is.na(x)) > 0)
  prop1 = rep(NA, length(D))
  if(cross_fitting == F){
    if(penalized){
      reg <- glmnet::cv.glmnet(y = as.factor(D)[which_is_na == F], x = as.matrix(Z)[which_is_na == F,], family = 'binomial')
      prop1[which_is_na == F] <- predict(reg, newx = as.matrix(Z)[which_is_na == F,], type = 'response')
      if(length(which_is_na) > 0){
        prop1[which_is_na == T] <- mean(predict(reg, newx = as.matrix(Z)[which_is_na == F,], type = 'response'))
      }
    } else{
      reg <- glmnet::glmnet(y = as.factor(D)[which_is_na == F], x = as.matrix(Z)[which_is_na == F,], family = 'binomial', lambda = exp(-7))
      prop1[which_is_na == F] <- predict(reg, newx = as.matrix(Z)[which_is_na == F,], type = 'response')
      if(length(which_is_na) > 0){
      prop1[which_is_na == T] <- mean(predict(reg, newx = as.matrix(Z)[which_is_na == F,], type = 'response'))
      }
      }
    return(prop1)
  } else {

    prop1 = rep(NA, length(D))
    for(j in sort(unique(cross_fitting_vector))){
      if(penalized){
        reg <- glmnet::cv.glmnet(y = as.factor(D)[cross_fitting_vector != j & which_is_na == F], x = as.matrix(Z)[cross_fitting_vector != j & which_is_na == F,], family = 'binomial')

        matrix_predict = as.matrix(Z)[cross_fitting_vector == j & which_is_na == F,]
        if(is.null(dim(matrix_predict))){
        matrix_predict = t(as.matrix(matrix_predict, ncol= 1))
        }
        prop1[cross_fitting_vector == j & which_is_na == F] <- predict(reg, newx = matrix_predict, type = 'response')
        if(length(which_is_na) > 0){
        prop1[cross_fitting_vector == j & which_is_na] <- mean(predict(reg, newx = matrix_predict, type = 'response'))
        }
        } else{
        reg <- glmnet::glmnet(y = as.factor(D)[cross_fitting_vector != j & which_is_na == F], x = as.matrix(Z)[cross_fitting_vector != j & which_is_na == F,], family = 'binomial', lambda = exp(-7))

        matrix_predict = as.matrix(Z)[cross_fitting_vector == j & which_is_na == F,]
        if(is.null(dim(matrix_predict))){
          matrix_predict = t(as.matrix(matrix_predict, ncol= 1))
        }

        prop1[cross_fitting_vector == j & which_is_na == F] <- predict(reg, newx = matrix_predict, type = 'response')
        if(length(which_is_na) > 0){
        prop1[cross_fitting_vector == j & which_is_na] <- mean(predict(reg, newx = matrix_predict, type = 'response'))
        }
        }
    }
    return(prop1)
  }

  }

## prop_Ni: propensity score of neighbors of individual of interest
## s: number of treated neighbors
## l: number of neighbors
## trimming : trimming value

compute_individual_specific_prob <- function(prop_Ni, s, l, trimming = 0.05){
  if(l > 1){
    ## Count the number of combinations:
    gg <- expand.grid(rep(list(0:1), l))
    ## Keep only rows that sum up to s
    gg <- gg[apply(gg, 1, function(x) sum(x) == s),]
    ## Compute products of probabilities
    probabilities <- apply(gg, 1, function(x) prod(prop_Ni[x == 1], 1 - prop_Ni[x == 0]))
    ### Sum such probabilities
    final_prob <- max(trimming, sum(probabilities))
    return(final_prob)} else if(l == 1){
      ifelse(s == 1, max(trimming, prop_Ni), max(trimming, 1 - prop_Ni))
    } else if (l == 0){
      return(1)
    }
}

## Compute the second term for the propensity score
## Note: Here is missing the stratification device (strata from the params list)

compute_ps_second_term <- function(
    i, h, ps, W_0)
{
  prop_Ni <- ps[as.vector(W_0[i,])]
  num_neighb <- sum(W_0[i,])
  if (num_neighb > 1) {
    gg <- expand.grid(rep(list(0:1), num_neighb)) # Count the number of combinations
    gg <- gg[apply(gg, 1, function(x) sum(x) == h),] # Keep only rows that sum up to s
    probabilities <- apply(gg, 1, function(x) prod(prop_Ni[x == 1], 1 - prop_Ni[x == 0])) # Compute products of probabilities
    final_prob <- sum(probabilities) # Sum such probabilities
    return(final_prob)
  } else if(num_neighb == 1){
    ifelse(h == 1, prop_Ni, 1 - prop_Ni)
  } else if (num_neighb == 0){
    return(1)
  }
}


## Compute the effects using the propensity score adjustment
## also_reg_adj: also include regression adjustment
## m1: trained model for the regression adjustment
## type: support of outcome
## penalized: use penalized regression for estimation of propensity score
## num_strata: either NA or number (if number it defines the number of different treatment status depending on number of treated units)
## Assumed: each individual has at least one neighbor
## no_atoms: if T every individual has at least one neighbor
##           if F then all individuals have no neighbors (the two cases are separated for estimation)
## W_0_train: adjacency matrix used for training (if different from the one for targeting)
## trimmed_subpopulation: use the regression adjustment but not the IPW component on the trimmed subpopulation (otherwise it uses the whole population)
## with bounded ipw estimator by trimming value
## units who are focal but with missing outcomes use the non-DR method (imputed outcome)
## reg: regression function (if NA is estimated)
## cross_fitting_vector vector of 0/1 for cross fitting
## if cross-fitting = T, m1 is a list of models, one for each unit

compute_DR_effects <- function(Y,D,W_0, Z, Z_int = NA, also_reg_adj = T, m1 = NA,
                               type = 'binomial', penalized = F, num_strata = NA,
                               trimming = 0.05, no_atoms = T, W_0_train= W_0, names_rf = NA,
                               focal_units = rep(1, length(Y)),  bound_outcomes = c(-Inf, Inf),
                               trimmed_subpopulation = F,
                               IPW = F,
                               Z_int_neighbors = Z_int, reg = NA, cross_fitting = F, cross_fitting_vector = NA,
                               pass_num_neighbors = NA, cost = 0){

  prop_D <- as.vector(estimate_individual_propensity_score(
    D, Z, penalized, reg = reg, cross_fitting = cross_fitting, cross_fitting_vector = cross_fitting_vector))
  no_neighb <- which(apply(W_0, 1, sum) == 0 & apply(W_0, 2, sum) == 0)

  if (no_atoms == F & length(no_neighb) > 1) {
    ## Remove those with neighbors
    Y <- Y[no_neighb]
    D <- D[no_neighb]
    W_0 <- W_0[no_neighb, no_neighb]
    W_0_train <- W_0_train[no_neighb, no_neighb]
    Z <- Z[no_neighb,]
    Z_int <- Z_int[no_neighb,]
    Z_int_neighbors <- Z_int_neighbors[no_neighb,]
    prop_D <- prop_D[no_neighb]
    focal_units  <- focal_units[no_neighb]
  }

  if (no_atoms & length(no_neighb) > 0) {
    ## Remove those without neighbors
    Y <- Y[-no_neighb]
    D <- D[-no_neighb]
    W_0 <- W_0[-no_neighb, -no_neighb]
    W_0_train <- W_0_train[-no_neighb, -no_neighb]
    Z <- Z[-no_neighb,]
    Z_int <- Z_int[-no_neighb,]
    Z_int_neighbors <- Z_int_neighbors[-no_neighb,]
    prop_D <- prop_D[-no_neighb]
    focal_units  <- focal_units[-no_neighb]
  }

  non_missing_Y <- focal_units
  nas_Y <- which(is.na(Y))

  ## kil the ipw component if Y is missing
  if(length(nas_Y) > 0) {
    non_missing_Y[nas_Y] <- 0
    Y[nas_Y] <- 0
  }

  ## Z_int and m1 should not be missing
  gs_condexp <- compute_gs(W_0, Z,Z_int, m1, type = type, no_atoms = no_atoms, names_rf = names_rf, focal_units = focal_units,
                           Z_int_neighbors = Z_int_neighbors, cross_fitting = cross_fitting, pass_num_neighbors = pass_num_neighbors,
                           cost = cost)
  final_g1 <- gs_condexp[[3]] ## Conditional mean for treated
  final_g0 <- gs_condexp[[2]] ## Conditional mean for untreated

  ## Apply the bound on the effects if outcome is a bounded random variable with known upper/lower bounds
  final_g1 <- sapply(final_g1, function(x) min(max(x, bound_outcomes[1]), bound_outcomes[2]))
  final_g0 <- sapply(final_g0, function(x) min(max(x, bound_outcomes[1]), bound_outcomes[2]))

  ## Shut down the regression adjustment for IPW
  if (IPW) {
    final_g1 <- final_g0 <- rep(0, length(final_g1))
  }

  neighbors <- apply(W_0, 1, sum) ## number of neighbors (recall # of neighbors is on the rows)
  if (all(is.na(num_strata))) num_strata <- max(neighbors) + 1 ## if num_strata == NA, each number of treated unit corresponds to a different treatment status

  acc <- 1   ## Set up accumulators
  saved_g1 <- NA ## say propensity score adjustments effects
  saved_g0 <- NA

  for (i in 1:length(Y)) {

    nn <- neighbors[i]
    non_missing_Y_i <- non_missing_Y[i]
    ## Define the strata of treatment depending on number of neighbors
    ## the elements with larger number of neighbors are collapsed together
    which_strat <- sapply(seq(from = 0, to = nn, length = num_strata + 1), floor)

    if (nn > 0) {

      res <-sapply(c(0:nn), function(x) {
        if(nn <= num_strata) {
          to_keep1 <- ifelse(W_0_train[i,]%*%D == x, 1, 0)
        } else {
          ## Define in which strata x belongs to
          save_stratum <- sapply(which_strat, function(y) ifelse(y <= x, 1, 0))
          save_stratum <- max(which(save_stratum == 1))
          ## indicator for the propensity score adjustment
          to_keep1 <- ifelse(W_0_train[i,]%*%D >= which_strat[save_stratum] & W_0_train[i,]%*%D <= which_strat[min(save_stratum + 1, length(which_strat))], 1, 0)
        }
        ## subctract the cost which is subtracted from final_g1 for the regression adjustment part
        to_return = non_missing_Y_i * D[i]*to_keep1*(Y[i] - final_g1[acc + x] - cost)
        ifelse(is.na(to_return) & non_missing_Y_i == 0, 0, to_return)
      })

      ## Repeat for the case under control
      res0 <- sapply(c(0:nn), function(x) {
        if (nn <= num_strata) {
          to_keep1 <- ifelse(W_0_train[i,]%*%D == x, 1, 0)
        } else {
          save_stratum <- sapply(which_strat, function(y) ifelse(y <= x, 1, 0))
          save_stratum <- max(which(save_stratum == 1))
          to_keep1 <- ifelse(W_0_train[i,]%*%D >= which_strat[save_stratum] & W_0_train[i,]%*%D <= which_strat[min(save_stratum + 1, length(which_strat))], 1, 0)
        }
        to_return = non_missing_Y_i * (1-D[i])*to_keep1*(Y[i] - final_g0[acc + x])
        ifelse(is.na(to_return) & non_missing_Y_i == 0, 0, to_return)
      })

      ## Compute propensity
      if (nn <= num_strata) {
        ## You look at your true network W_0, but you fix the number of treated units according to W_0_train
        pp <- compute_individual_specific_prob(prop_D[which(W_0[i,] == 1)], W_0_train[i,]%*%D, nn, trimming = trimming)
      } else {
        ## Compute probability for the stratum
        save_stratum <- sapply(which_strat, function(y) ifelse(y <= W_0_train[i,]%*%D, 1, 0))
        save_stratum <- max(which(save_stratum == 1))
        vec_to_sapply <- c(ceiling(which_strat[save_stratum]):floor(which_strat[min(save_stratum + 1, length(which_strat))]) )
        ## Probability to be in a given threshold , same as before here use W_0
        save_prob <- sapply(vec_to_sapply, function(x) compute_individual_specific_prob(prop_D[which(W_0[i,] == 1)], x, nn, trimming =trimming))
        pp <- sum(save_prob)
      }
      ## Save the propensity score adjustment under treatment and control with trimmed componet
      saved_g1 <- c(saved_g1 , res/(prop_D[i] * pp) * ifelse(trimmed_subpopulation & prop_D[i] * pp <= trimming, 0, 1) )
      saved_g0 <- c(saved_g0 , res0/((1 - prop_D[i]) * pp) * ifelse(trimmed_subpopulation & (1 - prop_D[i]) * pp <= trimming, 0, 1) )
      acc <- nn + 1 + acc ## update the accumulator
    }

    if (nn == 0 & no_atoms) {
      ## Get the people which are connection of somebody else but with no connections themselves
      ## subctract the cost here for the regression adjustment (final_g1 subctracts the cost)
      res <- non_missing_Y_i * D[i]*(Y[i] - final_g1[acc] - cost)
      res = ifelse(is.na(res) & non_missing_Y_i == 0, 0, res)
      res0 <-  non_missing_Y_i * (1-D[i])*(Y[i] - final_g0[acc])
      res0 = ifelse(is.na(res0) & non_missing_Y_i == 0, 0, res0)
      saved_g1 <- c(saved_g1 , res/(prop_D[i]) * ifelse(trimmed_subpopulation & prop_D[i] <= trimming, 0, 1) )
      saved_g0 <- c(saved_g0 , res0/((1 - prop_D[i])) * ifelse(trimmed_subpopulation & (1 - prop_D[i]) <= trimming, 0, 1) )
      acc <- 1 + acc ## update the accumulator
    }

  }

  if (no_atoms) {
    ## remove the intial starting value -- report those completely separated by everything else
    saved_g1 <- saved_g1[-1]
    saved_g0 <- saved_g0[-1]
    ## Return the effects for each given number of treated neighbors
    return(list(saved_g1 + final_g1 - saved_g0 - final_g0, saved_g0 + final_g0, saved_g1 + final_g1))
  } else {
    ## Case where all individuals have no neighbors
    ## subctract the cost (final_g1 also contains the cost)
    saved_g1 <-  non_missing_Y_i  * (Y -  final_g1 - cost)*D/prop_D
    saved_g1 = apply(cbind(saved_g1, non_missing_Y_i), 1, function(x) ifelse(x[2] == 0, 0, x[1]))
    saved_g0 <-  non_missing_Y_i  * (1 - D)*(Y - final_g0)/(1 - prop_D)
    saved_g0 = apply(cbind(saved_g0, non_missing_Y_i), 1, function(x) ifelse(x[2] == 0, 0, x[1]))
    return(list(saved_g1 + final_g1 - saved_g0 - final_g0, saved_g0 + final_g0, saved_g1 + final_g1))
  }

}


