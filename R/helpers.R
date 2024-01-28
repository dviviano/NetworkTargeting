

## Model checking function 
## The function check if the constraints are satisfied by the final model
##  Checks are performed only on units with neighbors for simplicity 

model_checking <- function(result, XX, Z, W_0, y_s_position, positions, beta_hat, final_gs, final_g0, cc, to_print = T, 
                           focal_units = rep(1, dim(XX)[1]), final_gs_noneighb = 0, final_g0_noneighb = 0){
  
  neighbors <- apply(W_0,1, sum)
  neighbors_column = apply(W_0, 2, sum)
  no_neighb <- which(neighbors == 0 & neighbors_column == 0)
  if(length(no_neighb) > 0){
    W_0 <- W_0[-no_neighb, -no_neighb]
    Z <- Z[-no_neighb,]
    neighbors <- neighbors[-no_neighb]
    neighbors_column = neighbors_column[-no_neighb]
  }
  p <- dim(XX)[2] - 1
  n <- length(y_s_position)
  ## remove the ones with no neighbors for model checking
  
  xx_result = result$x 
  xx_result_noneighb= 0
  if(length(no_neighb) > 1) { 
    xx_result_noneighb =  xx_result[c(1:length(no_neighb))] 
    xx_result <- xx_result[-c(1:length(no_neighb))]
   
    }
  if(to_print){
  print('obj program (excluding those with no neighbors)')
  print(xx_result%*%cc/n - sum(final_g0)/n + sum(final_g0_noneighb)/n + sum(xx_result_noneighb * final_gs_noneighb )/n)
  }
  #  xx_result <- result$x[-c(1:length(no_neighb))]
  predictions <- xx_result[y_s_position]
  z <- xx_result[(2*length(final_gs) + n + 1):(3*length(final_gs) + n)]
  if(to_print){
  print('diff in policies:')}
  
  r1 <- mean(sapply(XX%*%beta_hat, function(hh) ifelse(hh >= 0, 1, 0))!=predictions)
  
  if(to_print){
  print(mean(sapply(XX%*%beta_hat, function(hh) ifelse(hh >= 0, 1, 0))!=predictions))}
  diff <- rep(NA, n)
  for(i in 1:n){
    diff[i] <-  sum(z[positions == i]) == predictions[i]
  }
  r2 <- 1 -mean(diff)
  if(to_print){
  print(1 -mean(diff))}
  
  
  xs <- xx_result[-y_s_position]
  xs <- xs[1:(2*length(final_gs))]
  u <- 1
  zz <- rep(NA, length(z))
  acc <- 0
  for(i in 1:n){ 
    xi <- xs[(acc + 1):(acc + 2*neighbors[i] + 2) ]
    k <- 1
    acc <- acc + 2*neighbors[i] + 2
    for(j in 1:(neighbors[i] + 1)){
      zz[u]  <- prod(xi[k:(k+1)])*predictions[i]
      k <- k + 2 
      u <- u +1}
  }
  if(to_print){
  print('discrepancies in auxiliary zs variables and ys:')}
  r3 <- 1 - mean(z == zz  )
  if(to_print){
  print(1 - mean(z == zz  ))}
  if(to_print){
  print('Percentage of treated')
  print(mean(predictions))}
  ## Consider the xs
  treat_neigh <- W_0%*%predictions
  k <- 1
  save <- 0
  for(i in 1:n){
    xsi <- xs[k:(k + 2*neighbors[i] + 1)]
    xsi1 <- xsi[seq(from = 1, to = length(xsi), by = 2)]
    xsi2 <- xsi[seq(from = 2, to = length(xsi), by = 2)]
    xsii <- cbind(xsi1, xsi2)
    xsin <- apply(xsii, 1, sum)
    pp <- rep(1, length(xsi1))
    pp[treat_neigh[i] + 1] <- 2
    save <- sum(pp != xsin) + save
    k <- 2*neighbors[i] + 1 + k + 1
  }
  if(to_print){
  print('discrepancies between auxiliary variables (xs) and treatment assignments (ys)')}
  r4 <-mean(save)/n
  if(to_print){
  print(mean(save)/n)
  print('beta_hat')
  print(beta_hat)
  }
  zz_t <- xx_t1 <-xx_t2   <- NA
  for(i in 1:n){
    nn_save <- neighbors[i]
    zz <- rep(0, nn_save + 1)
    x1 <- x2 <- rep(1, nn_save + 1)
    treat_ne <- W_0[i,]%*%predictions + 1
    zz[treat_ne] <- 1*predictions[i]
    #x1[treat_ne] <- 1
    x2[-treat_ne] <- 0
    zz_t <- c(zz_t, zz)
    xx_t1 <- c(xx_t1, x1)
    xx_t2 <- c(xx_t2, x2)
  }
  zz_t <- zz_t[-1]
  xx_t1 <- xx_t1[-1]
  xx_t2 <- xx_t2[-1]
  if(to_print){
  print('objective function with auxialiary variables')}
  r5 <- (final_gs%*%z + rep(final_g0, each = 2)%*%xs - sum(final_g0))/n
  if(to_print){
  print((final_gs%*%z + rep(final_g0, each = 2)%*%xs - sum(final_g0))/n)
  print('objective function with actual treatment assignment (should match previous result)')
  }
  r6 <- (final_gs%*%zz_t +final_g0%*%xx_t1 + final_g0%*%xx_t2 - sum(final_g0)  )/n
  if(to_print){
  print(r6)}
  
  #print((final_gs%*%zz_t +final_g0%*%xx_t1 + final_g0%*%xx_t2 - sum(final_g0)  )/n)
  
  
  
  
  
  if(to_print){
  print('average number of binomial variables actually in 0-1')}
  r7 <- mean(result$x[1:(length(result$x) - p - 1)] %in% c(0,1))
  if(to_print){
  print(mean(result$x[1:(length(result$x) - p - 1)] %in% c(0,1)))}
  return(c(r1,r2,r3,r4,r5,r6,r7))
  
}



## Define function to compute the welfare on out of sample units 
## The function takes the predicted policies on such units, the model m1, and whether the method is the doubly robust method
## and it returns the welfare 
## names_rf: colnames of the dataframe used for training if a random forest is used
## policies: list of policies 
## bound_outcomes: lower and upper bounds on the outcome 
## reg: regression function for prop score (if NA is estmated)

compute_welfare <- function(my_policies, Y,D,W_0, 
                            Z, Z_int = NA, m1 = NA, 
                            family = 'binomial', 
                            penalized = F, num_strata = NA, 
                            trimming = 0.01, doubly_robust = F, 
                            names_rf, W_0_train = W_0, 
                            focal_units = rep(1, length(my_policies[[1]])),  
                            bound_outcomes = c(-Inf,Inf), trimmed_subpopulation = T, 
                            IPW = F, 
                            Z_int_neighbors = Z_int, reg  =NA, 
                            pass_num_neighbors = NA, 
                            cross_fitting = F, 
                            cost = 0, cross_fitting_vector = NA, 
                            n_folds = 10, 
                            parallel = F
                           ){
  
  n <- sum(focal_units)
  
  if(all(is.na(cross_fitting_vector)) & doubly_robust & cross_fitting){
    my_cut = opposite_min_cut_MILP_1_and_2_degree(W = W_0, n_folds, slack = 5, focal_units = focal_units)
    cross_fitting_vector = my_cut$cut1
  }
  ## Use doubly robust 

    neighbors <- apply(W_0, 1, sum)
    neighbors_column = apply(W_0, 2, sum)
    no_neighbor <- which(neighbors == 0 & neighbors_column == 0)
    if(all(is.na(pass_num_neighbors)) )pass_num_neighbors = neighbors 
    if(length(no_neighbor) > 0) {
      neighbors <- neighbors[-no_neighbor]
      neighbors_column = neighbors_column[-no_neighbor]
    } 
    ## Obtained the effect on those with at least one neighbor 
   if(doubly_robust) {  gs <- compute_DR_effects(Y,D,W_0, Z, Z_int, also_reg_adj = T, 
                             m1, family, penalized, num_strata, trimming, 
                             names_rf = names_rf, W_0_train = W_0_train, 
                             focal_units = focal_units, trimmed_subpopulation = trimmed_subpopulation, 
                             IPW = IPW, Z_int_neighbors = Z_int_neighbors, reg = reg, 
                             pass_num_neighbors = pass_num_neighbors, cross_fitting = cross_fitting, 
                             cross_fitting_vector = cross_fitting_vector)
   } else { 
     gs <- compute_gs(W_0, Z, Z_int, 
                              m1, family, focal_units = focal_units,
                              pass_num_neighbors = pass_num_neighbors, cross_fitting = cross_fitting, 
                      Z_int_neighbors = Z_int_neighbors, names_rf = names_rf)
     }
    gs1 <- gs[[3]] - cost ## Effect under treatment
    gs0 <- gs[[2]] ## Effect under control
    
    
    ## We obtained the effects gd for all those individuals with no neighbors 
    ## compute results for those with no neighbors 
    
    if(length(no_neighbor) > 1){
      ## It returns the DR scores: these are organized as effects_i_0_treated_neighbors, effect_i_1_treate_neighbor,... with size sum N_i + 1
   if(doubly_robust)  {  gs_noneighb <- compute_DR_effects(Y,D,W_0, Z, Z_int, 
                                        also_reg_adj = T, m1, 
                                        family, penalized, num_strata, 
                                        trimming, no_atoms = F, names_rf=names_rf, 
                                        W_0_train= W_0_train, focal_units = focal_units, 
                                        trimmed_subpopulation = trimmed_subpopulation, 
                                        IPW = IPW, 
                                        Z_int_neighbors = Z_int_neighbors, reg = reg, 
                                        pass_num_neighbors = pass_num_neighbors, 
                                        cross_fitting = cross_fitting, 
                                        cross_fitting_vector = cross_fitting_vector)
   
   } else {
     
     if(cross_fitting) m1 = m1[no_neighbor]
     gs_noneighb <-  compute_gs(W_0[no_neighbor, no_neighbor], Z[no_neighbor,], Z_int[no_neighbor,], 
                                      m1, family, focal_units = focal_units[no_neighbor],
                                      pass_num_neighbors = pass_num_neighbors[no_neighbor], cross_fitting = cross_fitting, no_atoms = F, 
                                Z_int_neighbors = Z_int_neighbors[no_neighbor,])
   }
      gs1_noneighb <- gs_noneighb[[3]] - cost
      gs0_noneighb <- gs_noneighb[[2]]
      
    }
    
    my_effect_final <- rep(NA, length(my_policies))
    if(length(no_neighbor) > 0) focal_neighb = focal_units[-no_neighbor]
    if(parallel == F){ 
    for(j in 1:length(my_policies)){
      
      acc <- 0 ## Set up the accumulator 
      policies_noneighb <- my_policies[[j]]
      policies <- my_policies[[j]]
      if(length(no_neighbor) > 0) policies_noneighb <- policies[-no_neighbor] ## Look at policies of those without any neighbor
      effect <- 0
      if(length(no_neighbor) == 0) num_treated_neighbors <- W_0%*%policies_noneighb
      
      if(length(no_neighbor) > 0) num_treated_neighbors <- W_0[-no_neighbor, -no_neighbor]%*%policies_noneighb
      
      for(i in 1:length(neighbors)){
        num_treated_neighb <- num_treated_neighbors[i]
        effect1 <- gs1[num_treated_neighb + 1 + acc]
        effect0 <- gs0[num_treated_neighb + 1 + acc]
        effect1 = ifelse(is.na(effect1), 0, effect1)
        effect0 = ifelse(is.na(effect0), 0, effect0)
        saved_effect1 =  focal_neighb[i] * (policies_noneighb[i]*effect1 + (1 - policies_noneighb[i])*effect0)
        effect <- effect + ifelse(is.na(saved_effect1),0, saved_effect1) ## look at the effect consistent with the direct effect of the policy
        acc <- acc + neighbors[i] + 1  ## Skip to next observation (after num_neighbor + 1 step you have the following unit)
        ## Store the effect (used for CI)
        if(i == 1){
          saved_effect <- focal_neighb[i] * (policies_noneighb[i]*effect1 + (1 - policies_noneighb[i])*effect0)
        } else{
        saved_effect <- c(saved_effect,(policies_noneighb[i]*effect1 + (1 - policies_noneighb[i])*effect0))
        }
      }
      if(length(no_neighbor) > 1){
        k <- 1
        for(i in no_neighbor){
          my_eff = focal_units[i] * (policies[i]*gs1_noneighb[k] + (1 - policies[i])*gs0_noneighb[k])
          effect <- effect + ifelse(is.na(my_eff),0, my_eff) 
          k <- k  + 1
        }
      }
      my_effect_final[j] <- effect/n
      
    }
    return(list(predictions = saved_effect, welfare =  my_effect_final))
    } else {
      my_effect_final =   foreach(j = 1:length(my_policies), .combine = c)%dopar%{
        
        acc <- 0 ## Set up the accumulator 
        policies_noneighb <- my_policies[[j]]
        policies <- my_policies[[j]]
        if(length(no_neighbor) > 0) policies_noneighb <- policies[-no_neighbor] ## Look at policies of those without any neighbor
        effect <- 0
        if(length(no_neighbor) == 0) num_treated_neighbors <- W_0%*%policies_noneighb
        
        if(length(no_neighbor) > 0) num_treated_neighbors <- W_0[-no_neighbor, -no_neighbor]%*%policies_noneighb
        
        for(i in 1:length(neighbors)){
          num_treated_neighb <- num_treated_neighbors[i]
          effect1 <- gs1[num_treated_neighb + 1 + acc]
          effect0 <- gs0[num_treated_neighb + 1 + acc]
          effect1 = ifelse(is.na(effect1), 0, effect1)
          effect0 = ifelse(is.na(effect0), 0, effect0)
          saved_effect1 =  focal_neighb[i] * (policies_noneighb[i]*effect1 + (1 - policies_noneighb[i])*effect0)
          effect <- effect + ifelse(is.na(saved_effect1),0, saved_effect1) ## look at the effect consistent with the direct effect of the policy
          acc <- acc + neighbors[i] + 1  ## Skip to next observation (after num_neighbor + 1 step you have the following unit)
          ## Store the effect (used for CI)
          if(i == 1){
            saved_effect <- focal_neighb[i] * (policies_noneighb[i]*effect1 + (1 - policies_noneighb[i])*effect0)
          } else{
            saved_effect <- c(saved_effect,(policies_noneighb[i]*effect1 + (1 - policies_noneighb[i])*effect0))
          }
        }
        if(length(no_neighbor) > 1){
          k <- 1
          for(i in no_neighbor){
            my_eff = focal_units[i] * (policies[i]*gs1_noneighb[k] + (1 - policies[i])*gs0_noneighb[k])
            effect <- effect + ifelse(is.na(my_eff),0, my_eff) 
            k <- k  + 1
          }
        }
        effect/n
        
      }
      return(list(predictions = NA, welfare =  my_effect_final))
    }
  }
  

