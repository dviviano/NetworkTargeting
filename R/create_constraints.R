

### Compute constraints for MILP program

## Pass output Y, X (matrix for policy evaluation: it must be the same of Z if m1 is NA)
## neighbors: vector with the number of neighbors per each observations
## final_gs: g_i(1, h) - g_i(0,h) per each observation ordered from h = 0 up to h = |N_i| sequentially over observations
## if final_gs = NA, computed from default model
## final_g0: the same only for g_i(0,h)
## B : constraint on betas
## m1: fitted model
## parms: parameter for gurobi
## tolerance constraint: see b constraint
## model checking: if T it checks discrepancies in the constraints .

## If missing values final gs, default method only uses the conditional mean function.

## Z is used to do predictions , X is used for policy evaluation
## Z_int: covariates used for interaction with treatment effects

## Cut_high_neighb: does not evaluate individuals with more than 15 neighbors in the objective function
## to speed up optimization
## include_reg_neighbors: include number of neighbors as covariate in the regression
## cost: relative cost of treatment
## also_include_atoms: use also atoms for estimating the policy function
## type : binomial or continuous: support of the Y
## if monotonicity == T
## First variable in X is the one for which monotonicity constraints are imposed
## The variables must be ordered increasing in X1

compute_constraints_function <- function(
    Y,X,Z, Z_int, D,W_0, final_gs =NA, final_g0=NA, final_gs_noneigh = NA,
    final_g0_noneighb = NA, m1 =NA, cost = 0, cut_high_neighb = F, B =1,
    tolerance_constraint1, tolerance_constraint2, monotonicity_c = F,
    direction = '>=', type= 'binomial', Z_int_neighbors = Z_int) {

  constraint_matrix_no_n <- NA
  neighbors <- apply(W_0 , 1, sum)
  neighbors_column = apply(W_0 , 2, sum)
  neighbors_todivide <- sapply(neighbors, function(x) ifelse(x > 0, x, 1))
  ## Generate gs

  ## Order variables
  ## First remove the ones with no neighbors, than add back at the end
  ## (y_snoneighbors, y_1, x_10, x_20, x11, x21, y2, ..., z10, z11, z12, ...., z21, ...)
  ## constraint matrix:
  ## first 2*n rows are constraints on ys
  ## second 1:4*length(final_gs) are constraints on xs
  ## last rows are constraints on zs

  ## Remove variables with no neighbors

  ## neighbors_column deal with asymmetric adjacency matrices
  no_neighb <- which(neighbors == 0 & neighbors_column == 0)
  any_neighb <- which(neighbors > 0 & neighbors_column > 0)
  also_include_atoms <- ifelse(length(no_neighb) > 1, TRUE, FALSE)

  if (length(no_neighb) > 0 & length(any_neighb) > 0) {

    W_0 <- W_0[-no_neighb, -no_neighb]
    ## Constraints on the ys of the ones with no neighbors
    n <- dim(X)[1]
    p <- dim(X)[2]

    if (also_include_atoms) {

      X_non  <- X[no_neighb,]
      ## scale from the maximum
      XX <- cbind(1, X_non)
      XX <- t(apply(XX, 1, function(h) h/max(abs(h))))

      ## C_i constant
      C <- B*apply(XX, 1, function(x) sum(abs(x)))
      ## rescale by C_i
      XX <- t(sapply(c(1:length(no_neighb)),  function (x)  XX[x,]/C[x]))

      ## Specify the constraints

      constraint_matrix_no_n <- matrix(0, nrow = 2*length(no_neighb), ncol = n + 3*length(final_gs) + p + 1  )

      ## Constraints on y
      constraint_matrix_no_n[1:length(no_neighb), (dim(constraint_matrix_no_n)[2] - p):dim(constraint_matrix_no_n)[2]] <- -XX
      constraint_matrix_no_n[(length(no_neighb) + 1):(2*length(no_neighb)), (dim(constraint_matrix_no_n)[2] - p):dim(constraint_matrix_no_n)[2]] <- -XX
      ## Pick the ys positions and fix them to 1 for the constraints
      for(i in 1:length(no_neighb)){
        constraint_matrix_no_n[c(i, i + length(no_neighb)), i] <- 1
      }

    }
    X <- X[-no_neighb,]

    neighbors <- neighbors[-no_neighb]
    neighbors_column = neighbors_column[-no_neighb]

  }

  n <- dim(X)[1]
  p <- dim(X)[2]

  ## Construct the vector c
  ## first repeat twice final_gs (g(1, h) - g(0, h))

  acc<- 1
  ## vector c
  c <- rep(0, 3*length(final_gs) + n + p +1)

  # Compute the position of the ys corresponding to the entries in final_gs
  positions <- NA
  for(i in 1:n){
    nn <- neighbors[i]
    nn_column = neighbors_column[i]
    if(nn > 0 | nn_column > 0){
      res <-  sapply(c(0:nn), function(z) i )
    } else {
      res <-  i
    }
    positions <- c(positions , res)
  }

  if (length(any_neighb)>0) {
    positions <- positions[-1]
  } else {
    positions <- 1:n
  }

  for(i in 1:n){
    ## entries for ys
    my_finals <- final_gs[positions == i]
    ## use c only to save the positions of the ys
    c[acc] <- -sum(my_finals, na.rm  = T) + exp(-10)*(sum(my_finals, na.rm = T) == 0)
    acc <- (2*neighbors[i]+2) + 1 + acc
  }
  y_s_position <- which(c!=0)
  ## Remove the - sum(finals) since it is wrong
  c[y_s_position] <- -cost
  ## entris for zs

  ### make high degree neighbots irrelevant to speed up optimization
  if(cut_high_neighb){
    which_high <- which(neighbors > 15)
    final_gs[positions %in% which_high] <-0
    final_g0[positions %in% which_high] <-0
  }

  c[(2*length(final_gs) + n + 1):(3*length(final_gs) + n)] <- final_gs
  rem_pos<- c(y_s_position, (2*length(final_gs) + n + 1):length(c))
  c[-rem_pos] <- rep(final_g0, each = 2)
  save_final_c <- c
  #c <- (10**(4))*c
  ## Constraints on ys for maximum score

  ## scale from the maximum
  XX <- cbind(1, X)
  XX <- t(apply(XX, 1, function(h) h/max(abs(h))))

  ## C_i constant
  C <- B*apply(XX, 1, function(x) sum(abs(x)))
  ## rescale by C_i
  XX <- t(sapply(c(1:n),  function (x)  XX[x,]/C[x]))

  ## Specify the constraints
  ## Consider the vector of (y_1, x_11, x_12, ... y_2, ...)








  constraint_matrix <- matrix(0, nrow = 2*length(c) - 2*p - 2, ncol = length(c) )

  ## Constraints on y
  constraint_matrix[1:n, (dim(constraint_matrix)[2] - p):dim(constraint_matrix)[2]] <- -XX
  constraint_matrix[(n + 1):(2*n), (dim(constraint_matrix)[2] - p):dim(constraint_matrix)[2]] <- -XX
  ## Pick the ys positions and fix them to 1 for the constraints
  k <- 1
  for(i in y_s_position){
    constraint_matrix[c(k, k + n), i] <- 1
    k <- k +1
  }

  ## Counts the number of neighbors for each entry of final_gs repeated twice

  neighb_count <- NA
  for(i in 1:n){
    nn <- neighbors[i]
    nn_column = neighbors_column[i]
    if(nn > 0 | nn_column > 0) {
      res <-  sapply(c(0:nn), function(z) c(nn,nn) )
    } else {
      res <- c(0,0)
    }
    neighb_count <- c(neighb_count , res)
  }
  neighb_count <- neighb_count[-1] + 1 ## add one  to each element !!! here the issue
  max_neigh <- max(neighbors)

  ### Constraints on xs
  ## pick the ones*|N_i| corresponding to the xs in the correct row of the constraint matrix

  ones <- diag(neighb_count)
  ones <- rbind(ones, ones)
  constraint_matrix[(2*length(y_s_position) + 1):(dim(constraint_matrix)[1] - length(final_gs)*2),
                    -c(y_s_position, (dim(constraint_matrix)[2] - p - length(final_gs)):dim(constraint_matrix)[2])] <- ones

  ## Pick the neighbors of ys that appear in the constraints of interest
  ## Save the positions of the xs
  pos_x <- (2*length(y_s_position) + 1):(dim(constraint_matrix)[1] - length(final_gs)*2)                 ### CHECK if problem
  ## positions of x_h1 saved in the row of the constraint matrix
  pos_x_1 <- pos_x[pos_x %% 2 > 0]
  ## positions of xh2 saved in the rows of the constraint matrix
  pos_x_2 <- pos_x[pos_x %% 2 == 0]

  ## pick the neighbors of each element i

  for(i in 1:n){
    pick_neigh <- W_0[i,]
    if(sum(pick_neigh) > 0){
      neighb_pos <- y_s_position[pick_neigh == 1]

      ## Set -1 for the neighbors on x1 and 1 for the neighbors in x2
      constraint_matrix[pos_x_1[which(c(positions,positions) == i)], neighb_pos] <- -1
      constraint_matrix[pos_x_2[which(c(positions,positions) == i)], neighb_pos] <- 1
    } else {
      constraint_matrix[pos_x_1[which(c(positions,positions) == i)], ] <- 0
      constraint_matrix[pos_x_2[which(c(positions,positions) == i)], ] <- 0
    }
  }

  ## Constraints on zs

  ## Pick the ys for the zs of interest
  k <- 1
  for(i in 1:n){
    constraint_matrix[((dim(constraint_matrix)[1] - length(final_gs)*2) + k):(dim(constraint_matrix)[1] - length(final_gs)*2 + k  + length(positions[positions == i]) - 1),
                      y_s_position[i]] <- -1
    k <- k + length(positions[positions == i])
  }
  ## Repeat the same for the second set of constraints on zs
  k <- 1
  for(i in 1:n){
    constraint_matrix[((dim(constraint_matrix)[1] - length(final_gs)) + k):(dim(constraint_matrix)[1] - length(final_gs) + k  + length(positions[positions == i]) - 1),
                      y_s_position[i]] <- -1
    k <- k + length(positions[positions == i])
  }

  ## Pick the zs for each constraint
  ones <- diag(rep(3, length(final_gs)))
  ones <- rbind(ones, ones)
  constraint_matrix[(dim(constraint_matrix)[1] - 2*length(final_gs) + 1):dim(constraint_matrix)[1],
                    (n + 2*length(final_gs) + 1):(dim(constraint_matrix)[2] - p -1 )] <- ones

  ## Pick the xs for each constraint:
  ## First set of constraints

  acc <- 2
  for(i in 1:length(final_gs)){
    constraint_matrix[dim(constraint_matrix)[1] - 2*length(final_gs) + i, (acc):(acc + 1)] <- -1
    acc <- acc + 2
    acc <- acc + sum(y_s_position == acc)
  }

  ## Second set of constraints
  acc <- 2
  for(i in 1:length(final_gs)){
    constraint_matrix[dim(constraint_matrix)[1] - length(final_gs) + i, (acc):(acc + 1)] <- -1
    acc <- acc + 2
    acc <- acc + sum(y_s_position == acc)
  }


  ## Add constraints on xs such that at least one neighbors is picked here the constraint is that sum_h x_1h + x_2h = |N_i| + 2

  ## n constraints
  acc <- 1
  additional_x_cons_x <- matrix(0, nrow = n, ncol = length(c))
  for(i in 1:n){
    additional_x_cons_x[i, (acc + 1):(2*neighbors[i] + acc + 2)] <- 1/(neighbors[i] +2)
    acc <- 2*neighbors[i] + 3 + acc
  }
  constraint_matrix <- rbind(constraint_matrix, additional_x_cons_x)

  ## add constraints on zs such that the sum of the z is equal to y (combact against slackness)
  additional_x_cons_yz <- matrix(0, nrow = n, ncol = length(c))
  k <- 1
  for(i in 1:n){
    additional_x_cons_yz[i, y_s_position[i]] <- -1
    additional_x_cons_yz[i, (n + 2*length(final_gs) + k):(n + 2*length(final_gs) + neighbors[i] + k)] <-1
    k <- neighbors[i] + 1 + k

  }
  constraint_matrix <- rbind(constraint_matrix, additional_x_cons_yz)

  ## Construct the b vector
  ## generate h for the xs
  my_h <- 0
  for (i in 1:n) {
    nn <- neighbors[i]
    nn_column = neighbors_column[i]
    if (nn > 0 | nn_column > 0) {
      my_h <- c(my_h, c(0:nn))
    } else {
      my_h <- c(my_h, c(0))
    }
  }
  my_h = my_h[-1]
  ## Count the neighbors (without repeating them twice)

  neighb_count <- NA
  for (i in 1:n) {
    nn <- neighbors[i]
    nn_column = neighbors_column[i]
    if (nn > 0 | nn_column > 0) {
      res <- sapply(c(0:nn), function(z) nn )
    } else {
      res <- c(0)
    }
    neighb_count <- c(neighb_count , res)
  }
  neighb_count <- neighb_count[-1] + 1

  b <- rep(NA, dim(constraint_matrix)[1])
  ## b for ys
  b[1:n] <- tolerance_constraint1 ## Ensure strict inequality with tolerance_constrain = 10**(-9)
  b[(n+1):(2*n)] <- 1
  ## b for xs
  b[pos_x_1[1:(length(pos_x_1)/2)]] <- -my_h + neighb_count #- tolerance_constraint
  b[pos_x_2[1:(length(pos_x_2)/2)]] <- my_h + neighb_count #- tolerance_constraint
  b[pos_x_1[(length(pos_x_1)/2 + 1):length(pos_x_1)]] <- -my_h + tolerance_constraint2 ## Here we would like slightle larger tolerance if feasibility permits
  b[pos_x_2[(length(pos_x_2)/2 + 1):length(pos_x_2)]] <- my_h + tolerance_constraint2  ## Same here : remove comment if feasible
  ## b for zs
  b[(2*n + 4*length(final_gs) + 1):(2*n + 4*length(final_gs) + length(final_gs))] <- 0 #- tolerance_constraint ## here tolerance: remove comment if feasible
  b[(2*n + 4*length(final_gs) + length(final_gs) + 1):(2*n + 4*length(final_gs) + 2*length(final_gs))] <- -3 + tolerance_constraint2

  ## Equality constraints
  b[(2*n + 4*length(final_gs) + 2*length(final_gs) + 1):(2*n + 4*length(final_gs) + 2*length(final_gs) + n)] <- 1
  b[(2*n + 4*length(final_gs) + 2*length(final_gs) + n + 1):(2*n + 4*length(final_gs) + 2*length(final_gs) + 2*n)] <- 0 ## equality sum of z_h being equal to y
  model <- list()
  model$sense<- c(rep('>', n), rep('<=', n), rep('<=', length(my_h)*2), rep('>', length(my_h)*2), rep('<=', length(my_h)), rep('>', length(my_h)), rep('=', 2*n))

  if (monotonicity_c) {
    X_minus1 <- X[,-1]
    aa <- unique(X_minus1)
    final_con <- rep(0, dim(constraint_matrix)[2])
    for(j in 1:dim(aa)[1]){
      keep_units <- X[,-1] == aa[j,]
      keep_units <- apply(keep_units,1, all)
      if(sum(keep_units) > 1){
        con <- matrix(0, nrow = sum(keep_units) - 1, ncol = dim(constraint_matrix)[2])
        keep_units <- which(keep_units == T)
        for(k in 1:(length(keep_units) - 1)){
          con[k, y_s_position[keep_units[k]]] <- 1
          con[k, y_s_position[keep_units[k + 1]]] <- -1
          final_con <- rbind(final_con, con)
        }
      }
    }
    final_con <- final_con[-1,]
    constraint_matrix <- rbind(constraint_matrix, final_con)
    model$sense<- c(model$sense, rep(direction, nrow(final_con)))
    b <- c(b, rep(0, nrow(final_con)))
  }

  return(list(b = b, model = model, constraint_matrix = constraint_matrix,
              no_neighb = no_neighb, any_neighb = any_neighb, c = c, XX = XX,
              y_s_position = y_s_position,
              positions = positions,
              pos_x = pos_x,
              pos_x_1 = pos_x_1,
              pos_x_2 = pos_x_2,
              neighbors = neighbors,
              final_gs_noneigh = final_gs_noneigh,
              my_h = my_h,
              final_g0_noneighb = final_g0_noneighb,
              constraint_matrix_no_n = constraint_matrix_no_n))
}



