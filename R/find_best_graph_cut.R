## Deal with high-dimensional graph via maximum cuts with recursive partitions



min_cut_constrained_MILP <- function(W, maxtime = 100, params=NA, numcores = 1,  focal_units = rep(1, dim(W)[1]), verbose = F){

  n <- sum(focal_units)
  model  <- list()
  model$obj <- 2*apply(W,1, sum)
  model$Q <- -2*W
  model$A<- rbind(focal_units, focal_units)
  model$modelsense<-'min'
  ## tolerance_constraint enters here
  model$sense<- c('>=', '<=')
  model$vtype<- rep('B', n)
  model$rhs <- c(floor(n/3), floor(n/2))
  # Put bounds on the parameter space (If commented, parameter space = real line)
  model$ub<- c(rep(1,n))
  model$lb<- c(rep(0,n))
  ## Run gurobi only on one thread by default
  if(is.na(params)){
    params<- list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit = maxtime, BarConvTol = exp(-2),
                  Threads = numcores, Disconnected=0,Heuristics=0,Cuts=0, OutputFlag = as.numeric(verbose))}

  result<- gurobi::gurobi(model, params = params)
  return(result$x)
}

## Find the best cut in the graph
recursive_min_cut <- function(
    W, maxtime = 100, params=NA, k = 1, num_cut = 3, vec, numcores = 1,
    cut_function = function(x, my_focal_units) min_cut_constrained_MILP(x, maxtime, params, numcores, focal_units = my_focal_units, verbose = F), minimum_size_stopping_condition = 30, focal_units = rep(1, dim(W)[1])){

  cuts <- cut_function(W, my_focal_units = focal_units)
  if(k == num_cut | dim(W)[1] < minimum_size_stopping_condition){list(vec[cuts == 1], vec[cuts == 0])}else {

    ## To do this function
    rec1 <- recursive_min_cut(W = W[cuts == 1, cuts == 1], maxtime, params, k + 1, num_cut, vec[cuts == 1],
                              numcores = numcores, cut_function = cut_function,
                              minimum_size_stopping_condition = minimum_size_stopping_condition,
                              focal_units = focal_units[cuts == 1])
    rec0 <- recursive_min_cut(W[cuts == 0, cuts == 0], maxtime, params, k + 1, num_cut, vec[cuts== 0],
                              numcores = numcores, cut_function = cut_function,
                              minimum_size_stopping_condition = minimum_size_stopping_condition,
                              focal_units = focal_units[cuts == 0])
    list(rec1, rec0)

  }
}




## Function to aggregate solutions from the MILP clustering
compute_aggregate_solution <- function(sol, X, indexes, impose_max_treated,
                                       maximum_time, params, numcores = 12, verbose = F){

  all_indexes = unlist(indexes)
  sample_sizes = sol[,dim(X)[2] + 3]
  k = 1
  acc = 1
  policies = rep(NA, dim(X)[1])
  for(j in sample_sizes){
    units = all_indexes[k:(k + j - 1)]
    betas = sol[acc, 1:(dim(X)[2] + 1)]
    policies[units] <- apply(cbind(1, X[units,]), 1, function(x) ifelse(x%*%as.numeric(betas) > 0, 1, 0))
    k = k + j
    acc = acc + 1
  }
  if(is.na(impose_max_treated)) impose_max_treated = dim(X)[1]
  g_i = 2 * policies - 1
  my_result <- Est_max_score(X,B=1,g_i,params = params, tolerance_constraint = 10**(-8), scale_Y = F,
                             max_treated_units = impose_max_treated, maxtime =  maximum_time,
                             numcores = numcores, verbose = verbose)
  objective = mean(apply(cbind(policies, my_result[[1]]), 1, function(x) x[1] == x[2]))
  return(list(objective = objective, all_res = my_result))

}


## Cuts for cross fitting



## Same as the max cut above but here also uses focal units
opposite_min_cut_MILP <- function(W, slack = 10,  maxtime = 100, params=NA, focal_units = rep(1, dim(W)[1]), verbose = F){

  n <- dim(W)[1]
  n_focal = sum(focal_units)
  model  <- list()
  model$obj <- 2*apply(W,1, sum)
  model$Q <- -2*W
  model$A<- rbind(as.numeric(focal_units), as.numeric(focal_units))
  model$modelsense<-'min'
  ## tolerance_constraint enters here
  model$sense<- c('>=', '<=')
  model$vtype<- rep('B', n)
  model$rhs <- c(n_focal/2 - slack, n_focal/2 + slack)
  # Put bounds on the parameter space (If commented, parameter space = real line)
  model$ub<- c(rep(1,n))
  model$lb<- c(rep(0,n))
  model$quadcon <- list()
  ## Run gurobi only on one thread by default
  if(is.na(params)){
    params<- list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9,
                  TimeLimit = maxtime, BarConvTol = exp(-2), Threads = 1,
                  Disconnected=0,Heuristics=0,Cuts=0, OutputFlag = as.numeric(verbose))}

  result<- gurobi::gurobi(model, params = params)
  return(result$x)
}


helper_unlist_list = function(my_list){

  unlisted = unlist(my_list, recursive = F)
  if(is.list(unlisted[[1]]) == F) return(unlisted)
  helper_unlist_list(unlisted)
}

## Compute the opposite min cut for one and two degree
## return a list with two vectors, for degree one and two respectively

opposite_min_cut_MILP_1_and_2_degree = function(W, num_folds = 5, slack = 10, maxtime = 100, params=NA,
                                                minimum_size_stopping_condition = 10, focal_units = rep(1, dim(W)[1]), verbose = F){

  W1 = W
  W2 = W%*%W
  W2 = W + W2
  W2[W2 > 1] = 1
  diag(W) = 0
  diag(W2) = 0

  my_cuts1 = recursive_min_cut(W, maxtime, params, k = 1, num_cut = floor(sqrt(num_folds)), vec = c(1:dim(W)[1]), numcores = 1, cut_function =
                      function(x, my_focal_units) opposite_min_cut_MILP(x, slack = slack, maxtime = maxtime, params = params, focal_units = my_focal_units, verbose = verbose),
                      minimum_size_stopping_condition = 5, focal_units = focal_units)

  my_cuts2 = recursive_min_cut(W2, maxtime, params, k = 1, num_cut = floor(sqrt(num_folds)), vec = c(1:dim(W)[1]), numcores = 1, cut_function =
                                 function(x, my_focal_units) opposite_min_cut_MILP(x, slack = slack, maxtime = maxtime, params = params, focal_units = my_focal_units, verbose = verbose),
                               minimum_size_stopping_condition = minimum_size_stopping_condition, focal_units = focal_units)

  vector_cuts1 = helper_unlist_list(my_cuts1)
  num_groups = length(vector_cuts1)
  vector_groups1 = rep(NA, dim(W)[1])
  for(j in 1:length(vector_cuts1)){
    vector_groups1[vector_cuts1[[j]]] = j
  }

  vector_cuts2 = helper_unlist_list(my_cuts2)
  num_groups = length(vector_cuts2)
  vector_groups2 = rep(NA, dim(W)[1])
  for(j in 1:length(vector_cuts2)){
    vector_groups2[vector_cuts2[[j]]] = j
  }

  return(list(cut1 = vector_groups1, cut2 = vector_groups2))
  }


# Optimally compute subgraphs and determine membership
#
# Perform recursive subgraph max-cut optimization and return vector of IDs for members of each subgraph.
# This function embeds the recursive_min_cut function into a wrapper

subgraph_membership <- function(
    W_0, cut_function,
    maxtime = 100, k = 1, num_cut = 3,
    minimum_size_stopping_condition = 30,
    numcores = 1, focal_units = rep(1, dim(W_0)[1]))
{

  subgraph_tree <- recursive_min_cut(
    W = W_0, vec = 1:dim(W_0)[1],
    cut_function = \(x, my_focal_units) min_cut_constrained_MILP(
      x, maxtime, params=NA, numcores, focal_units = my_focal_units, verbose = F),
    minimum_size_stopping_condition = minimum_size_stopping_condition,
    focal_units = focal_units)

  leaves <- lapply(rapply(subgraph_tree, function(x)  paste(x, collapse = " ")), \(x) as.numeric(strsplit(x, " ")[[1]]))

  sapply(1:dim(W_0)[1], \(z) which(sapply(leaves, function(y) z %in% y)))

}

