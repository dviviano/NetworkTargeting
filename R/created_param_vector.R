create_params = function(params){

  DEFAULT_PARAMS <- list(
    # newly introduced params relevant to sparse computation setup
    solver="gurobi",
    verbose=F,
    # already existing parameters
    ## if cross fitting = T, m1 needs to be a list of models for each obs
    m1 = NA,
    # remove effects on units with high degree from computations for increasing speed
    cut_high_neighb = F,
    # tolerance constraints for the MILP program (do not set less than 10**(-6))
    tolerance_constraint1 = 10**(-3),
    tolerance_constraint2 = 10**(-3),
    ## Upper and lower bound for the parameters
    B = 1, low_b = - 1,
    ## list of parameters to pass to gurobi
    params_gurobi = NA,
    ## introduce monotonicity constraint in the optimization (default FALSE)
    monotonicity_c = F,
    ## maximum time for the partition used to cut the network (active if MILP_graph_clustering is selected)
    maxtime_partitions = 100,
    ## size of each sub-partition of the network (active if MILP_graph_clustering is selected)
    length_small_samples = 50,
    ## use double robust estimation
    doubly_robust = T,
    ## penalized method for estimating the propensity score
    penalized= F,
    ## number of strata for propensity score estimation (aggregate individuals with similar percentage of treated neighbors)  If NA each number of neighbor correspond to a different treatment status
    num_strata = 4,
    ## trimming on the propensity score (extrapolation is performed)
    trimming = 0.01,
    ## if not NA, a matrix with nrow = dim(X)[2] + 1 indicating upper bound (left column) and lower bound (right column) for each parameter
    additional_param_constraints = NA,
    ## return the gurobi model without computations
    model_return = F,
    ## names of the column of the model passed to a random forest (active only if m1 != 0 and is a random forest)
    names_rf = NA,
    ## print checking values for model stability
    print_checking = F,
    ## return the estimated conditional mean functions
    return_m1 = F,
    ## upper and lower bounds on the outcme variable
    bound_outcomes = c(-Inf, Inf),
    ## removed the subpopulation that is trimmed
    trimmed_subpopulation = F,
    ## minimum number of treated units
    min_treated = 0,
    # if method = 'MILP_graph_clustering', list of indexes with each element containing the index of each cluster
    partitions  = NA,
    ## number of cores for computations
    numcores = 1,
    ## cross_fitting: boolean note: m1 if not NA has to be a list with models for each observation (done via mincut)
    cross_fitting = F,
    ## when computing partitions of the graph for cross fitting
    ## use a slack parameter to find approximately equal sized
    ## partitions up to +- slack number of obs
    slack = 10,
    ## n_folds: number of folds for cross fitting
    n_folds = 10,

    # added params for surrogate losses
    lambda = 1,
    sigma = 1,
    lr = 1e-2,
    epochs = 100,
    cross_validation = list(),
    min_size_stop_cv = 20
  )

  new_params = modifyList(DEFAULT_PARAMS, params)

  return(new_params)
}
