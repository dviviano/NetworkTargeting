% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper.R
\name{NetworkTargeting}
\alias{NetworkTargeting}
\title{Estimate optimal policy}
\usage{
NetworkTargeting(
  Y,
  X,
  Z,
  Z_int,
  D,
  W_0,
  W_0_train,
  cost = 0,
  method = "MILP",
  focal_units = rep(1, length(Y)),
  model = "penalized_regression",
  family = "gaussian",
  maximum_treated = NA,
  maximum_time,
  params = list()
)
}
\arguments{
\item{Y}{vector of outcomes. Y must contain missing values for those units whose neighbors are unknown but who are neighbors of focal points;}

\item{X}{variables for targeting;}

\item{Z}{variables used for estimating the conditional mean function and the propensity score;}

\item{Z_int}{variables included in the regression that interact with individual and neighbors' covariates;}

\item{D}{binary treatment assignments;}

\item{W_0}{adjacency matrix with 0/1 entries;}

\item{W_0_train}{adjacency matrix for training. Default is that W_0_train equals W_0.}

\item{cost}{cost of the treatment, default is zero;}

\item{method}{either "MILP" , "MILP_graph_clustering", "surrogate", or "MILP_sparse". The first return exact solutions but is computationally demanding; the second return approximate solutions obtained by first partitioning the graph into approximately independent components and then computing the solutions in each sub-graph. It both reports each separate solution and the best solution out of sample. 'surrogate' is strongly recommended for graphs with n >= 500, it solves an approximate surrogate objective function. 'MILP_sparse' solves an exact program, but it uses less RAM (using spars matrices).}

\item{focal_units}{vector of length n, indicating with 1 whether an individual is focal, and zero otherwise. An individual should be focal if all of her neighbors are observed. By default, all individuals are focal.}

\item{model}{model used for estimation. Any of 'penalized_regression', 'linear_regression', 'random_forest';}

\item{family}{family of the outcome, either gaussian or binary;}

\item{maximum_treated}{maximum number of treated units. NA indicates that no constraint is imposed;}

\item{maximum_time}{time limit;}

\item{params}{list with additional parameters (see references)}
}
\value{
final_coef: coefficients of the policy;

policies: individual treatment assignments

solutions_clustering: matrix with solutions when MILP_graph_clustering is selected

partitions: list of indexes with partitions of the graph if MILP_graph_clustering is selected
}
\description{
This function computes the optimal policy from a network.
The policy is the maximum score function.
Estimation uses the Mixed integer program solver Gurobi. Installation is available here: https://www.gurobi.com/ - free access for academic users.
}
