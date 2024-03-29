% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper.R
\name{PredictPolicies}
\alias{PredictPolicies}
\title{Predict the policy for new observations}
\usage{
PredictPolicies(opt, X)
}
\arguments{
\item{opt}{solutions to NetworkTargeting}

\item{X}{variables for targeting;}

\item{params}{list with additional parameters (see references)}
}
\value{
predictions: predicted effect on each individual;
}
\description{
This function predicts the policy for new observations
}
