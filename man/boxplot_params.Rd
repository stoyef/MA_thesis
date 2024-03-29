% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HMM_plots.R
\name{boxplot_params}
\alias{boxplot_params}
\title{Boxplot of different simulation runs}
\usage{
boxplot_params(data, name, true_value, ...)
}
\arguments{
\item{data}{Simulated values of the parameter.}

\item{name}{Name of the parameter.}

\item{true_value}{optional, true parameter value, 0 if no horizontal line should 
be drawn.}
}
\description{
Generate a boxplot of several fitted values, generated by fitting a model
to data. Includes a comparison to the true value.
}
