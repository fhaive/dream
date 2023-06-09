% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_avg_shortest_path.R
\name{get_avg_shortest_path}
\alias{get_avg_shortest_path}
\title{Creates a squared matrix of average shortes paths between couples of drugs. Note that each drug may have multiple targets, so the average shortest paths of all the targets of a certain drug is computed.}
\usage{
get_avg_shortest_path(Graph, drug_target_df, drugs_col, targets_col)
}
\arguments{
\item{Graph}{An object of class "igraph".}

\item{drug_target_df}{Dataframe of at least two columns containing drug names and corresponding targets.}

\item{drugs_col}{Column index of the dataframe drug_target_df where the drugs are stored.}

\item{targets_col}{Column index of the dataframe drug_target_df where the drugs are stored.}
}
\value{
A squared matrix of average shortest paths between couples of drug targets.
}
\description{
Creates a squared matrix of average shortes paths between couples of drugs. Note that each drug may have multiple targets, so the average shortest paths of all the targets of a certain drug is computed.
}
