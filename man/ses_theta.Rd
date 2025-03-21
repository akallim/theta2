% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theta_function.R
\name{ses_theta}
\alias{ses_theta}
\title{Standardized effect sizes for theta and zeta indices for phylogenetic signal}
\usage{
ses_theta(cat, tree, matching = "complete", nrep)
}
\arguments{
\item{cat}{A dataframe containing data for the categorical trait. The data preferably should be in different columns for each class of the trait and 0 indicating the species does not have this class, while 1 that it does. A species could have more than one categories.}

\item{tree}{The phylogenetic tree in phylo format. The tree must include all species in the trait dataset}

\item{matching}{The similarity criterion to be used "complete" for assuming two species as identical if they match in all classes, and "partial" if we assume that species are considered similar if they share at least one class of the trait}

\item{nrep}{Number of randomizations for estimating the distribution of theta and zeta under random expectations}
}
\value{
The observed values of theta and zeta indices. The values are between 0 and 1, and the smaller the value the stronger the signal. Also by comparison to null model it estimates the Standardized Effect Sizes and the significance level for how likely are the observed values compared to random expectations.
}
\description{
Standardized effect sizes for theta and zeta indices for phylogenetic signal
}
\examples{
data(tree)
data(cat)
ses_theta(cat, tree, "complete", 500)
ses_theta(cat, tree, "partial", 500)
}
