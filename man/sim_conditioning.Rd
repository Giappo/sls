% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sls_sim_utils.R
\name{sim_conditioning}
\alias{sim_conditioning}
\title{Check if the conditioning is fulfilled}
\usage{
sim_conditioning(data, l_2, cond)
}
\arguments{
\item{data}{contains all the information about the simulated process}

\item{l_2}{the matrix containing the information about how the subclades are
nested into the main clade. See \link{sim_get_standard_l_2} for more info.}

\item{cond}{type of conditioning:
\itemize{
  \item cond = 0 no conditiong;
  \item cond = 1 conditions on the survival of crown descendents;
  \item cond = 2 conditions on the survival of subclade and
  on the other crown descendents in the main clade;
  \item cond = 3 conditions on the survival of the subclade and
  both crown descendents in the main clade;
}}
}
\value{
a boolean
}
\description{
sls_sim module
}
\author{
Giovanni Laudanno
}
