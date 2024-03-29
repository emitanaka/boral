\name{get.enviro.cor}

\alias{get.enviro.cor}

\docType{package}

\title{Extract covariances and correlations due to shared environmental responses}

\description{Calculates the correlation between columns of the response matrix, due to similarities in the response to explanatory variables i.e., shared environmental response}

\usage{
get.enviro.cor(object, est = "median", prob = 0.95)
}

\arguments{
  \item{object}{An object for class "boral".}

  \item{est}{A choice of either the posterior median (\code{est = "median"}) or posterior mean (\code{est = "mean"}), which are then treated as estimates and the fitted values are calculated from. Default is posterior median.}

  \item{prob}{A numeric scalar in the interval (0,1) giving the target probability coverage of the intervals, by which to determine whether the correlations are "significant". Defaults to 0.95.}   
}

\details{
In both independent response and correlated response models, where the each of the columns of the response matrix \code{y} are fitted to a set of explanatory variables given by \code{X}, the covariance and thus between two columns \eqn{j} and \eqn{j'} due to similarities in their response to the model matrix is calculated based on the linear predictors \eqn{\bm{x}^\top_i\bm{\beta}_j} and \eqn{\bm{x}^\top_i\bm{\beta}_{j'})}, where \eqn{\bm{\beta}_j} are column-specific coefficients relating to the explanatory variables.

For multivariate abundance data, the correlation calculated by this function can be interpreted as the correlation attributable to similarities in the environmental response between species. Such correlation matrices are discussed and found in Ovaskainen et al., (2010), Pollock et al., 2014.
}

\value{
A list with the following components:
\item{cor, cor.lower, cor.upper}{A set of \eqn{p \times p} correlation matrices, containing either the posterior median or mean estimate plus lower and upper limits of the corresponding 95\% highest posterior interval.}

\item{sig.cor}{A \eqn{p \times p} correlation matrix containing only the ``significant" correlations whose 95\% highest posterior interval does not contain zero. All non-significant correlations are set to zero.}

\item{cov}{A \eqn{p \times p} covariance matrix.}
}

\references{
\itemize{
\item Ovaskainen et al. (2010). Modeling species co-occurrence by multivariate logistic regression generates new hypotheses on fungal interactions. Ecology, 91, 2514-2521.

\item Pollock et al. (2014). Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). Methods in Ecology and Evolution, 5, 397-406.
}
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\seealso{
\code{\link{get.residual.cor}}, which calculates the residual correlation matrix for models involving latent variables.
}

\examples{
\dontrun{
## NOTE: The values below MUST NOT be used in a real application;
## they are only used here to make the examples run quick!!!
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)
     
library(mvabund) ## Load a dataset from the mvabund package
library(corrplot) ## For plotting correlations
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y)
p <- ncol(y)
    
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
     save.model = TRUE, mcmc.control = example_mcmc_control)

enviro.cors <- get.enviro.cor(spiderfit_nb)

corrplot(enviro.cors$sig.cor, title = "Shared response correlations", 
	type = "lower", diag = FALSE, mar = c(3,0.5,2,1), tl.srt = 45)
}
}