
NULL





#' Correlation structure for latent variables
#' 
#' This help file provides more information how (non-independence) correlation
#' structures can be assumed for latent variables.
#' 
#' In the main boral function, when latent varaibles are included, the default
#' option is to assume that the latent variables are independent across the
#' rows (sites) of the response matrix i.e., \code{lv.type = "independent"}.
#' That is, \eqn{\bm{z}_i \sim N(\bm{0},\bm{I}_d)} where \code{d = num.lv}.
#' This is useful when we want to model between species correlations (is a
#' parsimonious manner), but it does make an assumption that sites are
#' independent.
#' 
#' If one \emph{a-priori} believes that the sites are, in fact, correlated
#' e.g., due to spatial correlation, and that it cannot be sufficiently well
#' accounted for by row effects (see comment below), then we can account for
#' this by assuming a non-independence correlation structure for the the latent
#' variables across sites. Note however we continue to assume that the \eqn{d}
#' latent variables are still independent of one another. That is, if we let
#' \eqn{\bm{z}_i = (z_{i1}, \ldots, z_{id})}, then we assume that for \eqn{l =
#' 1,\ldots,d},
#' 
#' \deqn{(z_{1l}, z_{2l}, \ldots, z_{nl}) \sim N(\bm{0}, \bm{\Sigma}),}
#' 
#' where \eqn{\bm{\Sigma}} is some correlation matrix. When \eqn{\bm{\Sigma} =
#' \bm{I}_n} then we are back in the independence case. However, if we allow
#' for the off-diagonals to be non-zero, then we the latent variables to be
#' correlated, \eqn{\Sigma_{ij} = Cov(z_{il}, z_{jl})}. This in turn induces
#' correlation across sites and species i.e., two species at two different
#' sites are now correlated because of the correlation across sites.
#' 
#' While there are fancier structures and attempts at accounting for
#' correlations between sites (Cressie and Wikle, 2015), in boral we assume
#' relatively simple structures. Specifically, we can assume that sites further
#' away are less correlated, and so \eqn{\Sigma} can be characterized based on
#' a distance matrix \code{distmat} and associated spatial covariance
#' parameters which require estimation. Indeed, such simple spatial latent
#' variable models have become rather popular in community ecology of late, at
#' least as a first attempt at accounting for spatial (and also temporal)
#' correlation e.g., Thorson et al., (2015, 2016); Ovaskainen et al., (2017).
#' 
#' At the moment, several correlation structures are permitted. Let
#' \eqn{D_{ij}} denote the distance between site \eqn{i} and j i.e., entry
#' \eqn{(i,j)} in \code{distmat}. Also, let \eqn{(\vartheta_1,\vartheta_2)}
#' denote the two spatial covariance parameters (noting that the second
#' parameter is not required for some of structures). Then we have: 1)
#' \code{lv.type = "exponential"} such that \eqn{\Sigma_{ij} =
#' \exp(-D_{ij}/\vartheta_1)}; 2) \code{lv.type = "squared.exponential"}, such
#' that \eqn{\Sigma_{ij} = \exp(-D_{ij}/\vartheta_1^2)}; 3) \code{lv.type =
#' "power.exponential"}, such that \eqn{\Sigma_{ij} =
#' \exp(-(D_{ij}/\vartheta_1)^{\vartheta_2})} where \eqn{\vartheta_1 \in (0,2]}
#' ; 4) \code{lv.type = "spherical"}, such that \eqn{(D_{ij} < \vartheta_1)*(1
#' - 1.5*D_{ij}/\vartheta_1 + 0.5*(D_{ij}/\vartheta_1)^3)}. We refer the reader
#' to the \code{geoR} and the function \code{cov.spatial} for more, simple
#' information on spatial covariance functions (Ribeiro Jr and Diggle, 2016).
#' 
#' It is important to keep in mind that moving away from an independence
#' correlation structure for the latent variables \emph{massively} increases
#' computation time for MCMC sampling (and indeed any estimation method for
#' latent variable models). Given JAGS is not the fastest of methods when it
#' comes to MCMC sampling, then one should be cautious about moving away from
#' indepndence. For example, if you \emph{a-priori} have a nested experimental
#' design which is inducing spatial correlation, then it is much faster and
#' more effective to include (multiple) row effects in the model to account for
#' this spatial correlation instead.
#' 
#' @name about.lvs
#' @docType package
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{boral}} for the main boral fitting function.
#' @references \itemize{ \item Cressie, N. and Wikle, C. K. (2015) Statistics
#' for spatio-temporal data. John Wiley & Sons.
#' 
#' \item Ovaskainen, O., Tikhonov, G. Norberg, A. Guillaume Blanchet, F. Duan,
#' L. Dunson, D. Roslin, T., and Abrego, N. (2017). How to make more out of
#' community data? A conceptual framework and its implementation as models and
#' software. Ecology Letters, 20, 561-576.
#' 
#' \item Ribeiro Jr, P. J., and Diggle P. J., (2016). geoR: Analysis of
#' Geostatistical Data. R package version 1.7-5.2.
#' \url{https://CRAN.R-project.org/package=geoR}.
#' 
#' \item Thorson, J. T., Ianelli, J. N., Larsen, E. A., Ries, L., Scheuerell,
#' M. D., Szuwalski, C., and Zipkin, E. F. (2016). Joint dynamic species
#' distribution models: a tool for community ordination and spatio-temporal
#' monitoring. Global Ecology and Biogeography, 25, 1144-1158
#' 
#' \item Thorson, J. T., Scheuerell, M. D., Shelton, A. O., See, K. E., Skaug,
#' H. J., and Kristensen, K. (2015). Spatial factor analysis: a new tool for
#' estimating joint species distributions and correlations in species range.
#' Methods in Ecology and Evolution, 6, 627-63 }
#' @examples
#' 
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' X <- scale(spider$x)
#' n <- nrow(y)
#' p <- ncol(y)
#' 
#' ## NOTE: The two examples below and taken directly from the boral help file
#' 
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#' 
#' \dontrun{
#' ## Example 2d - model with environmental covariates and 
#' ##  two structured latent variables using fake distance matrix
#' fakedistmat <- as.matrix(dist(1:n))
#' spiderfit_lvstruc <- boral(y, X = X, family = "negative.binomial", 
#'     lv.control = list(num.lv = 2, type = "powered.exponential", 
#'     distmat = fakedistmat))
#' 
#' summary(spiderfit_lvstruc)
#' 
#'     }
#' 
#' 
NULL





#' Stochastic search variable selection (SSVS) in boral
#' 
#' This help file provides more information regarding the implementation of the
#' stochastic search variable selection (SSVS, George and McCulloch, 1993) as
#' implemented in the boral package.
#' 
#' Stochastic search variable selection (SSVS, George and McCulloch, 1993) is a
#' approach for model selection, which is applicable specifically to the
#' Bayesian MCMC framework. As of boral version 1.5, SSVS is implemented in two
#' ways.
#' 
#' \bold{SSVS on coefficients in \code{X}:} SSVS is implemented on the
#' column-specific coefficients \eqn{\bm{\beta}_j}. Basically, SSVS works by
#' placing a spike-and-slab priors on these coefficients, such that the spike
#' is a narrow normal distribution concentrated around zero and the spike is a
#' normal distribution with a large variance.
#' 
#' \deqn{\rho(\beta) = I_{\beta = 1}\times\mathcal{N}(0,\sigma^2) + (1-I_{\beta
#' = 1})\times \mathcal{N}(0,g*\sigma^2),}
#' 
#' where \eqn{\sigma^2} is determined by \code{prior.control$hypparams[3]},
#' \eqn{g} is determined by \code{ssvs.g}, and \eqn{I_{\beta = 1} = P(\beta =
#' 1)} is an indicator function representing whether coefficient is included in
#' the model. It is given a Bernoulli prior with probability of inclusion 0.5.
#' After fitting, the posterior probability of \eqn{\beta} being included in
#' the model is returned based on posterior mean of the indicator function
#' \eqn{I_{\beta = 1}}. Note this is NOT the same as a \emph{p}-value seen in
#' maximum likelihood estimation: a \emph{p}-value provides an indication of
#' how much evidence there is against the null hypothesis of \eqn{\beta = 0},
#' while the posterior probability provides a measure of how likely it is for
#' \eqn{\beta \ne 0} given the data.
#' 
#' SSVS can be applied at a grouped or individual coefficient level, and this
#' is governed by \cr \code{prior.control$ssvs.index}: \itemize{ \item For
#' elements of \code{ssvs.index} equal to -1, SSVS is not applied on the
#' corresponding covariate of the \code{X}.
#' 
#' \item For elements equal to 0, SSVS is applied to each individual
#' coefficients of the corresponding covariate in \code{X}. That is, the fitted
#' model will return posterior probabilities for this covariate, one for each
#' column of \code{y}.
#' 
#' \item For elements taking positive integers 1,2,..., SSVS is applied to each
#' group of coefficients of the corresponding covariate in \code{X}. That is,
#' the fitted model will return a single posterior probability for this
#' covariate, indicating whether this covariate should be included for all
#' columns of \code{y}; see O'Hara and Sillanpaa (2009) and Tenan et al. (2014)
#' among many others for an discussion of Bayesian variable selection methods.
#' }
#' 
#' Note the last application of SSVS allows multiple covariates to be selected
#' \emph{simultaneously}. For example, suppose \code{X} consists of five
#' columns: the first two columns are environmental covariates, while the last
#' three correspond to quadratic terms of the two covariates as well as their
#' interaction. If we want to "test" whether any quadratic terms are required,
#' then we can set \cr \code{prior.control$ssvs.index = c(-1,-1,1,1,1)}, so a
#' single posterior probability of inclusion is returned for the last three
#' columns of \code{X}.
#' 
#' Finally, note that summaries such as posterior medians and HPD intervals of
#' the coefficients, as well as performing residual analysis, from a fitted
#' model that has implemented SSVS may be problematic because the posterior
#' distribution is by definition multi-modal. It may be advisable instead to
#' separate out their application of SSVS and posterior inference.
#' 
#' \bold{SSVS on trait coefficients:} If traits are included in boral, thereby
#' leading to a fourth corner model (see \code{\link{about.traits}} for more
#' details on this type of model), SSVS can also be performed on the associated
#' trait coefficients. That is, in such model we have
#' 
#' \deqn{\beta_{0j} \sim N(\kappa_{01} + \bm{traits}^\top_j\bm{\kappa}_1,
#' \sigma^2_1)}
#' 
#' for the column-specific intercepts, and
#' 
#' \deqn{\beta_{jk} \sim N(\kappa_{0k} + \bm{traits}^\top_j\bm{\kappa}_k,
#' \sigma^2_k)}
#' 
#' for \eqn{k = 1,\ldots,d} where \code{d = ncol(X)}. Then if the a particular
#' index in the argument \cr \code{prior.control$ssvs.traitsindex} is set to 0,
#' SSVS is performed on the corresponding element in \eqn{\bm{\kappa}_1} or
#' \eqn{\bm{\kappa}_k}. For example, suppose \code{which.traits[[2]] = c(2,3)},
#' meaning that the \eqn{\beta_{j1}}'s are drawn from a normal distribution
#' with mean depending only on the second and third columns of \code{traits}.
#' Then \cr \code{prior.control$ssvs.traitsindex[[2]] = c(0,1)}, then a
#' spike-and-slab prior is placed on the first coefficent in
#' \eqn{\bm{\kappa}_2}, while the second coefficient is assigned the
#' ``standard" prior governed by the \code{prior.control$hypparams}. That is,
#' SSVS is performed on the first but not the second coefficient in
#' \eqn{\bm{\kappa}_2}.
#' 
#' Please keep in mind that because boral allows the user to manually decide
#' which traits drive which covariates in \code{X}, then care must be taken
#' when setting up both \code{which.traits} and \cr
#' \code{prior.control$ssvs.traitsindex}. That is, when supplied then both
#' objects should be lists of have the same length, and the length of the
#' corresponding vectors comprising each element in the two lists should match
#' as well e.g., \code{which.traits[[2]]} and \cr
#' \code{prior.control$ssvs.traitsindex[[2]]} should be of the same length.
#' 
#' @name about.ssvs
#' @docType package
#' @section Warnings: \itemize{ \item Summaries of the coefficients such as
#' posterior medians and HPD intervals may also be problematic when SSVS is
#' being used, since the posterior distribution will be multi-modal.
#' 
#' \item If \code{save.model = TRUE}, the raw jags model is also returned. This
#' can be quite very memory-consuming, since it indirectly saves all the MCMC
#' samples. }
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{boral}} for the main boral fitting function which
#' implementing SSVS, and \code{\link{about.traits}} for how fourth corner
#' models work before applying SSVS to them.
#' @references \itemize{ \item George, E. I. and McCulloch, R. E. (1993).
#' Variable selection via Gibbs sampling. Journal of the American Statistical
#' Association, 85, 398-409.
#' 
#' \item O'Hara, B., and Sillianpaa, M.J. (2009). A Review of Bayesian Variable
#' Selection Methods: What, How and Which. Bayesian Analysis, 4, 85-118.
#' 
#' \item Tenan, S., et al. (2014). Bayesian model selection: The steepest
#' mountain to climb. Ecological Modelling, 283, 62-69.  }
#' @examples
#' 
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' X <- scale(spider$x)
#' n <- nrow(y)
#' p <- ncol(y)
#' 
#' ## NOTE: The two examples below and taken directly from the boral help file
#' 
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#' 
#' \dontrun{
#' ## Example 3a - Extend example 2 to demonstrate grouped covariate selection
#' ## on the last three covariates. 
#' example_prior_control <- list(type = c("normal","normal","normal","uniform"), 
#'      ssvs.index = c(-1,-1,-1,1,2,3))
#' spiderfit_nb2 <- boral(y, X = X, family = "negative.binomial", 
#'     mcmc.control = example_mcmc_control, prior.control = example_prior_control)
#'      
#' summary(spiderfit_nb2) 
#' 
#' 
#' ## Example 3b - Extend example 2 to demonstrate individual covariate selection
#' ## on the last three covariates. 
#' example_prior_control <- list(type = c("normal","normal","normal","uniform"), 
#'      ssvs.index = c(-1,-1,-1,0,0,0))
#' spiderfit_nb3 <- boral(y, X = X, family = "negative.binomial", 
#'     mcmc.control = example_mcmc_control, prior.control = example_prior_control)
#' summary(spiderfit_nb3) 
#' 
#' 
#' ## Example 5a - model fitted to count data, no site effects, and
#' ## two latent variables, plus traits included to explain environmental responses
#' data(antTraits)
#' y <- antTraits$abun
#' X <- as.matrix(scale(antTraits$env))
#' ## Include only traits 1, 2, and 5
#' traits <- as.matrix(antTraits$traits[,c(1,2,5)])
#' example_which_traits <- vector("list",ncol(X)+1)
#' for(i in 1:length(example_which_traits)) 
#'      example_which_traits[[i]] <- 1:ncol(traits)
#' ## Just for fun, the regression coefficients for the second column of X,
#' ## corresponding to the third element in the list example_which_traits,
#' ## will be estimated separately and not regressed against traits.
#' example_which_traits[[3]] <- 0
#' 
#' fit_traits <- boral(y, X = X, traits = traits, 
#'     which.traits = example_which_traits, family = "negative.binomial", 
#'     mcmc.control = example_mcmc_control, save.model = TRUE)
#' 
#' summary(fit_traits)
#' 
#' 
#' ## Example 5b - perform selection on trait coefficients
#' ssvs_traitsindex <- vector("list",ncol(X)+1)
#' for(i in 1:length(ssvs_traitsindex)) 
#'      ssvs_traitsindex[[i]] <- rep(0,ncol(traits))
#' ssvs_traitsindex[[3]] <- -1
#' fit_traits <- boral(y, X = X, traits = traits, which.traits = example_which_traits, 
#'     family = "negative.binomial", mcmc.control = example_mcmc_control, 
#'     save.model = TRUE, prior.control = list(ssvs.traitsindex = ssvs_traitsindex))
#' 
#' summary(fit_traits)
#' }
#' 
#' 
NULL





#' Including species traits in boral
#' 
#' This help file provides more information regarding the how species can be
#' included to help mediate environmental responses, analogous to the so-called
#' fourth corner problem.
#' 
#' In the main boral function, when covariates \code{X} are included i.e. both
#' the independent and correlated response models, one has the option of also
#' including traits to help explain differences in species environmental
#' responses to these covariates. Specifically, when \code{traits} and
#' \code{which.traits} are supplied, then the \eqn{\beta_{0j}}'s and
#' \eqn{\bm{\beta}_j}'s are then regarded as random effects drawn from a normal
#' distribution. For the column-specific intercepts, we have
#' 
#' \deqn{\beta_{0j} \sim N(\kappa_{01} + \bm{traits}^\top_j\bm{\kappa}_1,
#' \sigma^2_1),}
#' 
#' where \eqn{(\kappa_{01},\bm{\kappa}_1)} are the regression coefficients
#' relating to the traits to the intercepts and \eqn{\sigma_1} is the error
#' standard deviation. These are now the "parameters" in the model, in the
#' sense that priors are assigned to them and MCMC sampling is used to estimate
#' them (see the next section on estimation).
#' 
#' In an analogous manner, each of the elements in \eqn{\bm{\beta}_j =
#' (\beta_{j1},\ldots,\beta_{jd})} are now drawn as random effects from a
#' normal distribution. That is, for \eqn{k = 1,\ldots,d} where \code{d =
#' ncol(X)}, we have,
#' 
#' \deqn{\beta_{jk} \sim N(\kappa_{0k} + \bm{traits}^\top_j\bm{\kappa}_k,
#' \sigma^2_k),}
#' 
#' Which traits are to included (regressed) in the mean of the normal
#' distributions is determined by the list argument \code{which.traits} in the
#' main boral function. The first element in the list applies to
#' \eqn{beta_{0j}}, while the remaining elements apply to the the
#' \eqn{\bm{\beta}_j}. Each element of \code{which.traits} is a vector
#' indicating which traits are to be used. For example, if
#' \code{which.traits[[2]] = c(2,3)}, then the \eqn{\beta_{j1}}'s are drawn
#' from a normal distribution with mean depending only on the second and third
#' columns of \code{traits}. If \code{which.traits[[2]][1] = 0}, then the
#' regression coefficients are treated as independent, i.e. the values of
#' \eqn{\beta_{j1}} are given their own priors and estimated separately from
#' each other.
#' 
#' Including species traits in the model can be regarded as a method of
#' simplifying the model: rather than each to estimates \eqn{p} sets of
#' column-specific coefficients, we instead say that these coefficients are
#' linearly related to the corresponding values of their traits (Warton et al.,
#' 2015). That is, we are using trait data to help explain
#' similarities/differences in species responses to the environment. This idea
#' has close relations to the fourth corner problem in ecology (Brown et al.,
#' 2014). Unlike the models of Brown et al. (2014) however, which treat the
#' \eqn{\beta_{0j}}'s and \eqn{\beta_{jk}}'s are fixed effects and fully
#' explained by the traits, boral adopts a random effects approach similar to
#' Jamil et al. (2013) to "soak up" any additional between species differences
#' in environmental responses not explained by traits.
#' 
#' Finally, note that from boral version 1.5, stochastic search variable
#' selection (SSVS) can now be applied to the trait coefficients
#' \eqn{\bm{\kappa}_1} and \eqn{\bm{\kappa}_k}; please see
#' \code{\link{about.ssvs}} for more details.
#' 
#' @name about.traits
#' @docType package
#' @section Warnings: \itemize{ \item \emph{No} intercept column should be
#' required in \code{traits}, as it is included automatically. }
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{boral}} for the main boral fitting function, and
#' \code{\link{about.ssvs}} for implementing SSVS on fourth corner models.
#' @references \itemize{ \item Brown, et al. (2014). The fourth-corner solution
#' - using predictive models to understand how species traits interact with the
#' environment. Methods in Ecology and Evolution 5, 344-352.
#' 
#' \item Jamil, T. et al. (2013). Selecting traits that explain
#' species-environment relationships: a generalized linear mixed model
#' approach. Journal of Vegetation Science 24, 988-1000
#' 
#' \item Warton et al. (2015). So Many Variables: Joint Modeling in Community
#' Ecology. Trends in Ecology and Evolution, 30, 766-779. }
#' @examples
#' 
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' X <- scale(spider$x)
#' n <- nrow(y)
#' p <- ncol(y)
#' 
#' ## NOTE: The two examples below and taken directly from the boral help file
#' 
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#' 
#' \dontrun{
#' ## Example 5a - model fitted to count data, no site effects, and
#' ## two latent variables, plus traits included to explain environmental responses
#' data(antTraits)
#' y <- antTraits$abun
#' X <- as.matrix(scale(antTraits$env))
#' ## Include only traits 1, 2, and 5
#' traits <- as.matrix(antTraits$traits[,c(1,2,5)])
#' example_which_traits <- vector("list",ncol(X)+1)
#' for(i in 1:length(example_which_traits)) 
#'      example_which_traits[[i]] <- 1:ncol(traits)
#' ## Just for fun, the regression coefficients for the second column of X,
#' ## corresponding to the third element in the list example_which_traits,
#' ## will be estimated separately and not regressed against traits.
#' example_which_traits[[3]] <- 0
#' 
#' fit_traits <- boral(y, X = X, traits = traits, 
#'     which.traits = example_which_traits, family = "negative.binomial", 
#'     mcmc.control = example_mcmc_control, save.model = TRUE)
#' 
#' summary(fit_traits)
#' 
#' 
#' ## Example 5b - perform selection on trait coefficients
#' ssvs_traitsindex <- vector("list",ncol(X)+1)
#' for(i in 1:length(ssvs_traitsindex)) ssvs_traitsindex[[i]] <- rep(0,ncol(traits))
#' ssvs_traitsindex[[3]] <- -1
#' fit_traits <- boral(y, X = X, traits = traits, which.traits = example_which_traits, 
#'     family = "negative.binomial", mcmc.control = example_mcmc_control, 
#'     save.model = TRUE, prior.control = list(ssvs.traitsindex = ssvs_traitsindex))
#' 
#' summary(fit_traits)
#' 
#' 
#' ## Example 6 - simulate Bernoulli data, based on a model with two latent variables, 
#' ## no site variables, with two traits and one environmental covariates 
#' ## This example is a proof of concept that traits can used to 
#' ## explain environmental responses 
#' library(mvtnorm)
#' 
#' n <- 100; s <- 50
#' X <- as.matrix(scale(1:n))
#' colnames(X) <- c("elevation")
#' 
#' traits <- cbind(rbinom(s,1,0.5), rnorm(s)) 
#' ## one categorical and one continuous variable
#' colnames(traits) <- c("thorns-dummy","SLA")
#' 
#' simfit <- list(true.lv = rmvnorm(n, mean = rep(0,2)), 
#'     lv.coefs = cbind(rnorm(s), rmvnorm(s, mean = rep(0,2))), 
#'     traits.coefs = matrix(c(0.1,1,-0.5,1,0.5,0,-1,1), 2, byrow = TRUE))
#' rownames(simfit$traits.coefs) <- c("beta0","elevation")
#' colnames(simfit$traits.coefs) <- c("kappa0","thorns-dummy","SLA","sigma")
#' 
#' simy = create.life(true.lv = simfit$true.lv, lv.coefs = simfit$lv.coefs, X = X, 
#'     traits = traits, traits.coefs = simfit$traits.coefs, family = "binomial") 
#' 
#' 
#' example_which_traits <- vector("list",ncol(X)+1)
#' for(i in 1:length(example_which_traits)) 
#'      example_which_traits[[i]] <- 1:ncol(traits)
#' fit_traits <- boral(y = simy, X = X, traits = traits, 
#'     which.traits = example_which_traits, family = "binomial", 
#'     lv.control = list(num.lv = 2), save.model = TRUE, 
#'     mcmc.control = example_mcmc_control)	
#' }
#' 
#' 
NULL





#' Bayesian Ordination and Regression AnaLysis (boral)
#' 
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_description(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_description("boral")}
#' 
#' \tabular{ll}{ Package: \tab boral\cr Type: \tab Package\cr Version: \tab
#' 0.6\cr Date: \tab 2014-12-12\cr License: \tab GPL-2\cr }
#' 
#' @name boral-package
#' @docType package
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @references \itemize{ \item Hui et al. (2014). Model-based approaches to
#' unconstrained ordination. Methods in Ecology and Evolution, 6, 399-411.
#' 
#' \item Plummer, M. (2003). JAGS: A program for analysis of Bayesian graphical
#' models using Gibbs sampling. In Proceedings of the 3rd International
#' Workshop on Distributed Statistical Computing. March (pp. 20-22).
#' 
#' \item Skrondal, A., and Rabe-Hesketh, S. (2004). Generalized latent variable
#' modeling: Multilevel, longitudinal, and structural equation models. CRC
#' Press.
#' 
#' \item Warton et al. (2015). So Many Variables: Joint Modeling in Community
#' Ecology. Trends in Ecology and Evolution, 30, 766-779.
#' 
#' \item Yi W. et al. (2013). \code{mvabund}: statistical methods for analysing
#' multivariate abundance data. R package version 3.8.4. }
#' @examples
#' 
#' ## Please see main boral function for examples. 
#' 
NULL





#' Conditional log-likelihood for a fitted model
#' 
#' Calculates the conditional log-likelihood for a set of parameter estimates
#' from a fitted model, where everything is treated as "fixed effects"
#' including latent variables, row effects, and so on.
#' 
#' For an \eqn{n x p} response matrix \code{y}, suppose we fit a model with one
#' or more latent variables. If we denote the latent variables by
#' \eqn{\bm{z}_i; i = 1,\ldots,n}, then the conditional log-likelihood is given
#' by,
#' 
#' \deqn{ \log(f) = \sum_{i=1}^n \sum_{j=1}^p \log \{f(y_{ij} | \bm{z}_i,
#' \bm{\theta}_j, \beta_{0j}, \ldots)\}, }
#' 
#' where \eqn{f(y_{ij}|\cdot)} is the assumed distribution for column \eqn{j},
#' \eqn{\bm{z}_i} are the latent variables and \eqn{\bm{\theta}_j} are the
#' coefficients relating to them, \eqn{\beta_{0j}} are column-specific
#' intercepts, and \eqn{\ldots} denotes anything else included in the model,
#' such as row effects, regression coefficients related \code{X} and
#' \code{traits}, etc...
#' 
#' The key difference between this and the marginal likelihood (see
#' \code{\link{calc.marglogLik}}) is that the conditional likelihood treats
#' everything as "fixed effects" i.e., conditions on them. These include the
#' latent variables \eqn{\bm{z}_i} and other parameters that were included in
#' the model as random effects e.g., row effects if \code{row.eff = "random"},
#' regression coefficients related to \code{X} if traits were included in the
#' model, and so on.
#' 
#' The conditional DIC, WAIC, EAIC, and EBIC returned from
#' \code{\link{get.measures}} are based on the conditional likelihood
#' calculated from this function. Additionally, \code{\link{get.measures}}
#' returns the conditional likelihood evaluated at all MCMC samples of a fitted
#' model.
#' 
#' @name calc.condlogLik
#' @docType package
#' @param y The response matrix the model was fitted to.
#' @param X The model matrix used in the model. Defaults to \code{NULL}, in
#' which case it is assumed no model matrix was used.
#' @param family Either a single element, or a vector of length equal to the
#' number of columns in \code{y}. The former assumes all columns of \code{y}
#' come from this distribution. The latter option allows for different
#' distributions for each column of \code{y}. Elements can be one of "binomial"
#' (with probit link), "poisson" (with log link), "negative.binomial" (with log
#' link), "normal" (with identity link), "lnormal" for lognormal (with log
#' link), "tweedie" (with log link), "exponential" (with log link), "gamma"
#' (with log link), "beta" (with logit link), "ordinal" (cumulative probit
#' regression).
#' 
#' Please see \code{\link{about.distributions}} for information on
#' distributions available in boral overall.
#' @param trial.size Either equal to a single element, or a vector of length
#' equal to the number of columns in y. If a single element, then all columns
#' assumed to be binomially distributed will have trial size set to this. If a
#' vector, different trial sizes are allowed in each column of y. The argument
#' is ignored for all columns not assumed to be binomially distributed.
#' Defaults to 1, i.e. Bernoulli distribution.
#' @param lv.coefs The column-specific intercept, coefficient estimates
#' relating to the latent variables, and dispersion parameters from the fitted
#' model.
#' @param X.coefs The coefficients estimates relating to \code{X} from the
#' fitted model. Defaults to \code{NULL}, in which it is assumed there are no
#' covariates in the model.
#' @param row.coefs Row effect estimates for the fitted model. The conditional
#' likelihood is defined conditional on these estimates i.e., they are also
#' treated as ``fixed effects". Defaults to \code{NULL}, in which case it is
#' assumed there are no row effects in the model.
#' @param row.ids A matrix with the number of rows equal to the number of rows
#' in \code{y}, and the number of columns equal to the number of row effects to
#' be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of
#' row \eqn{i} in \code{y} for random effect eqnj; please see the
#' \code{\link{boral}} function for details. Defaults to \code{NULL}, so that
#' if \code{row.coefs = NULL} then the argument is ignored, otherwise if
#' \code{row.coefs} is supplied then \code{row.ids = matrix(1:nrow(y), ncol =
#' 1)} i.e., a single, row effect unique to each row. An internal check is done
#' to see \code{row.coefs} and \code{row.ids} are consistent in terms of
#' arguments supplied.
#' @param offset A matrix with the same dimensions as the response matrix
#' \code{y}, specifying an a-priori known component to be included in the
#' linear predictor during fitting. Defaults to \code{NULL}.
#' @param lv Latent variables "estimates" from the fitted model, which the
#' conditional likelihood is based on. Defaults to \code{NULL}, in which case
#' it is assumed no latent variables were included in the model.
#' @param cutoffs Common cutoff estimates from the fitted model when any of the
#' columns of \code{y} are ordinal responses. Defaults to \code{NULL}.
#' @param powerparam Common power parameter from the fitted model when any of
#' the columns of \code{y} are tweedie responses. Defaults to \code{NULL}.
#' @return A list with the following components: \item{logLik}{Value of the
#' conditional log-likelihood.} \item{logLik.comp}{A matrix of the
#' log-likelihood values for each element in \code{y}, \cr such that
#' \code{sum(logLik.comp) = logLik}.}
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{calc.logLik.lv0}} to calculate the conditional/marginal
#' log-likelihood for a model with no latent variables;
#' \code{\link{calc.marglogLik}} for calculation of the marginal
#' log-likelihood;
#' @examples
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#' 
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' n <- nrow(y)
#' p <- ncol(y)
#' 
#' ## Example 1 - model with 2 latent variables, site effects, 
#' ## 	and no environmental covariates
#' spiderfit_nb <- boral(y, family = "negative.binomial", 
#'     lv.control = list(num.lv = 2), row.eff = "fixed", 
#'     save.model = TRUE, mcmc.control = example_mcmc_control)
#' 
#' ## Extract all MCMC samples
#' fit_mcmc <- get.mcmcsamples(spiderfit_nb) 
#' mcmc_names <- colnames(fit_mcmc)
#' 
#' ## Find the posterior medians
#' coef_mat <- matrix(apply(fit_mcmc[,grep("lv.coefs",mcmc_names)],
#'     2,median),nrow=p)
#' site_coef <- list(ID1 = apply(fit_mcmc[,grep("row.coefs.ID1", mcmc_names)],
#'     2,median))
#' lvs_mat <- matrix(apply(fit_mcmc[,grep("lvs",mcmc_names)],2,median),nrow=n)
#' 
#' ## Calculate the conditional log-likelihood at the posterior median
#' calc.condlogLik(y, family = "negative.binomial", 
#'     lv.coefs = coef_mat, row.coefs = site_coef, lv = lvs_mat)
#' 
#' 
#' ## Example 2 - model with no latent variables and environmental covariates
#' X <- scale(spider$x)
#' spiderfit_nb2 <- boral(y, X = X, family = "negative.binomial", 
#'     save.model = TRUE, mcmc.control = example_mcmc_control)
#' 
#' ## Extract all MCMC samples
#' fit_mcmc <- get.mcmcsamples(spiderfit_nb2) 
#' mcmc_names <- colnames(fit_mcmc)
#' 
#' ## Find the posterior medians
#' coef_mat <- matrix(apply(fit_mcmc[,grep("lv.coefs",mcmc_names)],
#'     2,median),nrow=p)
#' X_coef_mat <- matrix(apply(fit_mcmc[,grep("X.coefs",mcmc_names)],
#'     2,median),nrow=p)
#' 
#' ## Calculate the log-likelihood at the posterior median
#' calc.condlogLik(y, X = X, family = "negative.binomial", 
#'     lv.coefs =  coef_mat, X.coefs = X_coef_mat)
#' }
#' 
NULL





#' Log-likelihood for a model fitted with no latent variables
#' 
#' Calculates the log-likelihood for a set of parameter estimates from a model
#' with no latent variables. If the row effects are assumed to be random, they
#' are integrated over using Monte Carlo integration.
#' 
#' For an \eqn{n x p} response matrix \code{y}, the log-likelihood for a model
#' with no latent variables included is given by,
#' 
#' \deqn{ \log(f) = \sum_{i=1}^n \sum_{j=1}^p \log \{f(y_{ij} | \beta_{0j},
#' \alpha_i, \ldots)\}, }
#' 
#' where \eqn{f(y_{ij}|\cdot)} is the assumed distribution for column \eqn{j},
#' \eqn{\beta_{0j}} is the column-specific intercepts, \eqn{\alpha_i} is the
#' row effect, and \eqn{\ldots} generically denotes anything else included in
#' the model, e.g. row effects, dispersion parameters etc...
#' 
#' Please note the function is written conditional on all regression
#' coefficients. Therefore, if traits are included in the model, in which case
#' the regression coefficients \eqn{\beta_{0j}, \bm{\beta}_j} become random
#' effects instead (please see \code{\link{about.traits}}), then the
#' calculation of the log-likelihood does NOT take this into account, i.e. does
#' not marginalize over them!
#' 
#' Likewise if more than two columns are ordinal responses, then the regression
#' coefficients \eqn{\beta_{0j}} corresponding to these columns become random
#' effects, and the calculation of the log-likelihood also does NOT take this
#' into account, i.e. does not marginalize over them!
#' 
#' When a single \eqn{\alpha_i} random row effect is inclued, then the
#' log-likelihood is calculated by integrating over this,
#' 
#' \deqn{ \log(f) = \sum_{i=1}^n \log ( \int \prod_{j=1}^p \{f(y_{ij} |
#' \beta_{0j}, \alpha_i, \ldots)\}f(\alpha_i) d\alpha_i ), }
#' 
#' where \eqn{f(\alpha_i)} is the random effects distribution with mean zero
#' and standard deviation given by the \code{row.params}. The integration is
#' performed using standard Monte Carlo integration. This naturally extends to
#' multiple random row effects structures.
#' 
#' @name calc.logLik.lv0
#' @docType package
#' @param y The response matrix the model was fitted to.
#' @param X The model matrix used in the model. Defaults to \code{NULL}, in
#' which case it is assumed no model matrix was used.
#' @param family Either a single element, or a vector of length equal to the
#' number of columns in \code{y}. The former assumes all columns of \code{y}
#' come from this distribution. The latter option allows for different
#' distributions for each column of \code{y}. Elements can be one of "binomial"
#' (with probit link), "poisson" (with log link), "negative.binomial" (with log
#' link), "normal" (with identity link), "lnormal" for lognormal (with log
#' link), "tweedie" (with log link), "exponential" (with log link), "gamma"
#' (with log link), "beta" (with logit link), "ordinal" (cumulative probit
#' regression).
#' 
#' Please see \code{\link{about.distributions}} for information on
#' distributions available in boral overall.
#' @param trial.size Either equal to a single element, or a vector of length
#' equal to the number of columns in y. If a single element, then all columns
#' assumed to be binomially distributed will have trial size set to this. If a
#' vector, different trial sizes are allowed in each column of y. The argument
#' is ignored for all columns not assumed to be binomially distributed.
#' Defaults to 1, i.e. Bernoulli distribution.
#' @param lv.coefs The column-specific intercept, coefficient estimates
#' relating to the latent variables, and dispersion parameters from the fitted
#' model.
#' @param X.coefs The coefficients estimates relating to \code{X} from the
#' fitted model. Defaults to \code{NULL}, in which it is assumed there are no
#' covariates in the model.
#' @param row.eff Single element indicating whether row effects are included as
#' fixed effects ("fixed"), random effects ("random") or not included ("none")
#' in the fitted model. If fixed effects, then for parameter identifiability
#' the first row effect is set to zero, which analogous to acting as a
#' reference level when dummy variables are used. If random effects, they are
#' drawn from a normal distribution with mean zero and standard deviation given
#' by \code{row.params}. Defaults to "none".
#' @param row.params Parameters corresponding to the row effect from the fitted
#' model. If \cr \code{row.eff = "fixed"}, then these are the fixed effects and
#' should have length equal to the number of columns in \code{y}. If
#' \code{row.eff = "random"}, then this is the standard deviation for the
#' random effects normal distribution. If \code{row.eff = "none"}, then this
#' argument is ignored.
#' @param row.ids A matrix with the number of rows equal to the number of rows
#' in \code{y}, and the number of columns equal to the number of row effects to
#' be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of
#' row \eqn{i} in \code{y} for random effect eqnj; please see
#' \code{\link{boral}} for details. Defaults to \code{NULL}, so that if
#' \code{row.params = NULL} then the argument is ignored, otherwise if
#' \code{row.params} is supplied then \cr \code{row.ids = matrix(1:nrow(y),
#' ncol = 1)} i.e., a single, row effect unique to each row. An internal check
#' is done to see \code{row.params} and \code{row.ids} are consistent in terms
#' of arguments supplied.
#' @param offset A matrix with the same dimensions as the response matrix
#' \code{y}, specifying an a-priori known component to be included in the
#' linear predictor during fitting. Defaults to \code{NULL}.
#' @param cutoffs Common cutoff estimates from the fitted model when any of the
#' columns of \code{y} are ordinal responses. Defaults to \code{NULL}.
#' @param powerparam Common power parameter from the fitted model when any of
#' the columns of \code{y} are tweedie responses. Defaults to \code{NULL}.
#' @return A list with the following components: \item{logLik}{Value of the
#' log-likelihood} \item{logLik.comp}{A vector of the log-likelihood values for
#' each row of \code{y}, \cr such that \code{sum(logLik.comp) = logLik}.}
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{calc.marglogLik}} for calculation of the log-likelihood
#' marginalizing over one or more latent variables, and
#' \code{\link{calc.condlogLik}} for calculation of the conditional
#' log-likelihood for models where everything is treated as "fixed effects",
#' including latent variables, row effects, and so on.
#' @examples
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#' 
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' n <- nrow(y)
#' p <- ncol(y)
#' 
#' ## Example 1 - NULL model with site effects only
#' spiderfit_nb <- boral(y, family = "negative.binomial", row.eff = "fixed", 
#'     save.model = TRUE, mcmc.control = example_mcmc_control)
#' 
#' ## Extract all MCMC samples
#' fit_mcmc <- get.mcmcsamples(spiderfit_nb) 
#' mcmc_names <- colnames(fit_mcmc)
#' 
#' ## Find the posterior medians
#' coef_mat <- matrix(apply(fit_mcmc[,grep("lv.coefs",mcmc_names)],
#'     2,median),nrow=p)
#' site_coef <- list(ID1 = apply(fit_mcmc[,grep("row.coefs.ID1", mcmc_names)],
#'     2,median))
#' 
#' ## Calculate the log-likelihood at the posterior median
#' calc.logLik.lv0(y, family = "negative.binomial",
#'     lv.coefs =  coef_mat, row.eff = "fixed", row.params = site_coef)
#' 
#' 
#' ## Example 2 - Model with environmental covariates and random row effects
#' X <- scale(spider$x)
#' spiderfit_nb2 <- boral(y, X = X, family = "negative.binomial", row.eff = "random",
#'     save.model = TRUE, mcmc.control = example_mcmc_control)
#' 
#' ## Extract all MCMC samples
#' fit_mcmc <- get.mcmcsamples(spiderfit_nb2) 
#' mcmc_names <- colnames(fit_mcmc)
#' 
#' ## Find the posterior medians
#' coef_mat <- matrix(apply(fit_mcmc[,grep("lv.coefs",mcmc_names)],
#'     2,median),nrow=p)
#' X_coef_mat <- matrix(apply(fit_mcmc[,grep("X.coefs",mcmc_names)],
#'     2,median),nrow=p)
#' site.sigma <- list(ID1 = 
#'     median(fit_mcmc[,grep("row.sigma.ID1", mcmc_names)]))
#' 
#' 	
#' ## Calculate the log-likelihood at the posterior median
#' calc.logLik.lv0(y, X = spider$x, family = "negative.binomial", 
#'     row.eff = "random",lv.coefs =  coef_mat, X.coefs = X_coef_mat, 
#'     row.params = site.sigma)
#' }
#' 
NULL





#' Marginal log-likelihood for a fitted model
#' 
#' Calculates the marginal log-likelihood for a set of parameter estimates from
#' a fitted model, whereby the latent variables and random effects (if
#' applicable) are integrated out. The integration is performed using Monte
#' Carlo integration. WARNING: As of version 1.6, this function will no longer
#' be updated...use at your own peril!!!
#' 
#' For an \eqn{n x p} response matrix \code{y}, suppose we fit a model with one
#' or more latent variables. If we denote the latent variables by
#' \eqn{\bm{z}_i; i = 1,\ldots,n}, then the marginal log-likelihood is given by
#' 
#' \deqn{ \log(f) = \sum_{i=1}^n \log ( \int \prod_{j=1}^p \{f(y_{ij} |
#' \bm{z}_i, \beta_{0j}, \bm{\theta}_j, \ldots) \} f(\bm{z}_i) d\bm{z}_i), }
#' 
#' where \eqn{f(y_{ij}|\cdot)} is the assumed distribution for column \eqn{j},
#' \eqn{\beta_{0j}} are the column-specific intercepts, \eqn{\bm{\theta}_j} are
#' the column-specific latent variable coefficients, and \eqn{\ldots}
#' generically denotes anything else included in the model, e.g. row effects,
#' dispersion parameters etc... The quantity \eqn{f(\bm{z}_i)} denotes the
#' distribution of the latent variable, which is assumed to be standard
#' multivariate Gaussian. Standard Monte Carlo integration is used for
#' calculating the marginal likelihood. If \code{lv.mc = NULL}, the function
#' automatically generates a matrix as \cr \code{lv.mc <- rmvnorm(1000,
#' rep(0,num.lv))}. If there is a need to apply this function numerous times,
#' we recommend a matrix be inserted into \code{lv.mc} to speed up computation.
#' 
#' The key difference between this and the conditional likelihood (see
#' \code{calc.condlogLik}) is that the marginal likelihood treats the latent
#' variables as "random effects" and integrates over them, whereas the
#' conditional likelihood treats the latent variables as "fixed effects".
#' 
#' Please note the function is written conditional on all regression
#' coefficients. Therefore, if traits are included in the model, in which case
#' the regression coefficients \eqn{\beta_{0j}, \bm{\beta}_j} become random
#' effects instead (please see \code{\link{about.traits}}), then the
#' calculation of the log-likelihood does NOT take this into account, i.e. does
#' not marginalize over them! Likewise if more than two columns are ordinal
#' responses, then the regression coefficients \eqn{\beta_{0j}} corresponding
#' to these columns become random effects, and the calculation of the
#' log-likelihood also does NOT take this into account, i.e. does not
#' marginalize over them!
#' 
#' When a single \eqn{\alpha_i} random row effect is inclued, then the
#' log-likelihood is calculated by integrating over this,
#' 
#' \deqn{ \log(f) = \sum_{i=1}^n \log ( \int \prod_{j=1}^p \{f(y_{ij} |
#' \bm{z}_i, \beta_{0j}, \alpha_i, \ldots)\} f(\bm{z}_i) f(\alpha_i) d\alpha_i
#' ), }
#' 
#' where \eqn{f(\alpha_i)} is the random effects distribution with mean zero
#' and standard deviation given by the \code{row.params}. The integration is
#' again performed using standard Monte Carlo integration. This naturally
#' extends to multiple random row effects structures.
#' 
#' @name calc.marglogLik
#' @docType package
#' @param y The response matrix that the model was fitted to.
#' @param X The model matrix used in the model. Defaults to \code{NULL}, in
#' which case it is assumed no model matrix was used.
#' @param family Either a single element, or a vector of length equal to the
#' number of columns in \code{y}. The former assumes all columns of \code{y}
#' come from this distribution. The latter option allows for different
#' distributions for each column of \code{y}. Elements can be one of "binomial"
#' (with probit link), "poisson" (with log link), "negative.binomial" (with log
#' link), "normal" (with identity link), "lnormal" for lognormal (with log
#' link), "tweedie" (with log link), "exponential" (with log link), "gamma"
#' (with log link), "beta" (with logit link), "ordinal" (cumulative probit
#' regression).
#' 
#' Please see \code{\link{about.distributions}} for information on
#' distributions available in boral overall.
#' @param trial.size Either equal to a single element, or a vector of length
#' equal to the number of columns in y. If a single element, then all columns
#' assumed to be binomially distributed will have trial size set to this. If a
#' vector, different trial sizes are allowed in each column of y. The argument
#' is ignored for all columns not assumed to be binomially distributed.
#' Defaults to 1, i.e. Bernoulli distribution.
#' @param lv.coefs The column-specific intercept, coefficient estimates
#' relating to the latent variables, and dispersion parameters from the fitted
#' model.
#' @param X.coefs The coefficients estimates relating to \code{X} from the
#' fitted model. Defaults to \code{NULL}, in which it is assumed there are no
#' covariates in the model.
#' @param row.eff Single element indicating whether row effects are included as
#' fixed effects ("fixed"), random effects ("random") or not included ("none")
#' in the fitted model. If fixed effects, then for parameter identifiability
#' the first row effect is set to zero, which analogous to acting as a
#' reference level when dummy variables are used. If random effects, they are
#' drawn from a normal distribution with mean zero and standard deviation given
#' by \code{row.params}. Defaults to "none".
#' @param row.params Parameters corresponding to the row effect from the fitted
#' model. If \cr \code{row.eff = "fixed"}, then these are the fixed effects and
#' should have length equal to the number of columns in \code{y}. If
#' \code{row.eff = "random"}, then this is standard deviation for the random
#' effects normal distribution. If \code{row.eff = "none"}, then this argument
#' is ignored.
#' @param row.ids A matrix with the number of rows equal to the number of rows
#' in \code{y}, and the number of columns equal to the number of row effects to
#' be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of
#' row \eqn{i} in \code{y} for random effect eqnj; please see
#' \code{\link{boral}} for details. Defaults to \code{NULL}, so that if
#' \code{row.params = NULL} then the argument is ignored, otherwise if
#' \code{row.params} is supplied then \cr \code{row.ids = matrix(1:nrow(y),
#' ncol = 1)} i.e., a single, row effect unique to each row. An internal check
#' is done to see \code{row.params} and \code{row.ids} are consistent in terms
#' of arguments supplied.
#' @param offset A matrix with the same dimensions as the response matrix
#' \code{y}, specifying an a-priori known component to be included in the
#' linear predictor during fitting. Defaults to \code{NULL}.
#' @param num.lv The number of latent variables used in the fitted model. For
#' models with no latent variables, please use \code{\link{calc.logLik.lv0}} to
#' calculate the log-likelihood.
#' @param lv.mc A matrix used for performing the Monte Carlo integration.
#' Defaults to \code{NULL}, in which case a matrix is generated within the
#' function.
#' @param cutoffs Common cutoff estimates from the fitted model when any of the
#' columns of \code{y} are ordinal responses. Defaults to \code{NULL}.
#' @param powerparam Common power parameter from the fitted model when any of
#' the columns of \code{y} are tweedie responses. Defaults to \code{NULL}.
#' @return A list with the following components: \item{logLik}{Value of the
#' marginal log-likelihood.} \item{logLik.comp}{A vector of the log-likelihood
#' values for each row of \code{y}, \cr such that \code{sum(logLik.comp) =
#' logLik}.}
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{calc.condlogLik}} for calculation of the conditional
#' log-likelihood; \code{\link{calc.logLik.lv0}} to calculate the
#' conditional/marginal log-likelihood for a model with no latent variables.
#' @examples
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#' 
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' n <- nrow(y)
#' p <- ncol(y)
#'     
#' ## Example 1 - model with two latent variables, site effects, 
#' ## 	and no environmental covariates
#' spiderfit_nb <- boral(y, family = "negative.binomial", 
#'     lv.control = list(num.lv = 2), row.eff = "fixed", save.model = TRUE, 
#'     mcmc.control = example_mcmc_control)
#' 
#' ## Extract all MCMC samples
#' fit_mcmc <- get.mcmcsamples(spiderfit_nb) 
#' mcmc_names <- colnames(fit_mcmc)
#' 
#' ## Find the posterior medians
#' coef_mat <- matrix(apply(fit_mcmc[,grep("lv.coefs",mcmc_names)],
#'     2,median),nrow=p)
#' site_coef <- list(ID1 = apply(fit_mcmc[,grep("row.coefs.ID1", mcmc_names)],
#'     2,median))
#'      
#' ## Calculate the marginal log-likelihood at the posterior median
#' calc.marglogLik(y, family = "negative.binomial",
#'     lv.coefs = coef_mat, row.eff = "fixed", row.params = site_coef, 
#'     num.lv = 2)
#' 
#' 	
#' ## Example 2 - model with one latent variable, no site effects, 
#' ## 	and environmental covariates
#' spiderfit_nb2 <- boral(y, X = spider$x, family = "negative.binomial", 
#'      lv.control = list(num.lv = 2), save.model = TRUE, 
#'      mcmc.control = example_mcmc_control)
#' 
#' ## Extract all MCMC samples
#' fit_mcmc <- get.mcmcsamples(spiderfit_nb2) 
#' mcmc_names <- colnames(fit_mcmc)
#' 
#' ## Find the posterior medians
#' coef_mat <- matrix(apply(fit_mcmc[,grep("lv.coefs",mcmc_names)],
#'     2,median),nrow=p)
#' X_coef_mat <- matrix(apply(fit_mcmc[,grep("X.coefs",mcmc_names)],
#'     2,median),nrow=p)
#' 
#' ## Calculate the log-likelihood at the posterior median
#' calc.marglogLik(y, X = spider$x, family = "negative.binomial", 
#'     lv.coefs = coef_mat, X.coefs = X_coef_mat, num.lv = 2)	
#' }
#' 
NULL





#' Caterpillar plots of the regression coefficients from a fitted model
#' 
#' Constructs horizontal line plot (point estimate and HPD intervals),
#' otherwise known as "caterpillar plots", for the column-specific regression
#' coefficients corresponding to a covariate in \code{X} fitted in the model.
#' 
#' For each species (column of \code{y}), the horizontal line or "caterpillar"
#' is constructed by first marking the point estimate (posterior mean or
#' median) with an "x" symbol. Then the line is construed based on the lower
#' and upper limits of the highest posterior density (HPD) intervals as found
#' in \code{x$hpdintervals}. By default these intervals of 95\% HPD intervals.
#' To complete the plot, a vertical dotted line is drawn to denote the zero
#' value. All HPD intervals that include zero are colored gray, while HPD
#' intervals that exclude zero are colored black.
#' 
#' The graph is probably better explained by, well, plotting it using the toy
#' example below =P
#' 
#' Thanks to Robert O'Hara for suggesting and providing the original code for
#' this function.
#' 
#' @name coefsplot
#' @docType package
#' @param covname The name of one of the covariates in the fitted model. That
#' is, it must be a character vector corresponding to one of the elements in
#' \code{colnames(x)$X.coefs.median}.
#' @param x An object for class "boral".
#' @param labely Controls the labels on the y-axis for the line plot. If it is
#' not \code{NULL}, then it must be a vector either of length 1 or the same
#' length as the number of columns in the \code{y} in the fitted boral object.
#' In the former, it is treated as the y-axis label. In the latter, it is used
#' in place of the column names of \code{y} to label each line. Defaults to
#' \code{NULL}, in which the each line in the plot is labeled according to the
#' columns of \code{y}, or equivalently \code{rownames(x$X.coefs.median)}.
#' @param est A choice of either the posterior median (\code{est = "median"})
#' or posterior mean (\code{est = "mean"}), which are then used as the point
#' estimates in the lines. Default is posterior median.
#' @param ... Additional graphical options to be included in. These include
#' values for \cr \code{cex, cex.lab, cex.axis, cex.main, lwd}, and so on.
#' @return If SSVS was applied individually to each coefficient of \code{X}
#' when fitting the model, then the posterior probabilities of including the
#' specified covariate are printed out i.e., \cr those from
#' \code{x$ssvs.indcoefs.mean}.
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{caterplot} from the \code{mcmcplots} package for other,
#' sexier caterpillar plots.
#' @examples
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#' 
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' n <- nrow(y)
#' p <- ncol(y)
#' 
#' X <- scale(spider$x)
#' spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
#'     lv.control = list(num.lv = 2), mcmc.control = example_mcmc_control)
#' 
#' 
#' ## Do separate line plots for all the coefficients of X
#' par(mfrow=c(2,3), mar = c(5,6,1,1))
#' sapply(colnames(spiderfit_nb$X), coefsplot, 
#'     spiderfit_nb)
#' }
#' 
#' 
NULL





#' Simulate a Multivariate Response Matrix
#' 
#' Simulate a multivariate response matrix, given parameters such as but not
#' necessarily all of: family, number of latent variables and related
#' coefficients, an matrix of explanatory variables and related coefficients,
#' row effects, cutoffs for cumulative probit regression of ordinal responses.
#' 
#' \code{create.life} gives the user full capacity to control the true
#' parameters of the model from which the multivariate responses matrices are
#' generated from. If \code{true.lv} is supplied, then the data generation
#' mechanism is based on this. If \code{true.lv = NULL}, then the function
#' looks to \code{lv.control} to determine whether and how the true latent
#' variables are to be simulated.
#' 
#' \code{simulate} makes use of the generic function of the same name in
#' \code{R}: it takes a fitted model, treats either the posterior medians and
#' mean estimates from the model as the true parameters, and generates response
#' matrices based off that.
#' 
#' @name create.life
#' @aliases create.life simulate.boral
#' @docType package
#' @param object An object of class "boral".
#' @param nsim Number of multivariate response matrices to simulate. Defaults
#' to 1.
#' @param seed Seed for dataset simulation. Defaults to \code{NULL}, in which
#' case no seed is set.
#' @param new.lvs If \code{FALSE}, then true latent variables are obtained from
#' the \code{object}. If \code{TRUE}, then new true latent variables are
#' generated.
#' @param distmat A distance matrix required to calculate correlations across
#' sites when a non-independent correlation structure on the latent variables
#' is imposed, when \code{new.lvs = TRUE} and \code{object$lv.type !=
#' "independent"}.
#' @param est A choice of either the posterior median (\code{est == "median"})
#' or posterior mean (\code{est == "mean"}), which are then treated as
#' estimates and the fitted values are calculated from. Default is posterior
#' median.
#' @param true.lv A matrix of true latent variables. With multivariate
#' abundance data in ecology for instance, each row corresponds to the true
#' site ordination coordinates. If supplied, then simulation is based of these
#' true latent variables. If \code{NULL}, then the function looks to the
#' argument \code{lv.control} to see what to do. Defaults to \code{NULL}.
#' @param lv.coefs A matrix containing column-specific intercepts, latent
#' variable coefficients relating to \code{true.lv}, and dispersion parameters.
#' @param lv.control This argument is utilized if \code{true.lv = NULL}, in
#' which case the function uses this argument to determine how to simulate new,
#' true latent variables. A list (currently) with the following arguments:
#' \itemize{ \item \emph{num.lv:} which specifies the number of true latent
#' variables to generate. Defaults to 0.
#' 
#' \item \emph{type:} which specifies the type the correlation structure of the
#' latent variables (across sites). Defaults to independence correlation
#' structure.
#' 
#' \item \emph{lv.covparams:} which is a vector containing one or two elements
#' required if parameterizing a non-independence spatial correlation structure
#' of the latent variables.
#' 
#' \item \emph{distmat:} which a distance matrix required to calculate
#' correlations across sites when a non-independence correlation structure on
#' the latent variables is imposed.  } Please see \code{\link{about.lvs}} for
#' more information.
#' @param X An model matrix of covariates, which can be included as part of the
#' data generation. Defaults to \code{NULL}, in which case no model matrix is
#' used. No intercept column should be included in \code{X}.
#' @param X.coefs The coefficients relating to \code{X}. Defaults to
#' \code{NULL}. This argument needs to be supplied if \code{X} is supplied and
#' no \code{traits} are supplied.
#' @param traits A model matrix of species traits, which can be included as
#' part of the model. Defaults to \code{NULL}, in which case no matrix was
#' used. No intercept column should be included in \code{traits}, as it is
#' included automatically.
#' @param traits.coefs A matrix of coefficients that are used to generate "new"
#' column-specific intercepts and \code{X.coefs}. The number of rows should
#' equal to (\code{ncol(X)+1}) and the number of columns should equal to
#' (\code{ncol(traits)}+2).
#' 
#' How this argument works is as follows: when both \code{traits} and
#' \code{traits.coefs} are supplied, then new column-specific intercepts (i.e.
#' the first column of \code{lv.coefs} is overwritten) are generated by
#' simulating from a normal distribution with mean equal to \cr
#' \code{crossprod(c(1,traits), traits.coefs[1,1:(ncol(traits.coefs)-1)])} and
#' standard deviation \cr \code{traits.coefs[1,ncol(traits.coefs)]}. In other
#' words, the last column of \code{trait.coefs} provides the standard deviation
#' of the normal distribution, with the other columns being the regression
#' coefficients in the mean of the normal distribution. Analogously, new
#' \code{X.coefs} are generated in the same manner using the remaining rows of
#' \code{trait.coefs}. Please see \code{\link{about.traits}} for more
#' information.
#' 
#' It is important that highlight then with in this data generation mechanism,
#' the new column-specific intercepts and \code{X.coefs} are now random
#' effects, being drawn from a normal distribution.
#' 
#' Defaults to \code{NULL}, in conjuction with \code{traits = NULL}.
#' @param family Either a single element, or a vector of length equal to the
#' number of columns in \code{y}. The former assumes all columns of \code{y}
#' come from this distribution. The latter option allows for different
#' distributions for each column of \code{y}. Elements can be one of "binomial"
#' (with probit link), "poisson" (with log link), "negative.binomial" (with log
#' link), "normal" (with identity link), "lnormal" for lognormal (with log
#' link), "tweedie" (with log link), "exponential" (with log link), "gamma"
#' (with log link), "beta" (with logit link), "ordinal" (cumulative probit
#' regression).
#' 
#' Please see \code{\link{about.distributions}} for information on
#' distributions available in boral overall.
#' @param row.eff Single element indicating whether row effects are included as
#' fixed effects ("fixed"), random effects ("random") or not included ("none")
#' in the fitted model. If fixed effects, then for parameter identifiability
#' the first row effect is set to zero, which analogous to acting as a
#' reference level when dummy variables are used. If random effects, they are
#' drawn from a normal distribution with mean zero and standard deviation given
#' by \code{row.params}. Defaults to "none".
#' @param row.params Parameters corresponding to the row effect from the fitted
#' model. If \cr \code{row.eff = "fixed"}, then these are the fixed effects and
#' should have length equal to the number of columns in \code{y}. If
#' \code{row.eff = "random"}, then this is the standard deviation for the
#' random effects normal distribution. \cr If \code{row.eff = "none"}, then
#' this argument is ignored.
#' @param row.ids A matrix with the number of rows equal to the number of rows
#' in \code{y}, and the number of columns equal to the number of row effects to
#' be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of
#' row \eqn{i} in \code{y} for random effect eqnj; please see
#' \code{\link{boral}} for details. Defaults to \code{NULL}, so that if
#' \code{row.params = NULL} then the argument is ignored, otherwise if
#' \code{row.params} is supplied then \cr \code{row.ids = matrix(1:nrow(y),
#' ncol = 1)} i.e., a single, row effect unique to each row. An internal check
#' is done to see \code{row.params} and \code{row.ids} are consistent in terms
#' of arguments supplied.
#' @param offset A matrix with the same dimensions as the response matrix
#' \code{y}, specifying an a-priori known component to be included in the
#' linear predictor during fitting. Defaults to \code{NULL}.
#' @param trial.size Either equal to a single element, or a vector of length
#' equal to the number of columns in y. If a single element, then all columns
#' assumed to be binomially distributed will have trial size set to this. If a
#' vector, different trial sizes are allowed in each column of y. The argument
#' is ignored for all columns not assumed to be binomially distributed.
#' Defaults to 1, i.e. Bernoulli distribution.
#' @param cutoffs A vector of common common cutoffs for proportional odds
#' regression when any of \code{family} is ordinal. They should be increasing
#' order. Defaults to \code{NULL}.
#' @param powerparam A common power parameter for tweedie regression when any
#' of \code{family} is tweedie. Defaults to \code{NULL}.
#' @param manual.dim A vector of length 2, containing the number of rows
#' (\eqn{n}) and columns (\eqn{p}) for the multivariate response matrix. This
#' is a "backup" argument only required when \code{create.life} can not
#' determine how many rows or columns the multivariate response matrix should
#' be.
#' @param save.params If \code{save.params = TRUE}, then all parameters
#' provided as input and/or generated (e.g. when \code{traits} and
#' \code{traits.coefs} are supplied then \code{X.coefs} is generated
#' internally; please see \code{traits.coefs} argument above) are returned, in
#' addition to the simulated multivariate response matrix. Defaults to
#' \code{FALSE}.
#' @param ... Not used.
#' @return If \code{create.life} is used, then: 1) if \code{save.params} =
#' FALSE, a \eqn{n} by \eqn{p} multivariate response matrix is returned only,
#' 2) if \code{save.params = TRUE}, then a list containing the element
#' \code{resp} which is a \eqn{n} times \eqn{p} multivariate response matrix,
#' as well as other elements for the parameters used in the true model, e.g.
#' \code{true.lv, lv.coefs = lv.coefs, traits.coef}, is returned.
#' 
#' If \code{simulate} is used, then a three dimensional array of dimension
#' \eqn{n} by \eqn{p} by \code{nsim} is returned. The same latent variables can
#' be used each time (\code{new.lvs = FALSE}), or new true latent variables can
#' be generated each time (\code{new.lvs = TRUE}).
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{boral}} for the default function for model fitting.
#' @examples
#' 
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#' 
#' ## Example 1a - Simulate a response matrix of normally distributed data
#' library(mvtnorm)
#' 
#' ## 30 rows (sites) with two latent variables 
#' true.lv <- rbind(rmvnorm(n=15,mean=c(1,2)),rmvnorm(n=15,mean=c(-3,-1))) 
#' ## 30 columns (species)
#' lv.coefs <- cbind(matrix(runif(30*3),30,3),1)
#' 
#' X <- matrix(rnorm(30*4),30,4) 
#' ## 4 explanatory variables
#' X.coefs <- matrix(rnorm(30*4),30,4)
#' 
#' simy <- create.life(true.lv = true.lv, lv.coefs = lv.coefs, 
#'     X = X, X.coefs = X.coefs, family = "normal")
#' 
#' \dontrun{
#' fit.boral <- boral(simy, X = X, family = "normal", lv.control = list(num.lv = 2),
#'     mcmc.control = example_mcmc_control)
#' 
#' summary(fit.boral)
#' }
#' 
#' ## Example 1b - Include a nested random row effect
#' ## 30 subregions nested within six regions
#' example_row_ids <- cbind(1:30, rep(1:6,each=5))
#' ## Subregion has a small std deviation; region has a larger one
#' true.row.sigma <- list(ID1 = 0.5, ID2 = 2)
#' 
#' simy <- create.life(true.lv = true.lv, lv.coefs = lv.coefs, 
#'     X = X, X.coefs = X.coefs, row.eff = "random",
#'     row.params = true.row.sigma, row.ids = example_row_ids, family = "normal",
#'     save.params = TRUE)
#' 
#' 	
#' ## Example 1c - Same as example 1b except new, true latent variables are generated
#' simy <- create.life(true.lv = NULL, lv.coefs = lv.coefs, 
#'     X = X, X.coefs = X.coefs, row.eff = "random",
#'     row.params = true.row.sigma, row.ids = example_row_ids, family = "normal",
#'     save.params = TRUE)
#' 
#'     
#' ## Example 1d - Same as example 1a except new, true latent variables are generated
#' ## with a non-independent correlation structure using a fake distance matrix
#' makedistmat <- as.matrix(dist(1:30))
#' simy <- create.life(true.lv = NULL, lv.coefs = lv.coefs, 
#'     lv.control = list(num.lv = 2, type = "exponential", lv.covparams = 5, distmat = makedistmat),
#'     X = X, X.coefs = X.coefs, row.eff = "random",
#'     row.params = true.row.sigma, row.ids = example_row_ids, family = "normal",
#'     save.params = TRUE)
#' 
#'     
#' ## Example 2 - Simulate a response matrix of ordinal data
#' ## 30 rows (sites) with two latent variables 
#' true.lv <- rbind(rmvnorm(15,mean=c(-2,-2)),rmvnorm(15,mean=c(2,2)))
#' ## 10 columns (species)
#' true.lv.coefs <- rmvnorm(10,mean = rep(0,3)); 
#' ## Cutoffs for proportional odds regression (must be in increasing order)
#' true.ordinal.cutoffs <- seq(-2,10,length=10-1)
#' 
#' simy <- create.life(true.lv = true.lv, lv.coefs = true.lv.coefs, 
#'      family = "ordinal", cutoffs = true.ordinal.cutoffs, save.params = TRUE) 
#' 
#' \dontrun{
#' fit.boral <- boral(y = simy$resp, family = "ordinal", lv.control = list(num.lv = 2),
#'       mcmc.control = example_mcmc_control)
#' }
#' 
#' \dontrun{
#' ## Example 3 - Simulate a response matrix of count data based off
#' ## a fitted model involving traits (ants data from mvabund)
#' library(mvabund)
#' data(antTraits)
#' 
#' y <- antTraits$abun
#' X <- as.matrix(antTraits$env)
#' ## Include only traits 1, 2, and 5, plus an intercept
#' traits <- as.matrix(antTraits$traits[,c(1,2,5)])
#' ## Please see help file for boral regarding the use of which.traits
#' example_which_traits <- vector("list",ncol(X)+1)
#' for(i in 1:length(example_which_traits)) 
#'      example_which_traits[[i]] <- 1:ncol(traits)
#' 
#' fit_traits <- boral(y, X = X, traits = traits, which.traits = example_which_traits, 
#'     family = "negative.binomial", lv.control = list(num.lv = 2),
#'      mcmc.control = example_mcmc_control)
#' 
#' ## The hard way
#' simy <- create.life(true.lv = fit_traits$lv.mean, 
#'      lv.coefs = fit_traits$lv.coefs.median, X = X, 
#'      X.coefs = fit_traits$X.coefs.median, traits = traits, 
#'      traits.coefs = fit_traits$traits.coefs.median, family = "negative.binomial")
#' 
#' ## The easy way, using the same latent variables as the fitted model
#' simy <- simulate(object = fit_traits)
#' 
#' ## The easy way, generating new latent variables
#' simy <- simulate(object = fit_traits, new.lvs = TRUE)
#' }
#' 
#' 
#' ## Example 4 - simulate Bernoulli data, based on a model with two latent variables, 
#' ## no site variables, with two traits and one environmental covariates 
#' ## This example is a proof of concept that traits can used 
#' ## to explain environmental responses 
#' library(mvtnorm)
#' 
#' n <- 100; s <- 50
#' X <- as.matrix(scale(1:n))
#' colnames(X) <- c("elevation")
#' 
#' traits <- cbind(rbinom(s,1,0.5), rnorm(s)) 
#' ## one categorical and one continuous variable
#' colnames(traits) <- c("thorns-dummy","SLA")
#' 
#' simfit <- list(lv.coefs = cbind(rnorm(s), rmvnorm(s, mean = rep(0,2))), 
#'      lv.control = list(num.lv = 2, type = "independent"),
#'     traits.coefs = matrix(c(0.1,1,-0.5,1,0.5,0,-1,1), 2, byrow = TRUE))
#' rownames(simfit$traits.coefs) <- c("beta0","elevation")
#' colnames(simfit$traits.coefs) <- c("kappa0","thorns-dummy","SLA","sigma")
#' 
#' simy <- create.life(lv.control = simfit$lv.control, lv.coefs = simfit$lv.coefs, 
#'     X = X, traits = traits, traits.coefs = simfit$traits.coefs, 
#'     family = "binomial") 
#' 
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'     n.thin = 1)
#' 
#' example_which_traits <- vector("list",ncol(X)+1); 
#' for(i in 1:length(example_which_traits)) 
#'     example_which_traits[[i]] <- 1:ncol(traits)
#' fit_traits <- boral(y = simy, X = X, traits = traits, 
#'     which.traits = example_which_traits, family = "binomial", 
#'     lv.control = list(num.lv = 2), save.model = TRUE, 
#'     mcmc.control = example_mcmc_control)
#' }
#' 
#' 
NULL





#' Extract Deviance Information Criterion for a fitted model
#' 
#' Calculates the Deviance Information Criterion (DIC) for a model fitted using
#' JAGS.
#' 
#' Details regarding the Deviance Information Criterion may be found in
#' (Spiegelhalter et al., 2002; Ntzoufras, 2011; Gelman et al., 2013). The DIC
#' here is based on the conditional log-likelihood i.e., the latent variables
#' (and row effects if applicable) are treated as "fixed effects". A DIC based
#' on the marginal likelihood is obtainable from
#' \code{\link{get.more.measures}}, although this requires a much longer time
#' to compute. For models with overdispered count data, conditional DIC may not
#' perform as well as marginal DIC (Millar, 2009)
#' 
#' @name get.dic
#' @docType package
#' @param jagsfit The \code{jags.model} component of the output, from a model
#' fitted using \code{boral} with \code{save.model = TRUE}.
#' @return DIC value for the jags model.
#' @note This function and consequently the DIC value is automatically returned
#' when a model is fitted using \code{\link{boral}} with \code{calc.ics =
#' TRUE}.
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @references \itemize{ \item Gelman et al. (2013). Bayesian data analysis.
#' CRC press.
#' 
#' \item Millar, R. B. (2009). Comparison of hierarchical Bayesian models for
#' overdispersed count data using DIC and Bayes' factors. Biometrics, 65,
#' 962-969.
#' 
#' \item Ntzoufras, I. (2011). Bayesian modeling using WinBUGS (Vol. 698). John
#' Wiley & Sons.
#' 
#' \item Spiegelhalter, et al. (2002). Bayesian measures of model complexity
#' and fit. Journal of the Royal Statistical Society: Series B (Statistical
#' Methodology), 64, 583-639. }
#' @examples
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#'      
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' n <- nrow(y)
#' p <- ncol(y)
#'     
#' spiderfit_nb <- boral(y, family = "negative.binomial", lv.control = list(num.lv = 2),
#'      save.model = TRUE, calc.ics = TRUE, mcmc.control = example_mcmc_control)
#' 
#' spiderfit_nb$ics ## DIC returned as one of several information criteria.
#' }
#' 
NULL





#' Extract covariances and correlations due to shared environmental responses
#' 
#' Calculates the correlation between columns of the response matrix, due to
#' similarities in the response to explanatory variables i.e., shared
#' environmental response
#' 
#' In both independent response and correlated response models, where the each
#' of the columns of the response matrix \code{y} are fitted to a set of
#' explanatory variables given by \code{X}, the covariance and thus between two
#' columns \eqn{j} and \eqn{j'} due to similarities in their response to the
#' model matrix is calculated based on the linear predictors
#' \eqn{\bm{x}^\top_i\bm{\beta}_j} and \eqn{\bm{x}^\top_i\bm{\beta}_{j'})},
#' where \eqn{\bm{\beta}_j} are column-specific coefficients relating to the
#' explanatory variables.
#' 
#' For multivariate abundance data, the correlation calculated by this function
#' can be interpreted as the correlation attributable to similarities in the
#' environmental response between species. Such correlation matrices are
#' discussed and found in Ovaskainen et al., (2010), Pollock et al., 2014.
#' 
#' @name get.enviro.cor
#' @docType package
#' @param object An object for class "boral".
#' @param est A choice of either the posterior median (\code{est = "median"})
#' or posterior mean (\code{est = "mean"}), which are then treated as estimates
#' and the fitted values are calculated from. Default is posterior median.
#' @param prob A numeric scalar in the interval (0,1) giving the target
#' probability coverage of the intervals, by which to determine whether the
#' correlations are "significant". Defaults to 0.95.
#' @return A list with the following components: \item{cor, cor.lower,
#' cor.upper}{A set of \eqn{p \times p} correlation matrices, containing either
#' the posterior median or mean estimate plus lower and upper limits of the
#' corresponding 95\% highest posterior interval.}
#' 
#' \item{sig.cor}{A \eqn{p \times p} correlation matrix containing only the
#' ``significant" correlations whose 95\% highest posterior interval does not
#' contain zero. All non-significant correlations are set to zero.}
#' 
#' \item{cov}{A \eqn{p \times p} covariance matrix.}
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{get.residual.cor}}, which calculates the residual
#' correlation matrix for models involving latent variables.
#' @references \itemize{ \item Ovaskainen et al. (2010). Modeling species
#' co-occurrence by multivariate logistic regression generates new hypotheses
#' on fungal interactions. Ecology, 91, 2514-2521.
#' 
#' \item Pollock et al. (2014). Understanding co-occurrence by modelling
#' species simultaneously with a Joint Species Distribution Model (JSDM).
#' Methods in Ecology and Evolution, 5, 397-406. }
#' @examples
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#'      
#' library(mvabund) ## Load a dataset from the mvabund package
#' library(corrplot) ## For plotting correlations
#' data(spider)
#' y <- spider$abun
#' X <- scale(spider$x)
#' n <- nrow(y)
#' p <- ncol(y)
#'     
#' spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
#'      save.model = TRUE, mcmc.control = example_mcmc_control)
#' 
#' enviro.cors <- get.enviro.cor(spiderfit_nb)
#' 
#' corrplot(enviro.cors$sig.cor, title = "Shared response correlations", 
#' 	type = "lower", diag = FALSE, mar = c(3,0.5,2,1), tl.srt = 45)
#' }
#' 
NULL





#' Highest posterior density intervals for a fitted model
#' 
#' Calculates the lower and upper bounds of the highest posterior density
#' intervals for parameters and latent variables in a fitted model.
#' 
#' The function uses the \code{HPDinterval} function from the \code{coda}
#' package to obtain the HPD intervals. See \code{HPDinterval} for details
#' regarding the definition of the HPD interval.
#' 
#' @name get.hpdintervals
#' @docType package
#' @param y The response matrix that the model was fitted to.
#' @param X The model matrix used in the model. Defaults to \code{NULL}, in
#' which case it is assumed no model matrix was used.
#' @param traits The matrix of species traits used in the model. Defaults to
#' \code{NULL}, in which case it is assumed no traits were included.
#' @param row.ids A matrix with the number of rows equal to the number of rows
#' in \code{y}, and the number of columns equal to the number of row effects to
#' be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of
#' row \eqn{i} in \code{y} for random effect eqnj; please see
#' \code{\link{boral}} for details. Defaults to \code{NULL}, in which case iti
#' assumed no random effects were included in the model.
#' @param fit.mcmc All MCMC samples for the fitted model. These can be
#' extracted by fitting a model using \code{boral} with \code{save.model =
#' TRUE}, and then applying \code{get.mcmcsamples(fit)}.
#' @param lv.control A list (currently) with the following arguments: \itemize{
#' \item \emph{num.lv:} which specifies the number of true latent variables to
#' generate. Defaults to 0.
#' 
#' \item \emph{type:} which specifies the type the correlation structure of the
#' latent variables (across sites). Defaults to independence correlation
#' structure.
#' 
#' \item \emph{distmat:} which a distance matrix required to calculate
#' correlations across sites when a non-independence correlation structure on
#' the latent variables is imposed.  } Please see \code{\link{about.lvs}} for
#' more information.
#' @param prob A numeric scalar in the interval (0,1) giving the target
#' probability coverage of the intervals. Defaults to 0.95.
#' @param num.lv Old argument superceded by \code{lv.control}. Defaults to
#' \code{NULL} and ignored.
#' @return A list containing the following components where applicable:
#' 
#' \item{lv.coefs}{A three dimensional array giving the lower
#' \code{lv.coefs[,,"lower"]} and upper \code{lv.coefs[,,"upper"]} bounds of
#' the HPD intervals for the column-specific intercepts, latent variable
#' coefficients, and dispersion parameters if appropriate.}
#' 
#' \item{lv}{A three dimensional array giving the \code{lv.coefs[,,"lower"]}
#' and upper \code{lv.coefs[,,"upper"]} bounds of the HPD intervals for the
#' latent variables.}
#' 
#' \item{lv.covparams}{A matrix giving the lower and upper bounds of the HPD
#' intervals for the parameters characterizing the correlation structure of the
#' latent variables when they are assumed to be non-independent across rows.}
#' 
#' \item{row.coefs}{A list with each element being a matrix giving the lower
#' and upper bounds of the HPD intervals for row effects. The number of
#' elements in the list should equal the number of row effects included in the
#' model i.e., \code{ncol(row.ids)}.}
#' 
#' \item{row.sigma}{A list with each element being a vector giving the lower
#' and upper bounds of the HPD interval for the standard deviation of the
#' normal distribution for the row effects. The number of elements in the list
#' should equal the number of row effects included in the model i.e.,
#' \code{ncol(row.ids)}.}
#' 
#' \item{X.coefs}{A three dimensional array giving the lower
#' \code{lv.coefs[,,"lower"]} and upper \code{lv.coefs[,,"upper"]} bounds of
#' the HPD intervals for coefficients relating to \code{X}.}
#' 
#' \item{traits.coefs}{A three dimensional array giving the lower
#' \code{lv.coefs[,,"lower"]} and upper \code{lv.coefs[,,"upper"]} bounds of
#' the HPD intervals for coefficients and standard deviation relating to the
#' traits matrix \code{traits}.}
#' 
#' \item{cutoffs}{A matrix giving the lower and upper bounds of the HPD
#' intervals for common cutoffs in proportional odds regression.}
#' 
#' \item{powerparam}{A vector giving the lower and upper bounds of the HPD
#' interval for common power parameter in tweedie regression.}
#' @note \code{\link{boral}} fits the model and returns the HPD intervals by
#' default.
#' @section Warnings: \itemize{ \item HPD intervals tend to be quite wide, and
#' inference is somewhat tricky with them. This is made more difficult by the
#' multiple comparison problem due to the construction one interval for each
#' parameter!
#' 
#' \item Be very careful with interpretation of coefficients and HPD intervals
#' if different columns of \code{y} have different distributions!
#' 
#' \item HPD intervals for the cutoffs in proportional odds regression may be
#' poorly estimated for levels with few data. }
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @examples
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'     n.thin = 1)
#'      
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' n <- nrow(y)
#' p <- ncol(y)
#'     
#' ## Example 1 - model with two latent variables, site effects, 
#' ## 	and no environmental covariates
#' spiderfit_nb <- boral(y, family = "negative.binomial", lv.control = list(num.lv = 2),
#'     row.eff = "fixed", save.model = TRUE, mcmc.control = example_mcmc_control)
#' 
#' ## Returns a list with components corresponding to values described above.
#' spiderfit_nb$hpdintervals 
#' 
#' ## Example 2 - model with two latent variable, site effects, 
#' ## 	and environmental covariates
#' spiderfit_nb2 <- boral(y, X = spider$x, family = "negative.binomial", 
#'     lv.control = list(num.lv = 2), row.eff = "fixed", save.model = TRUE,
#'     mcmc.control = example_mcmc_control)
#' 
#' ## Returns a list with components corresponding to values described above.
#' spiderfit_nb2$hpdintervals 
#' 
#' }
#' 
NULL





#' Information Criteria for models
#' 
#' Calculates some information criteria for a fitted model, which could be used
#' for model selection. WARNING: As of version 1.6, this function will no
#' longer be updated...use at your own peril!!!
#' 
#' The following information criteria are currently calculated, when permitted:
#' 1) Widely Applicable Information Criterion (WAIC, Watanabe, 2010) based on
#' the conditional log-likelihood; 2) expected AIC (EAIC, Carlin and Louis,
#' 2011); 3) expected BIC (EBIC, Carlin and Louis, 2011); 4) AIC (using the
#' marginal likelihood) evaluated at the posterior median; 5) BIC (using the
#' marginal likelihood) evaluated at the posterior median.
#' 
#' 1) WAIC has been argued to be more natural and extension of AIC to the
#' Bayesian and hierarchical modeling context (Gelman et al., 2013), and is
#' based on the conditional log-likelihood calculated at each of the MCMC
#' samples.
#' 
#' 2 & 3) EAIC and EBIC were suggested by (Carlin and Louis, 2011). Both
#' criteria are of the form -2*mean(conditional log-likelihood) + penalty*(no.
#' of parameters in the model), where the mean is averaged all the MCMC
#' samples. EAIC applies a penalty of 2, while EBIC applies a penalty of
#' \eqn{log(n)}.
#' 
#' 4 & 5) AIC and BIC take the form -2*(marginal log-likelihood) + penalty*(no.
#' of parameters in the model), where the log-likelihood is evaluated at the
#' posterior median. If the parameter-wise posterior distributions are unimodal
#' and approximately symmetric, these will produce similar results to an AIC
#' and BIC where the log-likelihood is evaluated at the posterior mode. EAIC
#' applies a penalty of 2, while EBIC applies a penalty of \eqn{log(n)}.
#' 
#' Intuitively, comparing models with and without latent variables (using
#' information criteria such as those returned) amounts to testing whether the
#' columns of the response matrix \code{y} are correlated. With multivariate
#' abundance data for example, where \code{y} is a matrix of \eqn{n} sites and
#' \eqn{p} species, comparing models with and without latent variables tests
#' whether there is any evidence of correlation between species.
#' 
#' Please note that criteria 4 and 5 are not calculated all the time. In models
#' where traits are included in the model (such that the regression
#' coefficients \eqn{\beta_{0j}, \bm{\beta}_j} are random effects), or more
#' than two columns are ordinal responses (such that the intercepts
#' \eqn{\beta_{0j}} for these columns are random effects), then criteria 4 and
#' 5 are will not calculated. This is because the calculation of the marginal
#' log-likelihood in such cases currently fail to marginalize over such random
#' effects; please see the details in \code{calc.logLik.lv0} and
#' \code{calc.marglogLik}.
#' 
#' @name get.measures
#' @docType package
#' @param y The response matrix that the model was fitted to.
#' @param X The model matrix used in the model. Defaults to \code{NULL}, in
#' which case it is assumed no model matrix was used.
#' @param family Either a single element, or a vector of length equal to the
#' number of columns in \code{y}. The former assumes all columns of \code{y}
#' come from this distribution. The latter option allows for different
#' distributions for each column of \code{y}. Elements can be one of "binomial"
#' (with probit link), "poisson" (with log link), "negative.binomial" (with log
#' link), "normal" (with identity link), "lnormal" for lognormal (with log
#' link), "tweedie" (with log link), "exponential" (with log link), "gamma"
#' (with log link), "beta" (with logit link), "ordinal" (cumulative probit
#' regression).
#' 
#' Please see \code{\link{about.distributions}} for information on
#' distributions available in boral overall.
#' @param trial.size Either equal to a single element, or a vector of length
#' equal to the number of columns in y. If a single element, then all columns
#' assumed to be binomially distributed will have trial size set to this. If a
#' vector, different trial sizes are allowed in each column of y. The argument
#' is ignored for all columns not assumed to be binomially distributed.
#' Defaults to 1, i.e. Bernoulli distribution.
#' @param row.eff Single element indicating whether row effects are included as
#' fixed effects ("fixed"), random effects ("random") or not included ("none")
#' in the fitted model. If fixed effects, then for parameter identifiability
#' the first row effect is set to zero, which analogous to acting as a
#' reference level when dummy variables are used. If random effects, they are
#' drawn from a normal distribution with mean zero and estimated standard
#' deviation. Defaults to "none".
#' @param row.ids A matrix with the number of rows equal to the number of rows
#' in \code{y}, and the number of columns equal to the number of row effects to
#' be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of
#' row \eqn{i} in \code{y} for random effect eqnj; please see
#' \code{\link{boral}} for details. Defaults to \code{NULL}, so that if
#' \code{row.eff = "none"} then the argument is ignored, otherwise if \cr
#' \code{row.eff = "fixed"} or \code{"random"}, \cr then \code{row.ids =
#' matrix(1:nrow(y), ncol = 1)} i.e., a single, row effect unique to each row.
#' @param offset A matrix with the same dimensions as the response matrix
#' \code{y}, specifying an a-priori known component to be included in the
#' linear predictor during fitting. Defaults to \code{NULL}.
#' @param num.lv The number of latent variables used in the model.
#' @param fit.mcmc All MCMC samples for the fitted model. These can be
#' extracted by fitting a model using \code{\link{boral}} with \code{save.model
#' = TRUE}, and then applying \code{get.mcmcsamples(fit)}.
#' @return A list with the following components: \item{waic}{WAIC based on the
#' conditional log-likelihood.} \item{eaic}{EAIC based on the mean of the
#' conditional log-likelihood.} \item{ebic}{EBIC based on the mean of the
#' conditional log-likelihood.} \item{all.cond.logLik}{The conditional
#' log-likelihood evaluated at all MCMC samples. This is done via repeated
#' application of \code{\link{calc.condlogLik}}.} \item{cond.num.params}{Number
#' of estimated parameters used in the fitted model, when all parameters are
#' treated as "fixed" effects.} \item{do.marglik.ics}{A boolean indicating
#' whether marginal log-likelihood based information criteria are calculated.}
#' 
#' If \code{do.marglik.ics = TRUE}, then we also have: \item{median.logLik}{The
#' marginal log-likelihood evaluated at the posterior median.}
#' \item{marg.num.params}{Number of estimated parameters used in the fitted
#' model, when all parameters are treated as "fixed" effects.}
#' \item{aic.median}{AIC (using the marginal log-likelihood) evaluated at the
#' posterior median.} \item{bic.median}{BIC (using the marginal log-likelihood)
#' evaluated at the posterior median.}
#' @note When a model is fitted using \code{\link{boral}} with \code{calc.ics =
#' TRUE}, then this function is applied and the information criteria are
#' returned as part of the model output.
#' @section Warning: As of version 1.5, this function will no longer be
#' updated...use at your own peril!!!
#' 
#' Using information criterion for variable selection should be done with
#' extreme caution, for two reasons: 1) The implementation of these criteria
#' are both \emph{heuristic} and experimental. 2) Deciding what model to fit
#' for ordination purposes should be driven by the science. For example, it may
#' be the case that a criterion suggests a model with 3 or 4 latent variables.
#' However, if we interested in visualizing the data for ordination purposes,
#' then models with 1 or 2 latent variables are far more appropriate. As an
#' another example, whether or not we include row effects when ordinating
#' multivariate abundance data depends on if we are interested in differences
#' between sites in terms of relative species abundance (\code{row.eff =
#' FALSE}) or in terms of species composition (\code{row.eff = "fixed"}).
#' 
#' Also, the use of information criterion in the presence of variable selection
#' using SSVS is questionable.
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{get.dic}} for calculating the Deviance Information
#' Criterion (DIC) based on the conditional log-likelihood;
#' \code{\link{get.more.measures}} for even more information criteria.
#' @references \itemize{ \item Carlin, B. P., and Louis, T. A. (2011). Bayesian
#' methods for data analysis. CRC Press.
#' 
#' \item Gelman et al. (2013). Understanding predictive information criteria
#' for Bayesian models. Statistics and Computing, 1-20.
#' 
#' \item Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
#' and widely applicable information criterion in singular learning theory. The
#' Journal of Machine Learning Research, 11, 3571-3594. }
#' @examples
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#'      
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' n <- nrow(y)
#' p <- ncol(y)
#' 
#' spiderfit_pois <- boral(y, family = "poisson", 
#'     lv.control = list(num.lv = 2), row.eff = "random",
#'     mcmc.control = example_mcmc_control)
#' 
#' spiderfit_pois$ics ## Returns information criteria
#' 
#' spiderfit_nb <- boral(y, family = "negative.binomial", 
#'     lv.control = list(num.lv = 2), row.eff = "random",
#'     mcmc.control = example_mcmc_control)
#' 
#' spiderfit_nb$ics ## Returns the information criteria 
#' }
#' 
NULL





#' Additional Information Criteria for models
#' 
#' Calculates some information criteria beyond those from
#' \code{\link{get.measures}} for a fitted model, although this set of criteria
#' takes much longer to compute! WARNING: As of version 1.6, this function will
#' no longer be updated...use at your own peril!!!
#' 
#' Currently, four information criteria are calculated using this function,
#' when permitted: 1) AIC (using the marginal likelihood) evaluated at the
#' posterior mode; 2) BIC (using the marginal likelihood) evaluated at the
#' posterior mode; 3) Deviance information criterion (DIC) based on the
#' marginal log-likelihood; 4) Widely Applicable Information Criterion (WAIC,
#' Watanabe, 2010) based on the marginal log-likelihood. When uninformative
#' priors are used in fitting models, then the posterior mode should be
#' approximately equal to the maximum likelihood estimates.
#' 
#' All four criteria require computing the marginal log-likelihood across all
#' MCMC samples. This takes a very long time to run, since Monte Carlo
#' integration needs to be performed for all MCMC samples. Consequently, this
#' function is currently not implemented as an argument in main
#' \code{\link{boral}} fitting function, unlike \code{\link{get.measures}}
#' which is available via the \code{calc.ics = TRUE} argument.
#' 
#' Moreover, note these criteria are not calculated all the time. In models
#' where traits are included in the model (such that the regression
#' coefficients \eqn{\beta_{0j}, \bm{\beta}_j} are random effects), or more
#' than two columns are ordinal responses (such that the intercepts
#' \eqn{\beta_{0j}} for these columns are random effects), then these extra
#' information criteria are will not calculated, and the function returns
#' nothing except a simple message. This is because the calculation of the
#' marginal log-likelihood in such cases currently fail to marginalize over
#' such random effects; please see the details in \code{calc.logLik.lv0} and
#' \code{calc.marglogLik}.
#' 
#' The two main differences between the criteria and those returned from
#' \code{\link{get.measures}} are: \itemize{ \item The AIC and BIC computed
#' here are based on the log-likelihood evaluated at the posterior mode,
#' whereas the AIC and BIC from \code{\link{get.measures}} are evaluated at the
#' posterior median. The posterior mode and median will be quite close to one
#' another if the component-wise posterior distributions are unimodal and
#' symmetric. Furthermore, given uninformative priors are used, then both will
#' be approximate maximum likelihood estimators.
#' 
#' \item The DIC and WAIC computed here are based on the marginal
#' log-likelihood, whereas the DIC and WAIC from \code{\link{get.measures}} are
#' based on the conditional log-likelihood. Criteria based on the two types of
#' log-likelihood are equally valid, and to a certain extent, which one to use
#' depends on the question being answered i.e., whether to condition on the
#' latent variables or treat them as "random effects" (see discussions in
#' Spiegelhalter et al. 2002, and Vaida and Blanchard, 2005).  }
#' 
#' @name get.more.measures
#' @docType package
#' @param y The response matrix that the model was fitted to.
#' @param X The model matrix used in the model. Defaults to \code{NULL}, in
#' which case it is assumed no model matrix was used.
#' @param family Either a single element, or a vector of length equal to the
#' number of columns in \code{y}. The former assumes all columns of \code{y}
#' come from this distribution. The latter option allows for different
#' distributions for each column of \code{y}. Elements can be one of "binomial"
#' (with probit link), "poisson" (with log link), "negative.binomial" (with log
#' link), "normal" (with identity link), "lnormal" for lognormal (with log
#' link), "tweedie" (with log link), "exponential" (with log link), "gamma"
#' (with log link), "beta" (with logit link), "ordinal" (cumulative probit
#' regression).
#' 
#' Please see \code{\link{about.distributions}} for information on
#' distributions available in boral overall.
#' @param trial.size Either equal to a single element, or a vector of length
#' equal to the number of columns in y. If a single element, then all columns
#' assumed to be binomially distributed will have trial size set to this. If a
#' vector, different trial sizes are allowed in each column of y. The argument
#' is ignored for all columns not assumed to be binomially distributed.
#' Defaults to 1, i.e. Bernoulli distribution.
#' @param row.eff Single element indicating whether row effects are included as
#' fixed effects ("fixed"), random effects ("random") or not included ("none")
#' in the fitted model. If fixed effects, then for parameter identifiability
#' the first row effect is set to zero, which analogous to acting as a
#' reference level when dummy variables are used. If random effects, they are
#' drawn from a normal distribution with mean zero and estimated standard
#' deviation. Defaults to "none".
#' @param row.ids A matrix with the number of rows equal to the number of rows
#' in \code{y}, and the number of columns equal to the number of row effects to
#' be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of
#' row \eqn{i} in \code{y} for random effect eqnj; please see
#' \code{\link{boral}} for details. Defaults to \code{NULL}, so that if
#' \code{row.eff = "none"} then the argument is ignored, otherwise if \cr
#' \code{row.eff = "fixed"} or \code{"random"}, \cr then \code{row.ids =
#' matrix(1:nrow(y), ncol = 1)} i.e., a single, row effect unique to each row.
#' @param offset A matrix with the same dimensions as the response matrix
#' \code{y}, specifying an a-priori known component to be included in the
#' linear predictor during fitting. Defaults to \code{NULL}.
#' @param num.lv The number of latent variables used in the fitted model.
#' @param fit.mcmc All MCMC samples for the fitted model. These can be
#' extracted by fitting a model using \code{\link{boral}} with \code{save.model
#' = TRUE}, and then applying \code{get.mcmcsamples(fit)}.
#' @param verbose If TRUE, a notice is printed every 100 samples indicating
#' progress in calculation of the marginal log-likelihood. Defaults to
#' \code{TRUE}.
#' @return If calculated, then a list with the following components:
#' \item{marg.aic}{AIC (using on the marginal log-likelihood) evaluated at
#' posterior mode.} \item{marg.bic}{BIC (using on the marginal log-likelihood)
#' evaluated at posterior mode.} \item{marg.dic}{DIC based on the marginal
#' log-likelihood.} \item{marg.waic}{WAIC based on the marginal
#' log-likelihood.} \item{all.marg.logLik}{The marginal log-likelihood
#' evaluated at all MCMC samples. This is done via repeated application of
#' \code{\link{calc.marglogLik}}.} \item{num.params}{Number of estimated
#' parameters used in the fitted model.}
#' @section Warning: As of version 1.6, this function will no longer be
#' updated...use at your own peril!!!
#' 
#' Using information criterion for variable selection should be done with
#' extreme caution, for two reasons: 1) The implementation of these criteria
#' are both \emph{heuristic} and experimental. 2) Deciding what model to fit
#' for ordination purposes should be driven by the science. For example, it may
#' be the case that a criterion suggests a model with 3 or 4 latent variables.
#' However, if we interested in visualizing the data for ordination purposes,
#' then models with 1 or 2 latent variables are far more appropriate. As an
#' another example, whether or not we include row effects when ordinating
#' multivariate abundance data depends on if we are interested in differences
#' between sites in terms of relative species abundance (\code{row.eff =
#' FALSE}) or in terms of species composition (\code{row.eff = "fixed"}).
#' 
#' Also, the use of information criterion in the presence of variable selection
#' using SSVS is questionable.
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{get.measures}} for several information criteria which
#' take considerably less time to compute, and are automatically implemented in
#' \code{\link{boral}} with \code{calc.ics = TRUE}.
#' @references \itemize{ \item Spiegelhalter et al. (2002). Bayesian measures
#' of model complexity and fit. Journal of the Royal Statistical Society:
#' Series B (Statistical Methodology), 64, 583-639.
#' 
#' \item Vaida, F., and Blanchard, S. (2005). Conditional Akaike information
#' for mixed-effects models. Biometrika, 92, 351-370.
#' 
#' \item Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
#' and widely applicable information criterion in singular learning theory. The
#' Journal of Machine Learning Research, 11, 3571-3594. }
#' @examples
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'     n.thin = 1)
#'      
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' n <- nrow(y)
#' p <- ncol(y)
#' 
#' spiderfit_nb <- boral(y, family = "negative.binomial", lv.control = list(num.lv = 2),
#'     row.eff = "fixed", save.model = TRUE, calc.ics = TRUE,
#'     mcmc.control = example_mcmc_control)
#' 
#' ## Extract MCMC samples
#' fit_mcmc <- get.mcmcsamples(spiderfit_nb)
#' 
#' ## WATCH OUT! The following takes a very long time to run!
#' get.more.measures(y, family = "negative.binomial", 
#'     num.lv = spiderfit_nb$num.lv, fit.mcmc = fit_mcmc, 
#'     row.eff = "fixed", row.ids = spiderfit_nb$row.ids)		
#' 
#'      
#' ## Illustrating what happens in a case where these criteria will 
#' ## 	not be calculated.
#' data(antTraits)
#' y <- antTraits$abun
#' X <- as.matrix(scale(antTraits$env))
#' ## Include only traits 1, 2, and 5
#' traits <- as.matrix(antTraits$traits[,c(1,2,5)])
#' example_which_traits <- vector("list",ncol(X)+1)
#' for(i in 1:length(example_which_traits))
#'     example_which_traits[[i]] <- 1:ncol(traits)
#' 
#' fit_traits <- boral(y, X = X, traits = traits, lv.control = list(num.lv = 2), 
#'     which.traits = example_which_traits, family = "negative.binomial", 
#'     save.model = TRUE, mcmc.control = example_mcmc_control)
#'      
#' ## Extract MCMC samples
#' fit_mcmc <- get.mcmcsamples(fit_traits)
#' 
#' get.more.measures(y, X = X, family = "negative.binomial", 
#'     num.lv = fit.traits$num.lv, fit.mcmc = fit_mcmc)	
#' }
#' 
NULL





#' Extract residual correlations and precisions from models
#' 
#' Calculates the residual correlation and precision matrices from models that
#' include latent variables.
#' 
#' In models with latent variables, the residual covariance matrix is
#' calculated based on the matrix of latent variables regression coefficients
#' formed by stacking the rows of \eqn{\bm{\theta}_j}. That is, if we denote
#' \eqn{\bm{\Theta} = (\bm{\theta}_1 \ldots \bm{\theta}_p)'}, then the residual
#' covariance and hence residual correlation matrix is calculated based on
#' \eqn{\bm{\Theta}\bm{\Theta}'}.
#' 
#' For multivariate abundance data, the inclusion of latent variables provides
#' a parsimonious method of accounting for correlation between species.
#' Specifically, the linear predictor,
#' 
#' \deqn{\beta_{0j} + \bm{x}^\top_i\bm{\beta}_j + \bm{z}^\top_i\bm{\theta}_j}
#' 
#' is normally distributed with a residual covariance matrix given by
#' \eqn{\bm{\Theta}\bm{\Theta}'}. A strong residual covariance/correlation
#' matrix between two species can then be interpreted as evidence of species
#' interaction (e.g., facilitation or competition), missing covariates, as well
#' as any additional species correlation not accounted for by shared
#' environmental responses (see also Pollock et al., 2014, for residual
#' correlation matrices in the context of Joint Species Distribution Models).
#' 
#' The residual precision matrix (also known as partial correlation matrix,
#' Ovaskainen et al., 2016) is defined as the inverse of the residual
#' correlation matrix. The precision matrix is often used to identify direct or
#' causal relationships between two species e.g., two species can have a zero
#' precision but still be correlated, which can be interpreted as saying that
#' two species do not directly interact, but they are still correlated through
#' other species. In other words, they are conditionally independent given the
#' other species. It is important that the precision matrix does not exhibit
#' the exact same properties of the correlation e.g., the diagonal elements are
#' not equal to 1. Nevertheless, relatively larger values of precision imply a
#' stronger direct relationships between two species.
#' 
#' In addition to the residual correlation and precision matrices, the median
#' or mean point estimator of trace of the residual covariance matrix is
#' returned, \eqn{\sum\limits_{j=1}^p [\bm{\Theta}\bm{\Theta}']_{jj}}. Often
#' used in other areas of multivariate statistics, the trace may be interpreted
#' as the amount of covariation explained by the latent variables. One
#' situation where the trace may be useful is when comparing a pure LVM versus
#' a model with latent variables and some predictors (correlated response
#' models) -- the proportional difference in trace between these two models may
#' be interpreted as the proportion of covariation between species explained by
#' the predictors. Of course, the trace itself is random due to the MCMC
#' sampling, and so it is not always guranteed to produce sensible answers!
#' 
#' @name get.residual.cor
#' @docType package
#' @param object An object for class "boral".
#' @param est A choice of either the posterior median (\code{est = "median"})
#' or posterior mean (\code{est = "mean"}), which are then treated as estimates
#' and the fitted values are calculated from. Default is posterior median.
#' @param prob A numeric scalar in the interval (0,1) giving the target
#' probability coverage of the intervals, by which to determine whether the
#' correlations and precisions are "significant". Defaults to 0.95.
#' @return A list with the following components: \item{cor, cor.lower,
#' cor.upper}{A set of \eqn{p \times p} correlation matrices, containing either
#' the posterior median or mean estimate plus lower and upper limits of the
#' corresponding 95\% highest posterior interval.}
#' 
#' \item{sig.cor}{A \eqn{p \times p} correlation matrix containing only the
#' ``significant" correlations whose 95\% highest posterior interval does not
#' contain zero. All non-significant correlations are set to zero.}
#' 
#' \item{cov}{A \eqn{p \times p} covariance matrix.}
#' 
#' \item{prec, prec.lower, prec.upper}{A set of \eqn{p \times p} precision
#' matrices, containing either the posterior median or mean estimate plus lower
#' and upper limits of the corresponding 95\% highest posterior interval.}
#' 
#' \item{sig.prec}{A \eqn{p \times p} residual precision matrix containing only
#' the ``significant" precisions whose 95\% highest posterior interval does not
#' contain zero. All non-significant precision are set to zero.}
#' 
#' \item{trace}{The median/mean point estimator of the trace (sum of the
#' diagonal elements) of the residual covariance matrix.}
#' @note Residual correlation and precision matrices are reliably modeled only
#' with two or more latent variables i.e., \code{num.lv > 1} when fitting the
#' model using \code{boral}.
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{get.enviro.cor}}, which calculates the correlation
#' matrix due to similarities in the response to the explanatory variables
#' (i.e., similarities due to a shared environmental response).
#' @references \itemize{ \item Ovaskainen et al. (2016). Using latent variable
#' models to identify large networks of species-to-species associations at
#' different spatial scales. Methods in Ecology and Evolution, 7, 549-555.
#' 
#' \item Pollock et al. (2014). Understanding co-occurrence by modelling
#' species simultaneously with a Joint Species Distribution Model (JSDM).
#' Methods in Ecology and Evolution, 5, 397-406.
#' 
#' }
#' @examples
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#' 
#' library(mvabund) ## Load a dataset from the mvabund package
#' library(corrplot) ## For plotting correlations
#' data(spider)
#' y <- spider$abun
#' n <- nrow(y)
#' p <- ncol(y)
#'     
#' spiderfit_nb <- boral(y, X = spider$x, family = "negative.binomial", 
#'     lv.control = list(num.lv = 2), save.model = TRUE, 
#'     mcmc.control = example_mcmc_control)
#' 
#' res.cors <- get.residual.cor(spiderfit_nb)
#' 
#' corrplot(res.cors$sig.cor, title = "Residual correlations", 
#'     type = "lower", diag = FALSE, mar = c(3,0.5,2,1), tl.srt = 45)
#' }
#' 
NULL





#' Plot the latent variables from a fitted model
#' 
#' Construct a 1-D index plot or 2-D scatterplot of the latent variables, and
#' their corresponding coefficients i.e., a biplot, from a fitted model.
#' 
#' This function allows an ordination plot to be constructed, based on either
#' the posterior medians and posterior means of the latent variables
#' respectively depending on the choice of \code{est}. The latent variables are
#' labeled using the row index of the response matrix \code{y}. If the fitted
#' model contains more than two latent variables, then one can specify which
#' latent variables i.e., ordination axes, to plot based on the
#' \code{which.lvs} argument. This can prove useful (to check) if certain sites
#' are outliers on one particular ordination axes.
#' 
#' If the fitted model did not contain any covariates, the ordination plot can
#' be interpreted in the exactly same manner as unconstrained ordination plots
#' constructed from methods such as Nonmetric Multi-dimensional Scaling (NMDS,
#' Kruskal, 1964) and Correspondence Analysis (CA, Hill, 1974). With
#' multivariate abundance data for instance, where the response matrix \code{y}
#' consists of \eqn{n} sites and \eqn{p} species, the ordination plots can be
#' studied to look for possible clustering of sites, location and/or dispersion
#' effects, an arch pattern indicative of some sort species succession over an
#' environmental gradient, and so on.
#' 
#' If the fitted model did include covariates, then a ``residual ordination"
#' plot is produced, which can be interpreted can offering a graphical
#' representation of the (main patterns of) residual covarations, i.e.
#' covariations after accounting for the covariates. With multivariate
#' abundance data for instance, these residual ordination plots represent could
#' represent residual species co-occurrence due to phylogency, species
#' competition and facilitation, missing covariates, and so on (Warton et al.,
#' 2015)
#' 
#' If \code{biplot = TRUE}, then a biplot is constructed so that both the
#' latent variables and their corresponding coefficients are included in their
#' plot (Gabriel, 1971). The latent variable coefficients are shown in red, and
#' are indexed by the column names of \code{y}. The number of latent variable
#' coefficients to plot is controlled by \code{ind.spp}. In ecology for
#' example, often we are only be interested in the "indicator" species, e.g.
#' the species with most represent a particular set of sites or species with
#' the strongest covariation (see Chapter 9, Legendre and Legendre, 2012, for
#' additional discussion). In such case, we can then biplot only the
#' \code{ind.spp} "most important" species, as indicated by the the L2-norm of
#' their latent variable coefficients.
#' 
#' As with correspondence analysis, the relative scaling of the latent
#' variables and the coefficients in a biplot is essentially arbitrary, and
#' could be adjusted to focus on the sites, species, or put even weight on both
#' (see Section 9.4, Legendre and Legendre, 2012). In \code{lvsplot}, this
#' relative scaling is controlled by the \code{alpha} argument, which basically
#' works by taking the latent variables to a power \code{alpha} and the latent
#' variable coefficients to a power \code{1-alpha}.
#' 
#' For latent variable models, we are generally interested in "symmetric plots"
#' that place the latent variables and their coefficients on the same scale. In
#' principle, this is achieved by setting \code{alpha = 0.5}, the default
#' value, although sometimes this needs to be tweaked slighlty to a value
#' between 0.45 and 0.55 (see also the \code{corresp} function in the
#' \code{MASS} package that also produces symmetric plots, as well as Section
#' 5.4, Borcard et al., 2011 for more details on scaling).
#' 
#' @name lvsplot
#' @docType package
#' @param x An object for class "boral".
#' @param jitter If \code{jitter = TRUE}, then some jittering is applied so
#' that points on the plots do not overlap exactly (which can often occur with
#' discrete data, small sample sizes, and if some sites are identical in terms
#' species co-occurence). Please see \code{\link{jitter}} for its
#' implementation. Defaults to \code{FALSE}.
#' @param biplot If \code{biplot = TRUE}, then a biplot is construct such that
#' both the latent variables \emph{and} their corresponding coefficients are
#' plotted. Otherwise, only the latent variable scores are plotted. Defaults to
#' \code{TRUE}.
#' @param ind.spp Controls the number of latent variable coefficients to plot
#' if \code{biplot = TRUE}. If \code{ind.spp} is an integer, then only the
#' first \code{ind.spp} "most important" latent variable coefficients are
#' included in the biplot, where "most important" means the latent variable
#' coefficients with the largests L2-norms. Defaults to \code{NULL}, in which
#' case all latent variable coefficients are included in the biplot.
#' @param alpha A numeric scalar between 0 and 1 that is used to control the
#' relative scaling of the latent variables and their coefficients, when
#' constructing a biplot. Defaults to 0.5, and we typically recommend between
#' 0.45 to 0.55 so that the latent variables and their coefficients are on
#' roughly the same scale.
#' @param main Title for resulting ordination plot. Defaults to \code{NULL}, in
#' which case a "standard" title is used.
#' @param est A choice of either the posterior median (\code{est = "median"})
#' or posterior mean (\code{est = "mean"}), which are then treated as estimates
#' and the ordinations based off. Default is posterior median.
#' @param which.lvs A vector of length two, indicating which latent variables
#' (ordination axes) to plot which \code{x} is an object with two or more
#' latent variables. The argument is ignored is \code{x} only contains one
#' latent variables. Defaults to \code{which.lvs = c(1,2)}.
#' @param return.vals If \code{TRUE}, then the \emph{scaled} latent variables
#' scores and corresponding scaled coefficients are returned (based on the
#' value of \code{alpha} used). This is useful if the user wants to construct
#' their own custom model-based ordinations. Defaults to \code{FALSE}.
#' @param ... Additional graphical options to be included in. These include
#' values for \cr \code{cex, cex.lab, cex.axis, cex.main, lwd}, and so on.
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @references \itemize{ \item Borcard et al. (2011). Numerical Ecology with R.
#' Springer.
#' 
#' \item Gabriel, K. R. (1971). The biplot graphic display of matrices with
#' application to principal component analysis. Biometrika, 58, 453-467.
#' 
#' \item Hill, M. O. (1974). Correspondence analysis: a neglected multivariate
#' method. Applied statistics, 23, 340-354.
#' 
#' \item Kruskal, J. B. (1964). Nonmetric multidimensional scaling: a numerical
#' method. Psychometrika, 29, 115-129.
#' 
#' \item Legendre, P. and Legendre, L. (2012). Numerical ecology, Volume 20.
#' Elsevier.
#' 
#' \item Warton et al. (2015). So Many Variables: Joint Modeling in Community
#' Ecology. Trends in Ecology and Evolution, to appear }
#' @examples
#' 
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#' 
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' n <- nrow(y)
#' p <- ncol(y)
#'      
#' spiderfit_nb <- boral(y, family = "negative.binomial", lv.control = list(num.lv = 2),
#'     row.eff = "fixed", mcmc.control = example_mcmc_control)
#' 
#' lvsplot(spiderfit_nb) 
#' 
NULL





#' Write a text file containing a model for use into JAGS
#' 
#' This function is designed to write models with one or more latent variables.
#' 
#' This function is automatically executed inside \code{\link{boral}}, and
#' therefore does not need to be run separately before fitting the model. It
#' can however be run independently if one is: 1) interested in what the actual
#' JAGS file for a particular model looks like, 2) wanting to modify a basic
#' JAGS model file to construct more complex model e.g., include environmental
#' variables.
#' 
#' Please note that \code{\link{boral}} currently does not allow the user to
#' manually enter a script to be run.
#' 
#' When running the main function \code{\link{boral}}, setting \code{save.model
#' = TRUE} which automatically save the JAGS model file as a text file (with
#' name based on the \code{model.name}) in the current working directory.
#' 
#' @name make.jagsboralmodel
#' @docType package
#' @param family Either a single element, or a vector of length equal to the
#' number of columns in \code{y}. The former assumes all columns of \code{y}
#' come from this distribution. The latter option allows for different
#' distributions for each column of \code{y}. Elements can be one of "binomial"
#' (with probit link), "poisson" (with log link), "negative.binomial" (with log
#' link), "normal" (with identity link), "lnormal" for lognormal (with log
#' link), "tweedie" (with log link), "exponential" (with log link), "gamma"
#' (with log link), "beta" (with logit link), "ordinal" (cumulative probit
#' regression).
#' 
#' Please see \code{\link{about.distributions}} for information on
#' distributions available in boral overall.
#' @param num.X Number of columns in \code{X}. Defaults to 0, in which case it
#' is assumed that no covariates are included in the model. Recall that no
#' intercept should be included in \code{X}.
#' @param X.ind An matrix of 1s and 0s, indicating whether a particular
#' covariate should be included (1) or excluded (0) in the mean structure of a
#' particular response. The matrix should the number of rows equal to the
#' number of columns in \code{y}, and the number of columns equal to the number
#' of columns in \code{X}. Defaults to \code{NULL}, in which case it is assumed
#' that all covariates are included in the mean structure of all responses
#' i.e., all 1s.
#' @param num.traits Number of columns in the model matrix \code{traits}.
#' Defaults to 0, in which case it is assumed no traits are included in model.
#' Recall that no intercept should be included in \code{traits}.
#' @param which.traits A list of length equal to (number of columns in \code{X}
#' + 1), informing which columns of \code{traits} the column-specific
#' intercepts and each of the column-specific regression coefficients should be
#' regressed against. The first element in the list applies to the
#' column-specific intercept, while the remaining elements apply to the
#' regression coefficients. Each element of \code{which.traits} is a vector
#' indicating which traits are to be used.
#' 
#' For example, if \code{which.traits[[2]] = c(2,3)}, then the regression
#' coefficients corresponding to the first column in \code{X} are regressed
#' against the second and third columns of \code{traits}. If
#' \code{which.traits[[2]][1] = 0}, then the regression coefficients for each
#' column are treated as independent. Please see \code{\link{about.traits}} for
#' more details.
#' 
#' Defaults to \code{NULL}, and used in conjunction with \code{traits} and \cr
#' \code{prior.control$ssvs.traitsindex}.
#' @param lv.control A list (currently) with the following arguments: \itemize{
#' \item \emph{num.lv:} which specifies the number of true latent variables to
#' generate. Defaults to 0.
#' 
#' \item \emph{type:} which specifies the type the correlation structure of the
#' latent variables (across sites). Defaults to independence correlation
#' structure.  } Please see \code{\link{about.lvs}} for more information.
#' @param row.eff Single element indicating whether row effects are included as
#' fixed effects ("fixed"), random effects ("random") or not included ("none")
#' in the fitted model. If fixed effects, then for parameter identifiability
#' the first row effect is set to zero, which analogous to acting as a
#' reference level when dummy variables are used. If random effects, they are
#' drawn from a normal distribution with mean zero and unknown standard
#' deviation. Defaults to "none".
#' @param row.ids A matrix with the number of rows equal to the number of rows
#' in \code{y}, and the number of columns equal to the number of row effects to
#' be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of
#' row \eqn{i} in \code{y} for random effect eqnj; please see
#' \code{\link{boral}} for details. Defaults to \code{NULL}, so that if
#' \code{row.eff = "none"} then the argument is ignored, otherwise if \cr
#' \code{row.eff = "fixed"} or \code{"random"}, \cr then \code{row.ids =
#' matrix(1:nrow(y), ncol = 1)} i.e., a single, row effect unique to each row.
#' @param offset A matrix with the same dimensions as the response matrix
#' \code{y}, specifying an a-priori known component to be included in the
#' linear predictor during fitting. Defaults to \code{NULL}.
#' @param trial.size Either equal to a single element, or a vector of length
#' equal to the number of columns in y. If a single element, then all columns
#' assumed to be binomially distributed will have trial size set to this. If a
#' vector, different trial sizes are allowed in each column of y. The argument
#' is ignored for all columns not assumed to be binomially distributed.
#' Defaults to 1, i.e. Bernoulli distribution.
#' @param n The number of rows in the response matrix \code{y}.
#' @param p The number of columns in the response matrix \code{y}.
#' @param model.name Name of the text file that the JAGS script is written to.
#' Defaults to \code{NULL}, in which case the default of "jagsboralmodel.txt"
#' is used.
#' @param prior.control A list of parameters for controlling the prior
#' distributions. These include: \itemize{ \item \emph{type:} Vector of four
#' strings indicating the type of prior distributions to use. In order, these
#' are: 1) priors for all column-specific intercepts, row effects, and cutoff
#' points for ordinal data; 2) priors for the latent variable coefficients and
#' covariance parameters. This is ignored if \code{lv.control$num.lv = 0}; 3)
#' priors for all column-specific coefficients relating to \code{X} (ignored if
#' \code{X = NULL}). When traits are included in the model, this is also the
#' prior for the trait regression coefficients (please see
#' \code{\link{about.traits}} for more information); 4) priors for any
#' dispersion parameters and variance (standard deviation, to be precise)
#' parameters in the model.
#' 
#' For elements 1-3, the prior distributions currently available include: I)
#' ``normal", which is a normal prior with the variance controlled by elements
#' 1-3 in \code{hypparams}; II) ``cauchy", which is a Cauchy prior with
#' variance controlled by elements 1-3 in \code{hypparams}. Gelman, et al.
#' (2008) considers using Cauchy priors with variance \eqn{2.5^2} as weakly
#' informative priors for coefficients in logistic and potentially other
#' generalized linear models; III) ``uniform", which is a symmetric uniform
#' prior with minimum and maximum values controlled by element 1-3 in
#' \code{hypparams}.
#' 
#' For element 4, the prior distributions currently available include: I)
#' ``uniform", which is uniform prior with minimum zero and maximum controlled
#' by element 4 in \code{hypparmas}; II) ``halfnormal", which is half-normal
#' prior with variance controlled by \code{hypparams}; III) ``halfcauchy",
#' which is a half-Cauchy prior with variance controlled by element 4 in
#' \code{hypparams}.
#' 
#' Defaults to the vector \code{c("normal","normal","normal","uniform")}.
#' 
#' \item \emph{hypparams:} Vector of four hyperparameters used in the set up of
#' prior distributions. In order, these are: 1) affects the prior distribution
#' for all column-specific intercepts, row effects, and cutoff points for
#' ordinal data; 2) affects the prior distribution for all latent variable
#' coefficients and correlation parameters. This is ignored if
#' \code{lv.control$num.lv = 0}; 3) affects the prior distribution for
#' column-specific coefficients relating to \code{X} (ignored if \code{X =
#' NULL}). When traits are included in the model, it also affects the prior
#' distribution for the trait regression coefficients; 4) affects the prior
#' distribution for any dispersion parameters, as well as the prior
#' distributions for the standard deviation of the random effects normal
#' distribution if \code{row.eff = "random"}, the standard deviation of the
#' column-specific random intercepts for these columns if more than two of the
#' columns are ordinal, and the standard deviation of the random effects normal
#' distribution for trait regression coefficients when traits are included in
#' the model.
#' 
#' Defaults to the vector \code{c(10, 10, 10, 30)}. The use of normal
#' distributions with mean zero and variance 10 as priors is seen as one type
#' of (very) weakly informative prior, according to
#' \href{https://github.com/stan-dev/stan/wiki/Prior-Choice-RecommendationsPrior
#' choice recommendations}.
#' 
#' \item \emph{ssvs.index:} Indices to be used for stochastic search variable
#' selection (SSVS, George and McCulloch, 1993). Either a single element or a
#' vector with length equal to the number of columns in the implied model
#' matrix \code{X}. Each element can take values of -1 (no SSVS is performed on
#' this covariate), 0 (SSVS is performed on individual coefficients for this
#' covariate), or any integer greater than 0 (SSVS is performed on collectively
#' all coefficients on this covariate/s.)
#' 
#' Please see \code{\link{about.ssvs}} for more information regarding the
#' implementation of SSVS. Defaults to -1, in which case SSVS is not performed
#' on \code{X} variables.
#' 
#' \item \emph{ssvs.g:} Multiplicative, shrinkage factor for SSVS, which
#' controls the strength of the "spike" in the SSVS mixture prior. In summary,
#' if the coefficient is included in the model, the "slab" prior is a normal
#' distribution with mean zero and variance given by element 3 in
#' \code{hypparams}, while if the coefficient is not included in the model, the
#' "spike" prior is normal distribution with mean zero and variance given by
#' element 3 in \code{hypparams} multiplied by \code{ssvs.g}. Please see
#' \code{\link{about.ssvs}} for more information regarding the implementation
#' of SSVS. Defaults to 1e-6.
#' 
#' \item \emph{ssvs.traitsindex:} Used in conjunction with \code{traits} and
#' \code{which.traits}, this is a list of indices to be used for performing
#' SSVS on the trait coefficients. Should be a list with the same length as
#' \code{which.traits}, and with each element a vector of indices with the same
#' length as the corresponding element in \code{which.traits}. Each index
#' either can take values of -1 (no SSVS on this trait coefficient) or 0 (no
#' SSVS on this trait coefficient).
#' 
#' Please see \code{\link{about.ssvs}} for more information regarding the
#' implementation of SSVS. Defaults to -1, in which case SSVS is not performed
#' on any of the trait coefficients, if they are included in the model.  }
#' @param num.lv Old argument superceded by \code{lv.control}. Defaults to
#' \code{NULL} and ignored.
#' @return A text file is created, containing the model to be called by the
#' boral function for entering into JAGS. This file is automatically deleted
#' once boral has finished running \code{save.model = TRUE}.
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{make.jagsboralnullmodel}} for writing JAGS scripts for
#' models with no latent variables i.e., so-called "null models".
#' @references \itemize{ \item Gelman, et al. (2008). A weakly informative
#' default prior distribution for logistic and other regression models. The
#' Annals of Applied Statistics, 2, 1360-1383. }
#' @examples
#' 
#' library(mvtnorm)
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' n <- nrow(y)
#' p <- ncol(y)
#' 
#' ## Example 1 - Create a JAGS model file, where distributions alternative 
#' ## between Poisson and negative binomial distributions 
#' ##   across the rows of y.
#' make.jagsboralmodel(family = rep(c("poisson","negative.binomial"),length=p), 
#'     row.eff = "fixed", num.X = 0, n = n, p = p)
#' 
#' ## Example 2 - Create a JAGS model file, where distributions are all 
#' ##	negative binomial distributions and covariates will be included.
#' make.jagsboralmodel(family = "negative.binomial", num.X = ncol(spider$x),
#'     n = n, p = p)
#' 
#' 	
#' ## Example 3 - Simulate some ordinal data and create a JAGS model file
#' ## 30 rows (sites) with two latent variables 
#' true.lv <- rbind(rmvnorm(15,mean=c(-2,-2)),rmvnorm(15,mean=c(2,2)))
#' ## 10 columns (species)
#' true.lv.coefs <- rmvnorm(10,mean = rep(0,3)); 
#' true.lv.coefs[nrow(true.lv.coefs),1] <- -sum(true.lv.coefs[-nrow(true.lv.coefs),1])
#' ## Impose a sum-to-zero constraint on the column effects
#' true.ordinal.cutoffs <- seq(-2,10,length=10-1)
#' 
#' simy <- create.life(true.lv = true.lv, lv.coefs = true.lv.coefs, 
#'     family = "ordinal", cutoffs = true.ordinal.cutoffs) 
#' 
#' make.jagsboralmodel(family = "ordinal", num.X = 0, 
#'     row.eff = FALSE, n=30, p=10, model.name = "myawesomeordmodel.txt")
#' 
#' 
#' ## Have a look at the JAGS model file for a model involving traits,
#' ## based on the ants data from mvabund.
#' library(mvabund)
#' data(antTraits)
#' 
#' y <- antTraits$abun
#' X <- as.matrix(antTraits$env)
#' ## Include only traits 1, 2, and 5, plus an intercept
#' traits <- as.matrix(antTraits$traits[,c(1,2,5)])
#' ## Please see help file for boral regarding the use of which.traits
#' example_which_traits <- vector("list",ncol(X)+1)
#' for(i in 1:length(example_which_traits)) 
#'     example_which_traits[[i]] <- 1:ncol(traits)
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'     n.thin = 1)
#' 
#' fit_traits <- boral(y, X = X, traits = traits, which.traits = example_which_traits, 
#'     family = "negative.binomial", lv.control = list(num.lv = 2),
#'     model.name = "anttraits.txt", mcmc.control = example_mcmc_control)
#' }
#' 
#' 
NULL





#' Write a text file containing a model for use into JAGS
#' 
#' This function is designed to write models with no latent variables i.e.,
#' so-called "null" models.
#' 
#' This function is automatically executed inside \code{\link{boral}}, and
#' therefore does not need to be run separately before fitting the model. It
#' can however be run independently if one is: 1) interested in what the actual
#' JAGS file for a particular model looks like, 2) wanting to modify a basic
#' JAGS model file to construct more complex model e.g., include environmental
#' variables.
#' 
#' Please note that \code{\link{boral}} currently does not allow the user to
#' manually enter a script to be run.
#' 
#' When running the main function \code{\link{boral}}, setting \code{save.model
#' = TRUE} which automatically save the JAGS model file as a text file (with
#' name based on the \code{model.name}) in the current working directory.
#' 
#' @name make.jagsboralnullmodel
#' @docType package
#' @param family Either a single element, or a vector of length equal to the
#' number of columns in \code{y}. The former assumes all columns of \code{y}
#' come from this distribution. The latter option allows for different
#' distributions for each column of \code{y}. Elements can be one of "binomial"
#' (with probit link), "poisson" (with log link), "negative.binomial" (with
#' lxog link), "normal" (with identity link), "lnormal" for lognormal (with log
#' link), "tweedie" (with log link), "exponential" (with log link), "gamma"
#' (with log link), "beta" (with logit link), "ordinal" (cumulative probit
#' regression).
#' 
#' Please see \code{\link{about.distributions}} for information on
#' distributions available in boral overall.
#' @param num.X Number of columns in \code{X}. Defaults to 0, in which case it
#' is assumed that no covariates are included in the model. Recall that no
#' intercept should be included in \code{X}.
#' @param X.ind An matrix of 1s and 0s, indicating whether a particular
#' covariate should be included (1) or excluded (0) in the mean structure of a
#' particular response. The matrix should the number of rows equal to the
#' number of columns in \code{y}, and the number of columns equal to the number
#' of columns in \code{X}. Defaults to \code{NULL}, in which case it is assumed
#' that all covariates are included in the mean structure of all responses
#' i.e., all 1s.
#' @param num.traits Number of columns in the model matrix \code{traits}.
#' Defaults to 0, in which case it is assumed no traits are included in model.
#' Recall that no intercept should be included in \code{traits}.
#' @param which.traits A list of length equal to (number of columns in \code{X}
#' + 1), informing which columns of \code{traits} the column-specific
#' intercepts and each of the column-specific regression coefficients should be
#' regressed against. The first element in the list applies to the
#' column-specific intercept, while the remaining elements apply to the
#' regression coefficients. Each element of \code{which.traits} is a vector
#' indicating which traits are to be used.
#' 
#' For example, if \code{which.traits[[2]] = c(2,3)}, then the regression
#' coefficients corresponding to the first column in \code{X} are regressed
#' against the second and third columns of \code{traits}. If
#' \code{which.traits[[2]][1] = 0}, then the regression coefficients for each
#' column are treated as independent. Please see \code{\link{about.traits}} for
#' more details.
#' 
#' Defaults to \code{NULL}, and used in conjunction with \code{traits} and \cr
#' \code{prior.control$ssvs.traitsindex}.
#' @param row.eff Single element indicating whether row effects are included as
#' fixed effects ("fixed"), random effects ("random") or not included ("none")
#' in the fitted model. If fixed effects, then for parameter identifiability
#' the first row effect is set to zero, which analogous to acting as a
#' reference level when dummy variables are used. If random effects, they are
#' drawn from a normal distribution with mean zero and unknown standard
#' deviation. Defaults to "none".
#' @param row.ids A matrix with the number of rows equal to the number of rows
#' in \code{y}, and the number of columns equal to the number of row effects to
#' be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of
#' row \eqn{i} in \code{y} for random effect eqnj; please see
#' \code{\link{boral}} for more details. for details. Defaults to \code{NULL},
#' so that if \code{row.eff = "none"} then the argument is ignored, otherwise
#' if \cr \code{row.eff = "fixed"} or \code{"random"}, \cr then \code{row.ids =
#' matrix(1:nrow(y), ncol = 1)} i.e., a single, row effect unique to each row.
#' @param offset A matrix with the same dimensions as the response matrix
#' \code{y}, specifying an a-priori known component to be included in the
#' linear predictor during fitting. Defaults to \code{NULL}.
#' @param trial.size Either equal to a single element, or a vector of length
#' equal to the number of columns in y. If a single element, then all columns
#' assumed to be binomially distributed will have trial size set to this. If a
#' vector, different trial sizes are allowed in each column of y. The argument
#' is ignored for all columns not assumed to be binomially distributed.
#' Defaults to 1, i.e. Bernoulli distribution.
#' @param n The number of rows in the response matrix \code{y}.
#' @param p The number of columns in the response matrix \code{y}.
#' @param model.name Name of the text file that the JAGS model is written to.
#' Defaults to \code{NULL}, in which case the default of "jagsboralmodel.txt"
#' is used.
#' @param prior.control A list of parameters for controlling the prior
#' distributions. These include: \itemize{ \item \emph{type:} Vector of four
#' strings indicating the type of prior distributions to use. In order, these
#' are: 1) priors for all column-specific intercepts, row effects, and cutoff
#' points for ordinal data; 2) priors for the latent variable coefficients and
#' correlation parameters. This is ignored for this function; 3) priors for all
#' column-specific coefficients relating to \code{X} (ignored if \code{X =
#' NULL}). When traits are included in the model, this is also the prior for
#' the trait regression coefficients (please see \code{\link{about.traits}} for
#' more information); 4) priors for any dispersion parameters and variance
#' (standard deviation, to be precise) parameters in the model.
#' 
#' For elements 1-3, the prior distributions currently available include: I)
#' ``normal", which is a normal prior with the variance controlled by elements
#' 1-3 in \code{hypparams}; II) ``cauchy", which is a Cauchy prior with
#' variance controlled by elements 1-3 in \code{hypparams}. Gelman, et al.
#' (2008) considers using Cauchy priors with variance \eqn{2.5^2} as weakly
#' informative priors for coefficients in logistic and potentially other
#' generalized linear models; III) ``uniform", which is a symmetric uniform
#' prior with minimum and maximum values controlled by element 1-3 in
#' \code{hypparams}.
#' 
#' For element 4, the prior distributions currently available include: I)
#' ``uniform", which is uniform prior with minimum zero and maximum controlled
#' by element 4 in \code{hypparmas}; II) ``halfnormal", which is half-normal
#' prior with variance controlled by \code{hypparams}; III) ``halfcauchy",
#' which is a half-Cauchy prior with variance controlled by element 4 in
#' \code{hypparams}.
#' 
#' Defaults to the vector \code{c("normal","normal","normal","uniform")}.
#' 
#' \item \emph{hypparams} Vector of four hyperparameters used in the set up of
#' prior distributions. In order, these are: 1) affects the prior distribution
#' for all column-specific intercepts, row effects, and cutoff points for
#' ordinal data; 2) affects the prior distribution for all latent variable
#' coefficients and correlation parameters. This is ignored for this function;
#' 3) affects the prior distribution for column-specific coefficients relating
#' to \code{X} (ignored if \code{X = NULL}). When traits are included in the
#' model, it also affects the prior distribution for the trait regression
#' coefficients; 4) affects the prior distribution for any dispersion
#' parameters, as well as the prior distributions for the standard deviation of
#' the random effects normal distribution if \code{row.eff = "random"}, the
#' standard deviation of the column-specific random intercepts for these
#' columns if more than two of the columns are ordinal, and the standard
#' deviation of the random effects normal distribution for trait regression
#' coefficients when traits are included in the model.
#' 
#' Defaults to the vector \code{c(10, 10, 10, 30)}. The use of normal
#' distributions with mean zero and variance 10 as priors is seen as one type
#' of (very) weakly informative prior, according to
#' \href{https://github.com/stan-dev/stan/wiki/Prior-Choice-RecommendationsPrior
#' choice recommendations}.
#' 
#' \item \emph{ssvs.index:} Indices to be used for stochastic search variable
#' selection (SSVS, George and McCulloch, 1993). Either a single element or a
#' vector with length equal to the number of columns in the implied model
#' matrix \code{X}. Each element can take values of -1 (no SSVS is performed on
#' this covariate), 0 (SSVS is performed on individual coefficients for this
#' covariate), or any integer greater than 0 (SSVS is performed on collectively
#' all coefficients on this covariate/s.)
#' 
#' Please see \code{\link{about.ssvs}} for more information regarding the
#' implementation of SSVS. Defaults to -1, in which case SSVS is not performed
#' on \code{X} variables.
#' 
#' \item \emph{ssvs.g:} Multiplicative, shrinkage factor for SSVS, which
#' controls the strength of the "spike" in the SSVS mixture prior. In summary,
#' if the coefficient is included in the model, the "slab" prior is a normal
#' distribution with mean zero and variance given by element 3 in
#' \code{hypparams}, while if the coefficient is not included in the model, the
#' "spike" prior is normal distribution with mean zero and variance given by
#' element 3 in \code{hypparams} multiplied by \code{ssvs.g}. Please see
#' \code{\link{about.ssvs}} for more information regarding the implementation
#' of SSVS. Defaults to 1e-6.
#' 
#' \item \emph{ssvs.traitsindex:} Used in conjunction with \code{traits} and
#' \code{which.traits}, this is a list of indices to be used for performing
#' SSVS on the trait coefficients. Should be a list with the same length as
#' \code{which.traits}, and with each element a vector of indices with the same
#' length as the corresponding element in \code{which.traits}. Each index
#' either can take values of -1 (no SSVS on this trait coefficient) or 0 (no
#' SSVS on this trait coefficient).
#' 
#' Please see \code{\link{about.ssvs}} for more information regarding the
#' implementation of SSVS. Defaults to -1, in which case SSVS is not performed
#' on any of the trait coefficients, if they are included in the model.  }
#' @return A text file is created, containing the JAGS model to be called by
#' the boral function for entering into jags. This file is automatically
#' deleted once boral has finished running unless \code{save.model = TRUE}.
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{make.jagsboralmodel}} for writing JAGS scripts for
#' models with one or more latent variables.
#' @references \itemize{ \item Gelman, et al. (2008). A weakly informative
#' default prior distribution for logistic and other regression models. The
#' Annals of Applied Statistics, 2, 1360-1383. }
#' @examples
#' 
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' n <- nrow(y)
#' p <- ncol(y)
#' 
#' ## Create a "null" model JAGS script, where distributions alternative 
#' ## between Poisson and negative distributions 
#' ##   across the rows of y.
#' make.jagsboralnullmodel(family = rep(c("poisson","negative.binomial"),length=p), 
#'     num.X = ncol(spider$x), row.eff = "fixed", n = n, p = p)
#' 
#'      
#' ## Create a "null" model JAGS script, where distributions are all negative 
#' ## 	binomial distributions and covariates will be included!
#' make.jagsboralnullmodel(family = "negative.binomial", 
#'     num.X = ncol(spider$x), n = n, p = p, 
#'     model.name = "myawesomeordnullmodel.txt")
#' 
#'      
#' ## Have a look at the JAGS model file for a model involving traits,
#' ## based on the ants data from mvabund.
#' library(mvabund)
#' data(antTraits)
#' 
#' y <- antTraits$abun
#' X <- as.matrix(antTraits$env)
#' ## Include only traits 1, 2, and 5, plus an intercept
#' traits <- as.matrix(antTraits$traits[,c(1,2,5)])
#' ## Please see help file for boral regarding the use of which.traits
#' example_which_traits <- vector("list",ncol(X)+1)
#' for(i in 1:length(example_which_traits)) 
#'      example_which_traits[[i]] <- 1:ncol(traits)
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#' 
#' fit_traits <- boral(y, X = X, traits = traits, which.traits = example_which_traits, 
#'     family = "negative.binomial", model.name = "anttraits.txt",
#'     mcmc.control = example_mcmc_control)
#' }
#' 
NULL





#' Plots of a fitted boral object
#' 
#' Produces four plots relating to the fitted boral object, which can be used
#' for residual analysis. If some of the columns are ordinal, then a single
#' confusion matrix is also produced.
#' 
#' Four plots are provided: \enumerate{ \item Plot of Dunn-Smyth residuals
#' against the linear predictors. This can be useful to assess whether the
#' assumed mean-variance relationship is adequately satisfied, as well as to
#' look for particular outliers. For ordinal responses things are more
#' ambiguous due to the lack of single definition for "linear predictor".
#' Therefore, instead of linear predictors the Dunn-Smyth residuals are plotted
#' against the fitted values (defined as the level with the highest fitted
#' probability). It is fully acknowledged that this makes things VERY hard to
#' interpret if only some of your columns are ordinal.  \item Plot of
#' Dunn-Smyth residuals against the row index/row names.  \item Plot of
#' Dunn-Smyth residuals against the column index/column names. Both this and
#' the previous plot are useful for assessing how well each row/column of the
#' response matrix is being modeled.  \item A normal quantile plot of the
#' Dunn-Smyth residuals, which can be used to assess the normality assumption
#' and overall goodness of fit.  } For ordinal responses, a single confusion
#' matrix between the predicted levels (as based on the class with the highest
#' probability) and true levels is aso returned. The table pools the results
#' over all columns assumed to be ordinal.
#' 
#' @name plot.boral
#' @docType package
#' @param x An object of class "boral".
#' @param est A choice of either the posterior median (\code{est == "median"})
#' or posterior mean (\code{est == "mean"}) of the parameters, which are then
#' treated as parameter estimates and the fitted values/residuals used in the
#' plots are calculated from. Default is posterior median.
#' @param jitter If \code{jitter = TRUE}, then some jittering is applied so
#' that points on the plots do not overlap exactly (which can often occur with
#' discrete data). Please see \code{\link{jitter}} for its implementation.
#' @param ... Additional graphical options to be included in. These include
#' values for \cr \code{cex, cex.lab, cex.axis, cex.main, lwd}, and so on.
#' @note Due the inherent stochasticity, Dunn-Smyth residuals and consequently
#' the plots will be slightly different time this function is run. Note also
#' the fitted values and residuals are calculated from point estimates of the
#' parameters, as opposed to a fully Bayesian approach (please see details in
#' \code{\link{fitted.boral}} and \code{\link{ds.residuals}}). Consequently, it
#' is recommended that this function is run several times to ensure that any
#' trends observed in the plots are consistent throughout the runs.
#' 
#' As mentioned above, for ordinal responses things are much more challenging
#' as there is no single definition for "linear predictor". Instead of linear
#' predictors then, for the first plot the Dunn-Smyth residuals are plotted
#' against the fitted values, defined as the level with the highest fitted
#' probability. It is fully acknowledged that this makes things VERY hard to
#' interpret if only some of your columns are ordinal though. Suggestions to
#' improve this are welcome!!!
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{fitted.boral}} to obtain the fitted values,
#' \code{\link{ds.residuals}} to obtain Dunn-Smyth residuals and details as to
#' what they are.
#' @examples
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'     n.thin = 1)
#' 
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' 
#' spider.fit.p <- boral(y, family = "poisson", lv.control = list(num.lv = 2),
#'     row.eff = "fixed", mcmc.control = example_mcmc_control)
#' 
#' par(mfrow = c(2,2))
#' plot(spider.fit.p) 
#' ## A distinct fan pattern is observed in the plot of residuals 
#' ## versus linear predictors plot. 
#' 
#' 
#' spiderfit_nb <- boral(y, family = "negative.binomial", lv.control = list(num.lv = 2), 
#'     row.eff = "fixed", mcmc.control = example_mcmc_control)
#' 
#' par(mfrow = c(2,2))
#' plot(spiderfit_nb) 
#' ## The fan shape is not as clear now, 
#' ## and the normal quantile plot also suggests a better fit to the data 
#' }
#' 
NULL





#' Predict using a model
#' 
#' Obtain predictions and associated intervals (lower and upper limits) on the
#' linear scale from a fitted boral object. Predictions can be made either
#' conditionally on the predicted latent variables and any random row effects
#' included in the model, or marginally (averaged) on the latent variables and
#' any random effects included in the model.
#' 
#' Due to the Bayesian MCMC framework, then predictive inference for models is
#' based around the posterior predictive distribution, which is the integral of
#' the quantity one wants to predict on, integrated or averaged over the
#' posterior distribution of the parameters and latent variables. Currently,
#' all predictions are made on the \emph{linear predictor scale} i.e.,
#' 
#' \deqn{\eta_{ij} = \alpha_i + \beta_{0j} + \bm{x}^\top_i\bm{\beta}_j +
#' \bm{z}^\top_i\bm{\theta}_j; \quad i = 1,\ldots,n; j = 1,\ldots,p,}
#' 
#' where \eqn{\bm{z}_i} are a vector of latent variables included in the model,
#' \eqn{\bm{\theta}_j} are the column-specific coefficients relating to these
#' latent variables, \eqn{\bm{x}_i} are covariates included in the model, and
#' \eqn{\bm{beta}_j} being the column-specific coefficients related to these
#' covariates. The quantity \eqn{beta_{0j}} denotes the column-specific
#' intercepts while \code{alpha_i} represents one or more optional row effects
#' that may be treated as a fixed or random effect.
#' 
#' Note that for the above to work, one must have saved the MCMC samples in the
#' fitted boral object, that is, set \code{save.model = TRUE} when fitting.
#' 
#' Two types of predictions are possible using this function: \itemize{ \item
#' The first type is \code{predict.type = "conditional"}, meaning predictions
#' are made conditionally on the predicted latent variables and any (random)
#' row effects in the model. This is mainly used when predictions are made onto
#' the \emph{same} set of sites that the model was fitted to, although a
#' \code{newX} can be supplied in this case if we want to extrapolate on to the
#' same set of sites but under different environmental conditions.
#' 
#' \item The second type of prediction is \code{predict.type = "marginal"},
#' meaning predictions are made marginally or averaging over the latent
#' variables and any (random) row effects in the model. This is mainly used
#' when predictions are made onto a \emph{new} set of sites where the latent
#' variables and/or row effects are unknown. A \code{newX} and/or
#' \code{newrow.ids} is often supplied since we are extrapolating to new sites.
#' The integration over the latent variables and random row effects is done via
#' Monte-Carlo integration. Please note however that, as mentioned before, the
#' integration will be done on the linear predictor scale.  }
#' 
#' More information on conditional versus marginal predictions in latent
#' variable models can be found in Warton et al., (2015). In both cases, the
#' function returns a point prediction (either the posterior mean or median
#' depending on \code{est}) and the lower and upper bounds of a
#' 100\eqn{\alpha\%} interval of the posterior prediction. All of these
#' quantities are calculated empirically based the MCMC samples e.g., the
#' posterior mean is the average of the predictions across the MCMC samples,
#' and the lower and upper bounds are based on quantiles.
#' 
#' @name predict.boral
#' @docType package
#' @param object An object of class "boral".
#' @param newX An optional model matrix of covariates for extrapolation to the
#' same sites (under different environmental conditions) or extrapolation to
#' new sites. No intercept column should be included in \code{newX}. Defaults
#' to \code{NULL}, in which case the model matrix of covariates is taken from
#' the fitted boral object if found.
#' @param newrow.ids An optional matrix with the number of columns equal to the
#' number of row effects to be included in the model. Element \eqn{(i,j)}
#' indicates to the cluster ID of row \eqn{i} in \code{y} for random effect
#' eqnj. Defaults to \code{NULL}, in which case row IDs are taken from the
#' fitted boral object itself (if required) i.e., from \code{object$row.ids}.
#' @param distmat A distance matrix required to calculate correlations across
#' sites when a non-independence correlation structure on the latent variables
#' is imposed.
#' @param predict.type The type of prediction to be made. Either takes value
#' \code{"conditional"} in which case the prediction is made conditionally on
#' the predicted latent variables and any random row effects in the model, or
#' \code{"marginal"} in which case the prediction marginalizes (averages) over
#' the latent variables and random row effects in the model. Defaults to
#' \code{"conditional"}.
#' @param scale The type of prediction required.  The default "link" is on the
#' scale of the linear predictors; the alternative \code{scale == "response"}
#' is on the scale of the response variable. Thus for a default binomial family
#' the default predictions provide probabilities on probit scale) and
#' \code{scale == "response"} gives the predicted probabilities.
#' @param est A choice of either whether to print the posterior median
#' (\code{est == "median"}) or posterior mean (\code{est == "mean"}) of the
#' parameters.
#' @param prob A numeric scalar in the interval (0,1) giving the target
#' probability coverage of the intervals. Defaults to 0.95.
#' @param lv.mc If the predictions are made marginalizing over the latent
#' variables, then number of Monte-Carlo samples to take when performing the
#' relevant integration.
#' @param return.alllinpred If \code{TRUE}, then the full array of predicted
#' linear predictions across all MCMC samples is predicted. This is useful if
#' the user wants to transform the predictions onto a different scale (say).
#' Defaults to \code{FALSE}.
#' @param ... Not used.
#' @return A list containing the following components: \item{linpred}{A matrix
#' containing posterior point predictions (either posterior mean or median
#' depending on \code{est}), on the linear predictor scale.}
#' 
#' \item{lower}{A matrix containing the lower bound of the 100\code{alpha}\%
#' interval of the posterior predictions, on the linear predictor scale.}
#' 
#' \item{upper}{A matrix containing the upper bound of the 100\code{alpha}\%
#' interval of the posterior predictions, on the linear predictor scale.}
#' 
#' \item{all.linpred}{If \code{return.alllinpred = TRUE}, then an array of
#' predicted linear predictions across all MCMC samples.}
#' 
#' @section Warnings: \itemize{ \item Marginal predictions can take quite a
#' while to construct due to the need to perform Monte-Carlo integration to
#' marginalize over the latent variables and any random row effects in the
#' model. }
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @references \itemize{ \item Gelman et al. (2013) Bayesian Data Analysis. CRC
#' Press.
#' 
#' \item Warton et al. (2015). So Many Variables: Joint Modeling in Community
#' Ecology. Trends in Ecology and Evolution, 30, 766-779. }
#' @examples
#' 
#' \dontrun{
#' library(mvabund) ## Load a dataset from the mvabund package
#' library(mvtnorm) 
#' data(spider)
#' y <- spider$abun
#' X <- scale(spider$x)
#' n <- nrow(y)
#' p <- ncol(y)
#' 
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'      n.thin = 1)
#' 
#' ## Example 1 - model with two latent variables, random site effects, 
#' ## 	and environmental covariates
#' spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
#'     row.eff = "random", lv.control = list(num.lv = 2), 
#'     mcmc.control = example_mcmc_control, save.model = TRUE)
#' 
#' 
#' ## Predictions conditional on predicted latent variables
#' getcondpreds <- predict(spiderfit_nb)
#' 
#' ## Predictions marginal on latent variables, random row effects
#' ## The intervals for these will generally be wider than the
#' ##   conditional intervals.
#' getmargpreds <- predict(spiderfit_nb, predict.type = "marginal")
#' 
#' 
#' ## Now suppose you extrpolate to new sites
#' newX <- rmvnorm(100, mean = rep(0,ncol(X)))
#' 
#' ## Below won't work since conditional predictions are made to the same sites
#' getcondpreds <- predict(spiderfit_nb, newX = newX)
#' 
#' ## Marginal predictions will work though provided newrow.ids is set up 
#' ## properly. For example,
#' new_row_ids <- matrix(sample(1:28,100,replace=TRUE), 100, 1)
#' getmargpreds <- predict(spiderfit_nb, newX = newX, predict.type = "marginal", 
#'      newrow.ids = new_row_ids)
#' 
#'      
#' ## Example 1b - Similar to 1a except with no random site effects, 
#' ## 	and a non-independence correlation structure for the latent variables
#' ##      based on a fake distance matrix
#' fakedistmat <- as.matrix(distmat(1:n))
#' spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
#'      lv.control = list(type = "squared.exponential", num.lv = 2, distmat = fakedistmat),
#'      mcmc.control = example_mcmc_control, save.model = TRUE)
#' 
#' getmargpreds <- predict(spiderfit_nb, predict.type = "marginal", distmat = fakedistmat)
#' 
#' ## Now suppose you extrpolate to new sites
#' newfakedistmat <- as.matrix(distmat(1:100))
#' 
#' getmargpreds <- predict(spiderfit_nb, newX = newX, predict.type = "marginal", 
#'      distmat = newfakedistmat)
#' 
#'      
#'      
#' ## Example 2 - simulate count data, based on a model with two latent variables, 
#' ## no site variables, with two traits and one environmental covariates 
#' library(mvtnorm)
#' 
#' n <- 100; s <- 50
#' X <- as.matrix(scale(1:n))
#' colnames(X) <- c("elevation")
#' 
#' traits <- cbind(rbinom(s,1,0.5), rnorm(s)) 
#' ## one categorical and one continuous variable
#' colnames(traits) <- c("thorns-dummy","SLA")
#' 
#' simfit <- list(true.lv = rmvnorm(n, mean = rep(0,2)), 
#' 	lv.coefs = cbind(rnorm(s), rmvnorm(s, mean = rep(0,2)), 1), 
#' 	traits.coefs = matrix(c(0.1,1,-0.5,0.1,0.5,0,-1,0.1), 2, byrow = TRUE))
#' rownames(simfit$traits.coefs) <- c("beta0","elevation")
#' colnames(simfit$traits.coefs) <- c("kappa0","thorns-dummy","SLA","sigma")
#' 
#' simy = create.life(true.lv = simfit$true.lv, lv.coefs = simfit$lv.coefs, X = X, 
#' 	traits = traits, traits.coefs = simfit$traits.coefs, family = "normal") 
#' 
#' 
#' example_which_traits <- vector("list",ncol(X)+1)
#' for(i in 1:length(example_which_traits)) 
#'      example_which_traits[[i]] <- 1:ncol(traits)
#' fit_traits <- boral(y = simy, X = X, traits = traits, 
#'      which.traits = example_which_traits, family = "normal", 
#'      lv.control = list(num.lv = 2), save.model = TRUE, 
#'      mcmc.control = example_mcmc_control)	
#' 
#'      
#' ## Predictions conditional on predicted latent variables   
#' getcondpreds <- predict(fit_traits)     
#'      
#' ## Predictions marginal on latent variables
#' ## The intervals for these will generally be wider than the
#' ##   conditional intervals.
#' getmargpreds <- predict(fit_traits, predict.type = "marginal")
#' }
#' 
#' 
NULL





#' Summary of fitted boral object
#' 
#' A summary of the fitted boral objects including the type of model fitted
#' e.g., error distribution, number of latent variables parameter estimates,
#' and so on.
#' 
#' 
#' @name summary.boral
#' @aliases summary.boral print.summary.boral
#' @docType package
#' @param object An object of class "boral".
#' @param x An object of class "boral".
#' @param est A choice of either whether to print the posterior median
#' (\code{est == "median"}) or posterior mean (\code{est == "mean"}) of the
#' parameters.
#' @param ... Not used.
#' @return Attributes of the model fitted, parameter estimates, and posterior
#' probabilities of including individual and/or grouped coefficients in the
#' model based on SSVS if appropriate.
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_author("boral")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "boral")\Sexpr{tools:::Rd_package_maintainer("boral")}
#' @seealso \code{\link{boral}} for the fitting function on which
#' \code{summary} is applied.
#' @examples
#' 
#' \dontrun{
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'     n.thin = 1)
#' 
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' 
#' spiderfit_nb <- boral(y, family = "negative.binomial", lv.control = list(num.lv = 2),
#'     row.eff = "fixed", mcmc.control = example_mcmc_control)
#' 
#' summary(spiderfit_nb)
#' }
#' 
NULL








