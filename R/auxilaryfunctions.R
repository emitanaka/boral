
# auxiliary functions -----------------------------------------------------

#' Variance partitioning for a latent variable model
#' 
#' For each response (species), partition the variance of the linear predictor 
#' into components associated with (groups of) the covariates, the latent 
#' variables, any row effects. If traits are also included in the model, then 
#' it also calculates an R-squared value for the proportion of the variance in 
#' the environmental response (due to the covariates) which can be explained by 
#' traits.
#' 
#' Variance partitioning on the linear predictor scale; the variance is 
#' calculated and and averaged wrt to the posterior distribution 
#' 
#' #' @usage  
#' \dontrun{calc.varpart(object, groupX = NULL)}
#' 
#' @param object An object of class "boral".
#' @param groupX A vector of group indicator variables, which allows the 
#' variance partitioning to be done for groups of covariates (including the 
#' intercept) i.e., how much of the total variation does a certain subset of 
#' the covariates explain. Defaults to \code{NULL}, in whih case all the 
#' covariates are treated as single group.
#' 
#' @details 
#' As an alternative to looking at differences in trace of the residual covariance matrix (Hui et al., 2014; Warton et al., 2015), an alternative way to quantify the amount of variance explained by covariates, traits, row effects, is to perform a variance decomposition of the linear predictor of a latent variable model. In particular, for a general model the linear predictor for response \eqn{j = 1,\ldots,p} at row \eqn{i = 1,\ldots,n} is given by
#' \deqn{\eta_{ij} = \alpha_i + \beta_{0j} + \bm{x}^\top_i\bm{\beta}_j, + \bm{z}^\top_i\bm{\theta}_j,}
#' where \eqn{\beta_{0j} + \bm{x}^\top_i\bm{\beta}_j} is the component of the linear predictor due to the covariates \code{X} plus an intercept, \eqn{\bm{z}^\top_i\bm{\theta}_j} is the component due to the latent variables, and \eqn{\alpha_i} is the component due to one or more fixed or random row effects. The regression coefficients \eqn{\bm{\beta}_j} may be further as random effects and regressed against traits; please see \code{\link{about.traits}} for further information on this. 
#' For the response, a variation partitioning of the linear is performed by calculating the variance due to component in \eqn{\eta_{ij}} and then rescaling them to ensure that they sum to one. The general details of this type of variation partitioning is given in Ovaskainen et al., (2017); see also Nakagawa and Schielzeth (2013) for R-squared and proportion of variance explained in the case of generalized linear mixed model. In brief, for response \eqn{j = 1,\ldots,p}: 1) the variance due to the \code{X} covariates and intercept is given by the variance of \eqn{\beta_{0j} + \bm{x}^\top_i\bm{\beta}_j} calculated across the \eqn{n} rows; 2) the variance due the latent variables is given by the diagonal elements of \eqn{\bm{\theta}^\top_j\bm{\theta}_j}; 3) the variance due to (all) the random effects is given by variance of \eqn{\alpha_i} calculated across the \eqn{n} rows for fixed row effects (\code{row.eff = "fixed"}), and given by the (sum of the) variance \eqn{\sigma^2_{\alpha}} for random row effects (\code{row.eff = "random"}). After scaling, we can then obtain the proportion of variance for each response which is explained by the variance components. These proportions are calculated for each MCMC sample and then acrossed them to calculate a posterior mean variance partitioning. 
#' If \code{groupX} is supplied, the variance due to the \code{X} covariates is done based on subsets of \code{X} variables (including the intercept) as identified by code{groupX}, and then rescaled correspondingly. This is useful if one was to, for example, quantify the proportion of variation in each species which is explained by each \code{X} covariate.
#' If a fitted model also containing traits, which are included to help explain/mediate differences in species environmental responses, then the function calculates \eqn{R^2} value for the proportion of variance in the \code{X} variables which is explained by the traits. In brief, this is calculated based the correlation between \eqn{\beta_{0j} + \bm{x}^\top_i\bm{\beta}_j} and \eqn{\tau_{0j} + \bm{x}^\top_i\bm{\tau}_j}, where \eqn{\tau_{0j}} and \eqn{\bm{\tau}_j} are the ``predicted" values of the species coefficients based on values i.e., \eqn{\tau_{0j} = \kappa_{01} + \bm{traits}^\top_j\bm{\kappa}_1} and \eqn{\tau_{jk} = \kappa_{0k} + \bm{traits}^\top_j\bm{\kappa}_k} for element \eqn{k} in \eqn{\bm{\tau}_j}. 
#' 
#' @return 
#' A list containing the following components if applicable:
#' \describe{
#' \item{varpart.X}{Vector containing the proportion of variance (in the linear predictor) for each response which is explained by \code{X}.}
#' \item{varpart.lv}{Vector containing the proportion of variance (in the linear predictor) for each response which is explained by the latent variables.}
#' \item{varpart.row}{Vector containing the proportion of variance (in the linear predictor) for each response which is explained by the row effects.}
#' \item{R2.traits}{Vector containing the proportion of variance due to the covariates for each response, which can be explained by traits for each response.}
#' }
#' @references 
#' \itemize{
#' \item Nakagawa, S., and Schielzeth, H. (2013). A general and simple method for obtaining R2 from generalized linear mixed-effects models. Methods in Ecology and Evolution 4, 133-142.
#' \item Ovaskainen, et al. (2017). How to make more out of community data? A conceptual framework and its implementation as models and software. Ecology Letters 20, 561-576.
#' \item Hui et al. (2014). Model-based approaches to unconstrained ordination. Methods in Ecology and Evolution, 6, 399-411.
#' \item Warton et al. (2015). So Many Variables: Joint Modeling in Community Ecology. Trends in Ecology and Evolution, 30, 766-779.}
#' 
#' \section{Warnings}{
#' There is considerable controversy over exactly what quantities such as R-squared and proportion of variance explained are in the case mixed models and latent variable models, and how they can interpreted e.g., what is considered a high value for the proportion of variance by \code{X} variables, is it consistent with whether the coefficients are significantly different from zero or not; see for instance \href{https://stats.stackexchange.com/questions/111150/calculating-r2-in-mixed-models-using-nakagawa-schielzeths-2013-r2glmm-me}{R2 controversy}.
#' 
#' When reporting these values, researchers should be at least aware of this and that there are multiple ways of manufacturing such quantities with no single best approach e.g., using relative changes in trace of the residual covariance matrix, relative changes in marginal and conditional log-likelihoods are other possible approaches. 
#' }
#' @examples 
#' \dontrun{
#' library(mvabund) ## Load a dataset from the mvabund package
#' data(spider)
#' y <- spider$abun
#' X <- scale(spider$x)
#' n <- nrow(y)
#' p <- ncol(y)
#' 
#' ## NOTE: The values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
#'                           n.thin = 1)
#'                           
#' ## Example 1 - model with X variables, two latent variables, and no row effects
#' spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
#' lv.control = list(num.lv = 2), 
#' save.model = TRUE, mcmc.control = example_mcmc_control)
#' ## Partition variance for each species into that explained by covariates 
#' ## and by the latent variables
#' dovar <- calc.varpart(spiderfit_nb)
#' 
#' ## Consider the intercept and first two covariates in X as one group, 
#' ## and remaining four covariates in X as another group, 
#' ## then partition variance for each species based on these groups.
#' dovar <- calc.varpart(spiderfit_nb, groupX = c(1,1,1,2,2,2,2))
#' 
#' 
#' ## Example 2 - model fitted to count data, no site effects, and
#' ## two latent variables, plus traits included to explain environmental responses
#' data(antTraits)
#' y <- antTraits$abun
#' X <- as.matrix(scale(antTraits$env))
#' ## Include only traits 1, 2, and 5
#' traits <- as.matrix(antTraits$traits[,c(1,2,5)])
#' example_which_traits <- vector("list",ncol(X)+1)
#' for(i in 1:length(example_which_traits)) 
#'   example_which_traits[[i]] <- 1:ncol(traits)
#'   ## Just for fun, the regression coefficients for the second column of X,
#'   ## corresponding to the third element in the list example_which_traits,
#'   ## will be estimated separately and not regressed against traits.
#'   example_which_traits[[3]] <- 0
#'   
#'   fit_traits <- boral(y, X = X, traits = traits, which.traits = example_which_traits, 
#'                  family = "negative.binomial", mcmc.control = example_mcmc_control, 
#'                                save.model = TRUE)
#'                                
#' ## Partition variance for each species due to covariates in X 
#' ## and latent variables. Also calculate proportion of variance 
#' ## due to the covariates which can be explained by traits 
#' dovar <- calc.varpart(fit_traits)}
#' 
calc.varpart <- function(object, groupX = NULL) {
    fit.mcmc <- get.mcmcsamples(object)
    num.lv <- object$num.lv
    if(is.null(object$X))
        stop("Variance partitioning is not particularly useful when there are no covariates X included in the model")
    if(is.null(object$jags.model)) 
        stop("MCMC samples not found")
    if(!is.null(groupX)) { if(length(groupX) != (object$num.X+1)) 
        stop("If groupX is supplied, then it must be a vector with the length (object$num.X + 1)") }
       
      row_var <- lv_var <- X_var <- matrix(0, nrow(fit.mcmc), object$p)
  if(!is.null(groupX))
  groupX_var <- array(0, dim = c(length(unique(groupX)), nrow(fit.mcmc), object$p))          
  if(!is.null(object$traits)) {
  trait.linpred <- array(0, dim = c(nrow(fit.mcmc), object$n, object$p))
  all_cors_spp <- matrix(0, nrow(fit.mcmc), object$p)
  }

        
    for(k in 1:nrow(fit.mcmc)) {
        cw_X_coefs <- matrix(fit.mcmc[k, grep("X.coefs", colnames(fit.mcmc))], nrow = object$p)
        cw_lv_coefs <- matrix(fit.mcmc[k, grep("lv.coefs", colnames(fit.mcmc))], nrow = object$p) ## Need for spp intercept

        fullX <- cbind(1,object$X)
        full.coefs <- cbind(cw_lv_coefs[,1],as.matrix(cw_X_coefs))
        cw.X.linpred <- tcrossprod(fullX, full.coefs)
        if(!is.null(groupX)) {
            cw.groupX.linpred <- vector("list", length(unique(groupX)))
            for(k2 in 1:length(unique(groupX))) 
                cw.groupX.linpred[[k2]] <- tcrossprod(fullX[,which(groupX==k2)], full.coefs[,which(groupX==k2)])
            }
        if(object$row.eff == "random") {
            for(k2 in 1:ncol(object$row.ids)) 
                row_var[k,] <- row_var[k,] + rep(fit.mcmc[k, grep(paste0("row.sigma.ID",k2,"$"), colnames(fit.mcmc))]^2,object$p) 
            }
        if(object$row.eff == "fixed") {
            cw.row.linpred <- matrix(0, object$n, object$p)
            for(k2 in 1:ncol(object$row.ids)) 
                cw.row.linpred <- cw.row.linpred + matrix(fit.mcmc[k, grep(paste0("row.coefs.ID",k2,"\\["), colnames(fit.mcmc))][object$row.ids[,k2]], nrow = object$n, ncol = object$p, byrow = FALSE)
            row_var[k,] <- apply(cw.row.linpred, 2, var)
            }

        X_var[k,] <- apply(cw.X.linpred,2,var)
        if(!is.null(groupX)) {
            for(k2 in 1:length(unique(groupX))) 
                groupX_var[k2,k,] <- apply(cw.groupX.linpred[[k2]],2,var)
            }
        if(num.lv > 0) 
            lv_var[k,] <- rowSums(cw_lv_coefs[,2:(1+object$num.lv)]^2) 
        if(!is.null(object$traits)) {
            cw.traits.coefs <- cbind(fit.mcmc[k, grep("traits.int",colnames(fit.mcmc))], matrix(fit.mcmc[k, grep("traits.coefs",colnames(fit.mcmc))], nrow = ncol(object$X)+1))
            rownames(cw.traits.coefs) <- c("beta0", colnames(object$X))
            trait.X.coefs <- tcrossprod(cbind(1,object$traits), cw.traits.coefs) ## beta = intercept + trait %*% trait.coefs
            cw.trait.linpred <- tcrossprod(cbind(1,object$X), trait.X.coefs)
            all_cors_spp[k,] <-  sapply(1:object$p, function(i) cor(cw.X.linpred[,i], cw.trait.linpred[,i])^2)
            }
        }
    

    total.var <- X_var + row_var + lv_var
    var.X <- colMeans(X_var/total.var); names(var.X) <- colnames(object$y)
    var_lv <- NULL
    if(num.lv > 0) { 
        var_lv <- colMeans(lv_var/total.var)
        names(var_lv) <- colnames(object$y) 
        }
    var_row <- NULL
    if(object$row.eff != "none") { 
        var_row <- colMeans(row_var/total.var)
        names(var_row) <- colnames(object$y) 
        }

    ## As soon as groupX is supplied, change the way variance decomposition is done
    if(!is.null(groupX)) {
        total.var <- apply(groupX_var, c(2,3), sum) + row_var + lv_var ## Note this is not equal to total.var
        var.X <- matrix(0, length(unique(groupX)), object$p) 
        for(k2 in 1:length(unique(groupX))) 
            var.X[k2,] <- colMeans(groupX_var[k2,,]/total.var)
        rownames(var.X) <- unique(groupX)
        colnames(var.X) <- colnames(object$y)
        var_lv <- NULL
        if(num.lv > 0) { 
            var_lv <- colMeans(lv_var/total.var)
            names(var_lv) <- colnames(object$y) 
            }
        var_row <- NULL
        if(object$row.eff != "none") { 
            var_row <- colMeans(row_var/total.var)
            names(var_row) <- colnames(object$y) 
            }          
        }
    
    out <- list(varpart.X = var.X, varpart.lv = var_lv, varpart.row = var_row)
    if(!is.null(object$traits)) {
        out$R2.traits <- colMeans(all_cors_spp)
        names(all_cors_spp) <- colnames(object$y) 
        }
    return(out)
    }


## Dunn-Smyth residuals
## Also create a confusion matrix for ordinal and multinomial data
ds.residuals <- function(object, est = "median") {  
    n <- object$n; p <- object$p; 
    num.lv <- object$num.lv
    num.ord.levels <- object$num.ord.levels; 
    X <- object$X; y <- object$y
    mus <- fitted.boral(object, est = est)

    if(any(object$family == "ordinal")) {
        message("One or more columns of y have ordinal responses. Constructing a single confusion matrix for these")
        true_resp <- as.matrix(y[,which(object$family == "ordinal")])
        pred_resp <- matrix(NA,n,ncol(true_resp)) 
        }
# 	if(any(object$family == "multinom")) {
# 		print("One or more columns of y have multinomial responses. Constructing a single confusion matrix for these")
# 		true.multinom.resp <- as.matrix(y[,which(object$family == "multinom")])
# 		pred.multinom.resp <- matrix(NA,n,ncol(true.multinom.resp)) }
    if(any(object$family == "tweedie")) {
        if(est == "median") 
            powerparam <- object$powerparam.median
        if(est == "mean") 
            powerparam <- object$powerparam.mean 
        }

    dsres_out <- matrix(NA,n,p)
    rownames(dsres_out) <- rownames(y)
    colnames(dsres_out) <- colnames(y)
    for(i in 1:n) { for(j in 1:p) {
        if(object$family[j] == "poisson") { 
            a <- ppois(as.vector(unlist(y[i,j]))-1, mus$out[i,j]); 
            b <- ppois(as.vector(unlist(y[i,j])), mus$out[i,j]); 		
            u <- runif(n = 1, min = a, max = b); dsres_out[i,j] <- qnorm(u) 
            }
        if(object$family[j] == "negative.binomial") {
            if(est == "median") 
                phis <- object$lv.coefs.median[,num.lv+2]+1e-5
            if(est == "mean") 
                phis <- object$lv.coefs.mean[,num.lv+2]+1e-5
            a <- pnbinom(as.vector(unlist(y[i,j]))-1, mu=mus$out[i,j], size=1/phis[j]); 
            b <- pnbinom(as.vector(unlist(y[i,j])), mu=mus$out[i,j], size=1/phis[j])
            u <- runif(n = 1, min = a, max = b); dsres_out[i,j] <- qnorm(u) 
            }
        if(object$family[j] == "binomial") { 
            a <- pbinom(as.vector(unlist(y[i,j]))-1, object$trial.size[j], prob = mus$out[i,j]); 
            b <- pbinom(as.vector(unlist(y[i,j])), object$trial.size[j], prob = mus$out[i,j])
            u <- runif(n = 1, min = a, max = b); dsres_out[i,j] <- qnorm(u) 
            }
        if(object$family[j] == "exponential") { 
            a <- pexp(as.vector(unlist(y[i,j])), rate=1/mus$out[i,j]); 
            dsres_out[i,j] <- qnorm(a) 
            }
        if(object$family[j] == "gamma") { 
            if(est == "median") 
                phis <- object$lv.coefs.median[,num.lv+2]
            if(est == "mean") 
                phis <- object$lv.coefs.mean[,num.lv+2]
            a <- pgamma(as.vector(unlist(y[i,j])), shape=mus$out[i,j]*phis[j], rate=phis[j]); 
            dsres_out[i,j] <- qnorm(a) 
            }
        if(object$family[j] == "beta") { 
            if(est == "median") 
                phis <- object$lv.coefs.median[,num.lv+2]
            if(est == "mean") 
                phis <- object$lv.coefs.mean[,num.lv+2]
            a <- pbeta(as.vector(unlist(y[i,j])), shape1=phis[j]*mus$out[i,j], shape2=phis[j]*(1-mus$out[i,j]))
            dsres_out[i,j] <- qnorm(a) 
            }
        if(object$family[j] == "normal") { 
            if(est == "median") 
                phis <- object$lv.coefs.median[,num.lv+2]
            if(est == "mean") 
                phis <- object$lv.coefs.mean[,num.lv+2]
            a <- pnorm(as.vector(unlist(y[i,j])), mus$out[i,j], sd = (phis[j])); 
            dsres_out[i,j] <- qnorm(a) 
            }
# 			X2 <- cbind(1,X); hatmat <- X2%*%solve(crossprod(X2))%*%t(X2)
# 			dsres_out[i,j] <- (y[i,j]-mus$out[i,j])/(sqrt(phis[j])*sqrt(1-hatmat[i,i])) }
        if(object$family[j] == "lnormal") { 
            if(est == "median") 
                phis <- object$lv.coefs.median[,num.lv+2]
            if(est == "mean") 
                phis <- object$lv.coefs.mean[,num.lv+2]
            a <- plnorm(as.vector(unlist(y[i,j])), log(mus$out[i,j]), sdlog = (phis[j])); dsres_out[i,j] <- qnorm(a) 
            }
        if(object$family[j] == "tweedie") { 
            if(est == "median") 
                phis <- object$lv.coefs.median[,num.lv+2]
            if(est == "mean") 
                phis <- object$lv.coefs.mean[,num.lv+2]
            a <- pTweedie(as.vector(unlist(y[i,j])), mu = mus$out[i,j], phi = phis[j], p = powerparam); dsres_out[i,j] <- qnorm(a) 
            }
        if(object$family[j] == "ordinal") { 
            pred_resp[,which(object$family == "ordinal")==j] <- mus$out[,which(object$family == "ordinal")==j] ## get max predicted probability
            cumsum.b <- sum(mus$ordinal.probs[i,j,1:(y[i,j])])
            cumsum.a <- sum(mus$ordinal.probs[i,j,1:(y[i,j]-1)])
            u <- runif(n = 1, min = cumsum.a, max = cumsum.b); 
            if(abs(u-1) < 1e-5) 
                u <- 1
            if(abs(u-0) < 1e-5) 
                u <- 0
            dsres_out[i,j] <- qnorm(u) 
            }
# 		if(object$family[j] == "multinom") { ## get max predicted probability
# 			pred_resp[i,which(object$family == "multinom")==j] <- which.max(mus$multinom.probs[i,j,]) }
        } }

    if(sum(object$family == "ordinal") > 0) {
        agree_tab <- table(as.vector(pred_resp), as.vector(true_resp))
        }
    else { 
        agree_tab <- NULL 
        }
    #if(sum(object$family == "multinom") > 0) { agree.multinom.tab <- table(as.vector(pred.multinom.resp), as.vector(true.multinom.resp)); }	else { agree.multinom.tab <- NULL }
    
    return(list(agree.ordinal = agree_tab, residuals = dsres_out))
    }

	
## Fitted values
## For ordinal and multinomial data, returns a matrix of probabilities for each vector of rows
fitted.boral <- function(object, est = "median",...) {
    n <- object$n; p <- object$p
    num.lv <- object$num.lv; 
    X <- object$X
    y <- object$y
    fitted_out <- matrix(NA,n,p)
    rownames(fitted_out) <- rownames(y)
    colnames(fitted_out) <- colnames(y)

    if(any(object$family == "ordinal")) { 
        fitted_ordinal_probs <- array(NA, dim=c(n,p,object$num.ord.levels)) 
        dimnames(fitted_ordinal_probs) <- list(r = rownames(y), c = colnames(y), levels = 1:object$num.ord.levels) 
        } 
    else { 
        fitted_ordinal_probs <- NULL 
        }
# 	if(any(object$family == "multinom")) { 
# 		fitted_multinom_probs <- array(NA,dim=c(n,p,object$num.multinom.levels))
# 		dimnames(fitted_multinom_probs) <- list(r = rownames(y), c = colnames(y), levels = 1:object$num.multinom.levels) } 
# 	else { fitted_multinom_probs <- NULL }
    if(length(object$lv.coefs.median) == p) {
        object$lv.coefs.median <- matrix(object$lv.coefs.median, ncol = 1)
        object$lv.coefs.mean <- matrix(object$lv.coefs.mean, ncol = 1)
        }
    if(length(object$X.coefs.median) == p) {
        object$X.coefs.median <- matrix(object$X.coefs.median, ncol = 1)
        object$X.coefs.mean <- matrix(object$X.coefs.mean, ncol = 1)
        }

        
    if(is.null(object$lv.median)) 
        eta <- tcrossprod(matrix(1,n,1), matrix(object$lv.coefs.median[,1],p,1))
    if(!is.null(object$lv.median)) 
        eta <- tcrossprod(cbind(1,object$lv.median), object$lv.coefs.median[,1:(num.lv+1)])
    if(!is.null(object$X.coefs.median)) 
        eta <- eta + tcrossprod(as.matrix(X), object$X.coefs.median)
    if(!is.null(object$offset)) 
        eta <- eta + object$offset
    
    if(est == "mean") {
        if(is.null(object$lv.mean)) 
            eta <- tcrossprod(matrix(1,n,1), object$lv.coefs.mean[,1:(num.lv+1)])
        if(!is.null(object$lv.mean)) 
            eta <- tcrossprod(cbind(1,object$lv.mean), object$lv.coefs.mean[,1:(num.lv+1)]) 
        if(!is.null(object$X.coefs.mean)) 
            eta <- eta + tcrossprod(as.matrix(X), object$X.coefs.mean) 
            }

    if(!is.null(object$row.ids) && est == "median") {
        for(j in 1:p) 
            for(k in 1:ncol(object$row.ids)) 
                eta[,j] <- eta[,j] + object$row.coefs[[k]]$median[object$row.ids[,k]] 
        }
    if(!is.null(object$row.ids) && est == "mean") {
        for(j in 1:p) 
            for(k in 1:ncol(object$row.ids)) 
                eta[,j] <- eta[,j] + object$row.coefs[[k]]$mean[object$row.ids[,k]] 
        }
    
    index_multinom_cols <- which(object$family == "multinom")
    for(j in 1:p) {
        if(object$family[j] %in% c("binomial")) 
            fitted_out[,j] <- pnorm(eta[,j])
        if(object$family[j] %in% c("beta")) 
            fitted_out[,j] <- exp(eta[,j])/(1+exp(eta[,j]))
        if(object$family[j] %in% c("poisson","lnormal","negative.binomial","tweedie","exponential","gamma")) 
            fitted_out[,j] <- exp(eta[,j])
        if(object$family[j] == "normal") 
            fitted_out[,j] <- (eta[,j]) 
# 		if(object$family[j] == "multinom") {
# 			if(est == "median") { if(!is.null(object$X.multinom.coefs.median)) eta2 <- eta[,j] + as.matrix(X)%*%object$X.multinom.coefs.median[which(index_multinom_cols == j),,] }
# 			if(est == "mean") { if(!is.null(object$X.multinom.coefs.mean)) eta2 <- eta[,j] + as.matrix(X)%*%object$X.multinom.coefs.mean[which(index_multinom_cols == j),,] }
# 			get_probs <- exp(eta2)/apply(exp(eta2),1,sum)	
# 			fitted_multinom_probs[,j,] <- get_probs
# 			}

        if(object$family[j] == "ordinal") {
            if(est == "median")
                fitted_ordinal_probs[,j,] <- ordinal.conversion.spp(n = n, lv = object$lv.median, lv.coefs.j = object$lv.coefs.median[j,], num.lv = num.lv, row.coefs = object$row.coefs, row.ids = object$row.ids, X = X, X.coefs.j = object$X.coefs.median[j,], cutoffs = object$cutoffs.median, est = "median")
            if(est == "mean")
                fitted_ordinal_probs[,j,] <- ordinal.conversion.spp(n = n, lv = object$lv.mean, lv.coefs.j = object$lv.coefs.mean[j,], num.lv = num.lv, row.coefs = object$row.coefs, row.ids = object$row.ids, X = X, X.coefs.j = object$X.coefs.mean[j,], cutoffs = object$cutoffs.mean, est = "mean")
            fitted_out[,j] <- apply(fitted_ordinal_probs[,j,],1,which.max) ## get max predicted probability
            }
        }	

    return(list(ordinal.probs = fitted_ordinal_probs, out = fitted_out))
    }

	
## Calculates DIC based on the conditional log-likelihood
get.dic <- function(jagsfit) 
     { 
     jagsfit$BUGSoutput$DIC
     }
	

## Simple extraction of MCMC samples from fitted boral object
get.mcmcsamples <- function(object) {
    fit.mcmc <- object$jags.model$BUGSoutput
    if(is.null(fit.mcmc)) 
    stop("MCMC samples not found. Please use save.model = TRUE to save MCMC samples when using boral")
    fit.mcmc <- mcmc(fit.mcmc$sims.matrix, start = 1, thin = object$mcmc.control$n.thin) ## Thanks to Guilliaume Blanchet for the original formatting!

    return(fit.mcmc)
    }

     
