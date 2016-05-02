#####################___SCRIPT BODY___#######################################
#  Purpose: R script for grm() function for doing generalised linear modelling
#   with a gamma distributed response and canonical link function, using the 
#   iterative IWLS algorithm.
#
#       ========================================================
#       Structure of script: Each numbered section contains functionality for:
#       --------------------
#       1. Error handling
#       2. Initialisation of IWLS algorithm
#       3. IWLS algorithm loop for iterively estimating beta values
#       4. Computing results estimates and statistics
#       5. Computing diagnostic statistics
#       6. Creating results object, which contains:
            # y: The observed responses.
            # fitted: The fitted values.
            # betahat: The estimated regression coefficients.
            # sebeta: The standard errors of the estimated regression coefficients. 
            # cov.beta: The covariance matrix of the estimated regression coefficients.
            # phihat: The estimated dispersion parameter.
            # df.model: The degrees of freedom for the model.
            # df.residual: The residual degrees of freedom.
            # deviance: The deviance for the model.
            # residuals.raw: The residuals of the model.
            # The Pearson standardised residuals of the model.
            # The standardised deviance residuals of the model.
#       7. Create, format and output summary results in the console for inspection, including coefficient
#          estimates with accompanying t statistcs and p values
#       8. Produce diagnostic plots using diagnostic statistics, which are:
            # Residuals vs fitted
            # Normal Q-Q plot for the theoretical quantiles vs standardised deviance residuals
            # Scale location for sqrt(Std.deviance res) Vrs predicted valuess
            # Std. Pearson resid vs leverage

#       finally returns list object containing the results
#
# 
#####################___SCRIPT BODY___#######################################


grm<- function(y, X, startval= c(-0.1, -0.1)) {  
  
  # --- SECTION 1 ERROR HANDLING ---
  
  # if any of the the following error handling checks are violated (i.e. the condition in the if 
  # statements is TRUE, then the execution of the function will be stopped using stop() command 
  #  with a descriptive error message to aid with the users investigation)
 
  n <- length(y) # number of observations
  p <- ncol(X) # number of predictors for the model
  X = as.matrix(X) # convert X to matrix format, needed for matrix operations within the IWLS algorithm
  
  # check for missing data
  if (any(is.na(y))|| any(is.na(X))|| any(is.na(startval)))
          stop("Issue with input data: NAs detected. Please fix the missing data present in the inputs.")
  
  # check to ensure there are no dimensionality between the response and the predictor inputs
  if (n != nrow(X)){
    stop("Issue with input data: dimensionality mismatch detected. Please ensure that X and y have the same 
         number of rows.")
  }
  
  # check to ensure the problem is well posed, in the sense there are enough observations present for the number
  # of predictors, to allow for a useful model.
  if (p >=  n){
    stop("Issue with input data: ill posed problem detected. Please ensure the number of observations is greater 
         than the number of explanatory variables in X.")
  }
  
  # check formatting of y, ensure it is a vector with only positive values with quantitative response
  # (required by the definition of gamma density for the responses provided)
  if (is.vector(y)!=TRUE){
    stop("Issue with input data: multi-variate response detected. Please ensure y is univariate represented as a 
         vector of length the number of obervations in X")
  }
  if (is.factor(y)){
    stop("Issue with input data: categoical response deteced. GRM supports a quantitative response variable only")
  }
  if (any(y<0)){
    stop("Issue with input data: negative response deteced. Gamma density requires response to be positive")
  }
  
  # check to ensure the initial guesses for the linear predictor coefficients are negative, ensuring expectation will
  # be positive using. This constraint is required due to the choice of the canonical link.
  if (any(startval>0)){
    stop("Issue with input data: positive coefficient start values detected. 
          Because the canonical link is utilised in this implementeation of grm(), and because of the 
          restriction that the expectation of the response must be positive, please ensure the initial 
          guesses for the linear coefficient are negative valued.")
  }
  
  # --- SECTION 2 INITIALISATION FOR IWLS ALGORITHM ---
  
  betahat <- startval  # initial value for beta is provided by the function input startval  
  U <- 10   #Define a value for the score, U (this is just so that the test for convergence
  # on the next line doesn't fail on the first attempt because U has not been defined yet)
  
  # --- SECTION 3 IWLS ALGORITHM FOR ESTIMATING BETA ---
  
  # the while loop condition specifies the convergence criteria in terms of the score vector, that is
  # if all the values of U are sufficiently close to zero, take the current Beta values as our estimated
  # parameter value. Otherwise, continue to iterate, reducing U further and hence updating our Beta values. 
  # We utilise a modified version of the Newton-Raphson method where the Fisher Scoring method replaces 
  # the Hessian matrix by its expectation, which is the negative of the information matrix.
  
  while(all(abs(U) > 1e-6)) { 
    
    # get first linear predictor value
    eta <- as.vector(X%*%betahat) 
   
    # distribution specific parameter derivations:
    mu <- -1/eta   # compute the current estimate for the expectd value, using the canonical link function for gamma.
    
    if (any(mu <= 0)){ # check to ensure that mu stays positive throughout, required by the gamma distribution
      stop("The expectation of the response value has been set to 0 during the fisher 
           scoring method for estimation of linear predictor parameters, which is an 
           invalid value for the gamma density assumptions") 
    }
    
    V <- (mu^2) # compute the current variance function value           
    W <- (mu^2) # compute the current  weights
    z <- eta + (y-mu)*(1/(mu^2))# compute the current adjusted dependendent variate
    
    # generic IWLS steps, using the fisher scoring method:
    
    XW <- t(W*X) # calculation of X'W (uses elementwise multiplication, exploiting the fact that W will be 
                 # recycled to match the number of elements in X)     
    XWX <- solve(XW%*%X) # calculation of [X'WX]^-1    
    XWz <- XW%*%z # calculation of  X'Wz              
    U <- XW%*%(z-eta) # calculation of score of current iteration. This should convege towards zero.
    betahat <- XWX%*%XWz # update betahat, and go back if necessary
   
  }
  
  # When the while loop is broken on convergence, the latest value for betahat is taken as the final estimate,
  # and used to compute the results:
  
  # --- SECTION 4 COMPUTE RESULTS STATISTICS ---
  
  phihat <- sum(((y-mu)/mu)^2/(n-p)) # compute estimate for dispersion, based on the methods of moments
  betahat.se <- sqrt(diag(XWX)) # compute standard error of beta estimate     
  betahat.cov <- XWX *phihat # compute variance/covariance of beta estimate,
                             # using inverse of the information matrix.
  deviance <- 2*(-1 + (y/mu) + log(mu/y)) # compute deviance for each fitted value estimate
  D <- sum(deviance) # compute model deviance
  tstat = betahat/betahat.se # compute t statistic for beta estimate, to allow for p-value production,
                             # which is used in results to provide evidence for statistical signaficance of the beta 
                             # coeffcient estimates within the model.
  
  # --- SECTION 5 COMPUTE DIAGNOSTIC STATISTICS ---
  
  # The derivations and interpretations for these diagnostics refernece STG001 class notes and the 
  # report accompanying the script.
  
  residuals.raw = (y-mu) # compute the residuals between the expected value of the response estimates and 
                         # the actual response values, helping to show evidence for how well the raw data are fitted.
  
  residuals.pearson = (y-mu)/sqrt((V*phihat)) # compute pearson method residuals, standardising by standard 
                                              # deviation estimate of the response, making them comparable size.
  
  W.root = W^0.5 # compute weights required by the hat matrix, for use in standardising the pearson residuals
  hat.matrix = (W.root*X)%*%XWX%*%t(W.root*X) # compute hat matrix
  
  hat.diag = diag(hat.matrix) # compute diagonal elements of hat matrix, which represents the leverage of the ith observation
  hat.stddiag = (sqrt(1-hat.diag))
  residuals.pearson.std = residuals.pearson/hat.stddiag  # compute standardised pearson residuals using diagonal 
                                                         # elements of hat matrix, which results in variance of 1.
  
  residuals.deviance = sign(residuals.raw)*deviance # compute deviance residuals, which show how strongly the observation
                                                    # contributes to deviance.
  
  residuals.deviance.std = residuals.deviance/hat.stddiag # compute standardised deviance residuals, resulting in unified variance.
  cooks.stat = (1/p)*(hat.diag/(1-hat.diag))*(residuals.deviance.std^2) # compute cook statistic, which illustrates the effect 
                                                                        # the beta estimate of omitting observation i, which can
                                                                        # be used in investigating effects of leverage.
  
  # --- SECTION 6 CREATE RESULTS OBJECT ---
  
  results <- list(y=y,
                  fitted= mu,
                  betahat = betahat,
                  sebeta = betahat.se,
                  cov.beta = betahat.cov,
                  phihat = phihat,
                  df.model = p-1, 
                  df.residual = (n-1) - (p-1), 
                  deviance =  D,
                  residuals.raw = residuals.raw,
                  residuals.pearson.std = residuals.pearson.std,
                  residuals.deviance.std  = residuals.deviance.std)
  
  # --- SECTION 7 CONSOLE OUTPUT  ---
  # all output derived above rounded to 5 decimal places, except p value, which is formated to 3 significant digits
  # with scientific format.

  mle.table <- data.frame(Estimate=round(betahat,5), 
               S.E.= round(betahat.se,5) , 
               t = round(tstat,5),
               pvalue = formatC(pt(tstat, df = n-1),3)) # compute p value for beta coefficient estimates using 
                                                        # t statistics at 5% significance level.
  
  # print output to console 
  cat("\n Coefficients\n") # print output heading
  print( mle.table) # coefficients estimate stat table
  
  cat(" \n \nModel Details\n") # print output heading
  cat("\nLink function: Canonical") # in future could extend to have link function as an input, instead of hard coded
  # print model level stats, rounded to 2 decimal places
  cat(paste("\nDispersion estimate: ", round(phihat,5)))
  cat(paste("\nNumber of coefficients: ", round(p,5))) 
  cat(paste("\nResidual degrees of freedom: ", (n-1)-(p-1)))
  cat(paste("\nDeviance: ", round(D,5)))
      
  # --- SECTION 8 PRODUCE PLOTS  ---

  par(ask = TRUE) # command prompt to force user to cycle through plots one by one in the console
 
  # predictions vs residuals
  plot(mu, residuals.pearson, ylab = "Residuals", 
       xlab = "Predicted values", main = "Residuals vs Fitted")
  
  # Normal Q-Q: theoretical quantiles against standardised deviance residual
  qqnorm(residuals.deviance.std, ylab = "Std deviance resid", 
         xlab = "Theoretical Quantiles", main = "Normal Q-Q")

  # Scale-Location: sqrt(standardised deviance res) vs pred vals
  plot(mu, sqrt(residuals.deviance.std), ylab = "predicted values", 
       xlab = "sqrt(Std. deviance resid)", main = "Scale-Location")
  
  # Residuals Vs Leverage
  plot(hat.diag, residuals.pearson, ylab = "Std. Pearson resid", 
       xlab = "Leverage", main = "Std. Pearson Residuals vs Leverage")
  
  # return results list object for the function call
  return(results)  
  
}


