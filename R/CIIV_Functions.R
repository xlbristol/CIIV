#Depends: AER, sandwich
library(AER);
library(sandwich);

#'The Confidence Interval Method for Selecting Valid Instrumental Variables.
#'
#'Implements the confidence interval (CI) downward testing procedure based on
#'the Sargan/Hansen-J test for overidentifying restrictions. The algorithm
#'selects the valid instrumental variables from a set of potential instruments
#'that may contain invalid ones. It also provides post-selection IV estimation
#'results using the selected valid instruments as IV and controlling for the
#'selected invalid instruments as explanatory variables.
#'
#'@param Y A numeric vector of outcomes.
#'@param D A numeric vector of exposures/treatments.
#'@param Z A numeric matrix of instrumental variables, with each column
#'  referring to one instrument.
#'@param X An optional numeric matrix of exogenous explanatory variables, with
#'  each column referring to one variable.
#'@param alpha A numeric scalar between 0 and 1 specifying the significance
#'  level for the confidence interval for the causal effect estimate (default =
#'  0.05).
#'@param tuning A numeric scalar specifiying the threshold p-value for the
#'  Saran/Hansen-J test (default = 0.1/log(n)).
#'@param intercept Logical. If intercept = TRUE, an intercept term is included
#'  in the linear model (default = TRUE).
#'@param robust Logical. If robust = TRUE, the linear model is robust to
#'  heteroskedasticity.
#'@return Invalid instruments:
#'    Identities of the invalid instrumental variables
#'  selected by the algorithm.
#'@return Number of Invalid Instruments:
#'    The number of the selected invalid
#'  instrumental variables.
#'@return Coefficients:
#'    The matrix for the post-selection IV estimation results
#'  for the coefficients of the exposure/treatment variable and exogenous
#'  explanatory variables using the selected valid instruments as IV and
#'  controlling for the selected invalid instruments. The first two columns are
#'  2SLS estimates and their standard errors. If robust = TRUE, the second
#'  column is heteroskedasticity-robust stand errors. The two-step GMM estimates
#'  and their standard errors are also reported. If intercept = TRUE, the
#'  estimates for the intercept are reported in the last row.
#'@return Confidence Interval:
#'    The confidence interval for the 2SLS estimates for
#'  the coefficient of the exposure/treatment variable with significance level
#'  specified by alpha (default = 0.05). If robust = TRUE, the confidence
#'  interval for the two-step GMM estimate is also reported.
#'@return p-value of Sargan:
#'    The p-value for the Sargan overidentifying test for
#'  the selected valid instruments. If robust = TRUE, the p-value for the
#'  Hansen-J test is reported.
#'@examples
#'library(AER) library(sandwich)
#'# the MASS package is only needed to
#'run the working example
#'library(MASS)
#'#Generate data
#'n = 2000; L = 10; s = 3
#'pi = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,L))
#'epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
#'epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
#'Z = matrix(rnorm(n*L),n,L)
#'D = 0.5 + Z %*% gamma + epsilon[,1]
#'Y = -0.5 + Z %*% pi + D * beta + epsilon[,2]
#'result = CIIV(Y,D,Z,robust=FALSE)
#'result
#' @export
CIIV <- function(Y, D, Z, X, alpha = 0.05, tuning = 0.1/log(length(Y)), intercept = TRUE, robust = TRUE) {
  #Check Input
  #Check Data

  if(is.data.frame(Y)){Y <- as.matrix(Y)}
  if(!is.vector(Y) && !is.matrix(Y) | !is.numeric(Y) | ncol(as.matrix(Y))!=1)
    stop("Y must be a numeric vector.");
  Y = as.numeric(Y);

  if(!is.null(colnames(D))){Dname = colnames(D)}else{Dname = "D"};
  if(is.data.frame(D)){D <- as.matrix(D)}
  if(!is.vector(D) && !is.matrix(D) | !is.numeric(D) | ncol(as.matrix(D))!=1)
    stop("D must be a numeric vector.");
  D = as.numeric(D);

  if(is.data.frame(Z)){Z <- as.matrix(Z)}
  if(!is.matrix(Z) | !is.numeric(Z))
    stop("Z must be a numerical matrix.");
  if(ncol(Z)<2)
    stop("The number of instruments must be greater than 1.");
  if( ncol(Z) > nrow(Z))
    stop("The number of instruments must be smaller than the sample size.");

  if(!missing(X) && is.data.frame(X)){X <- as.matrix(X)}
  if(!missing(X) && (!is.matrix(X) | !is.numeric(X)))
    stop("X must be a numerical matrix.");

  stopifnot(length(Y) == length(D),length(Y) == nrow(Z));

  #Other Arguments
  stopifnot(is.logical(intercept), is.logical(robust));
  stopifnot(is.numeric(alpha), length(alpha) == 1,alpha <= 1,alpha >= 0);
  stopifnot(is.numeric(tuning), length(tuning) == 1);

  # Define Constants
  n <- length(Y);
  pz <- ncol(Z);

  # Preserve Data and Variable Names
  if(!missing(X)){
    Covariates = cbind(D,X);
    Exogenous = cbind(Z,X);
    if(!is.null(colnames(X))){colnames(Covariates)[2:(ncol(X)+1)] = colnames(X)}else{colnames(Covariates)[2:(ncol(X)+1)] = paste0('X', 1:ncol(X))};
  }else{
    Covariates = cbind(D);
    Exogenous = cbind(Z);
  }
  colnames(Covariates)[1] <- Dname;
  if(intercept){CovariatesNames <- c(colnames(Covariates), "intercept")}else{CovariatesNames <- colnames(Covariates)}
  if(!is.null(colnames(Z))){InstrumentNames = colnames(Z)}else{InstrumentNames = paste0('Z', 1:ncol(Z))};
  y = Y;
  z = Z;

  # Centralize
  if(!missing(X)) {
    if(intercept){X <- cbind(1, X)}
    Y <- qr.resid(qr(X), Y);
    D <- qr.resid(qr(X), D);
    Z <- qr.resid(qr(X), Z);
  } else {
    Y <- Y - mean(Y);
    D <- D - mean(D);
    Z <- scale(Z, center = TRUE, scale = TRUE);
  }

  # OLS Estimation for the IV-Specific Estimates and Standard Error
  lm_Reduced <- lm(cbind(Y, D) ~ Z - 1);

  # IV-Specific Estimates
  gamma_Y <- coef(lm_Reduced)[, 1];
  gamma_D <- coef(lm_Reduced)[, 2];
  betaIV <- gamma_Y / gamma_D;

  # Variances by Delta Method
  U <- solve(crossprod(Z));
  Jacobian_betaIV <- cbind(
    diag(c(1 / gamma_D)), diag(c(-gamma_Y / gamma_D^2))
  );
  #Covmat_gamma <- if(robust) vcovHC(lm_Reduced, "HC0") else vcov(lm_Reduced);
  Covmat_gamma <- if(robust)  vcovHAC(lm_Reduced) else vcov(lm_Reduced);
  sdIV <- sqrt(diag(Jacobian_betaIV %*% Covmat_gamma %*% t(Jacobian_betaIV)));

  # CI Downward Testing Sargen/Hansen-J Procedure and Post-Selection Estimation
  CIIV.TestingSelectionEstimation(
    Y, D, Z, U,
    Covariates,
    Exogenous,
    y, z,
    CovariatesNames,
    InstrumentNames,
    betaIV = betaIV,
    sdIV = sdIV,
    alpha = alpha,
    tuning = tuning,
    intercept = intercept,
    robust = robust
  );
}


#Selection by the Sargan/Hansen Downward Testing Procedure

#'Internal CIIV Functions
#'
#'These are not to be called by the user.
CIIV.TestingSelectionEstimation <- function(
  Y, D, Z, U, Covariates, Exogenous, y, z, CovariatesNames, InstrumentNames, betaIV, sdIV, alpha = 0.05, tuning = 0.1/log(length(Y)), intercept = TRUE, robust = TRUE
) {

  # Function for Sargan Test
  non_robust_sar_CI <- function(res_CI, Z, U, n) {
    (t(res_CI) %*% Z %*% U %*% t(Z) %*% res_CI) /
      (t(res_CI) %*% res_CI / n)
  }

  # Function for Two Step GMM and Hansen-J Test
  CIM.HansenJTest <- function(Y,X,Z){
    res_FirstStep <- residuals(AER::ivreg(Y ~ X - 1 | Z));

    Weight_SecondStep <- crossprod(res_FirstStep * Z);

    Coef_SecondStep <- solve(
      t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
    ) %*% t(X) %*% Z %*% solve(Weight_SecondStep) %*% t(Z) %*% Y;

    res_SecondStep <- as.vector(Y - X %*% Coef_SecondStep);

    sd_SecondStep <- sqrt(diag(solve(
      t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
    ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)%*%crossprod(res_SecondStep * Z)%*%t(
      solve(
        t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
      ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)
    )));

    HansenJ_Stat <- t(res_SecondStep) %*% Z %*% solve(Weight_SecondStep) %*%
      t(Z) %*% res_SecondStep;

    list(HansenJ_Stat, Coef_SecondStep,sd_SecondStep)
  }

  # Define Constants
  n <- length(Y);
  pz <- ncol(Z);

  # Break Points
  crit <- matrix(0, pz, pz);
  for (i in 1:pz) {
    for (j in 1:pz) {
      crit[i, j] <- abs(betaIV[i] - betaIV[j]) / (sdIV[i] + sdIV[j]);
    }
  }

  # The Downward Testing Procedure
  maxs <- pz;
  Psar <- 0;
  epsi <- 10^(-7);
  Psi <- max(crit) + epsi;

  while (maxs > 0 && Psar < tuning) {
    # The confidence intervals
    ciIV <- matrix(0, pz, 3);
    ciIV[, 1] <- betaIV - Psi * sdIV;
    ciIV[, 2] <- betaIV + Psi * sdIV;
    ciIV[, 3] <- 1:pz;

    # Ordering the Confidence Intervals by the Lower Ends
    ciIV <- ciIV[order(ciIV[, 1]), ];

    # Check Overlap
    grid_CI <- matrix(0, pz, pz)
    for (i in 2:pz) {
      for (j in 1:i) {
        grid_CI[i, j] <- as.numeric(ciIV[i, 1] <= ciIV[j, 2]);
      }
    }

    sumCI <- apply(grid_CI, 1, sum);
    selCI <- as.matrix(grid_CI[which(sumCI == max(sumCI)), ]);
    maxs <- max(sumCI);

    #Selection
    if(maxs >= 2 && length(selCI) < pz + 1) {  # No tie
      selCI <- cbind(ciIV[, 3], selCI);
      selVec <- sort(selCI[selCI[, 2] == 1, 1]);

      if(maxs == pz) {  # All valid
        if(robust) {  # Do GMM and Hansen-J test
          sar_CI <- CIM.HansenJTest(Y, D, Z)[[1]];
        } else {
          res_CI <- resid(AER::ivreg(Y ~ D - 1 | Z));
          sar_CI <- non_robust_sar_CI(res_CI, Z, U, n);
        }
        Psi <- max(crit) - epsi; # updated psi
      } else {
        Zv <- Z[, -selVec];
        if(robust) {
          sar_CI <- CIM.HansenJTest(Y, cbind(D, Zv), Z)[[1]];
        } else {
          res_CI <- resid(AER::ivreg(Y ~ D + Zv - 1 | Z));
          sar_CI <- non_robust_sar_CI(res_CI, Z, U, n);
        }
        Psi <- max(crit[selVec, selVec]) - epsi; # updated psi
      }

      Psar <- pchisq(sar_CI, maxs - 1, lower.tail = FALSE);
    } else if(maxs >= 2) {  # Tie
      selCI <- cbind(ciIV[, 3], t(selCI));
      selVec <- matrix(0, maxs, ncol(selCI) - 1);
      Psar <- Psi <- rep(0, ncol(selCI) - 1);

      for (i in 1:(ncol(selCI) - 1)) {
        selVec[, i] <- sort(selCI[selCI[, i+1] == 1, 1]);

        if(robust){
          sar_CI <- CIM.HansenJTest(Y, cbind(D, Z[, -selVec[, i]]), Z)[[1]];
        } else{
          res_CI <- resid(AER::ivreg(Y ~ D + Z[, -selVec[, i]] - 1 | Z));
          sar_CI <- non_robust_sar_CI(res_CI, Z, U, n);
        }

        Psar[i] <- pchisq(sar_CI, maxs - 1, lower.tail = FALSE);
        Psi[i] <- max(crit[selVec[, i], selVec[, i]]);
      }

      index.tie <- match(max(Psar), Psar);
      selVec <- selVec[,index.tie];
      Psar <- Psar[index.tie];
      Psi <- min(Psi) - epsi;  # updated psi
    } else {  # None Valid
      maxs <- 0;
      Psar <- 0;
      selVec <- NULL;
    }
  }

  Nr_valid <- length(selVec);
  Nr_invalid <- pz - Nr_valid;

  # Post-Selection Estimation
  if (Nr_valid == 0) {
    print("None of the instruments is selected as valid, do OLS.")

    if(intercept){
      regressor_CIM_temp = regressor_CIM <- cbind(Covariates, 1, z);
    }else{
      regressor_CIM_temp = regressor_CIM <- cbind(Covariates, z);
    }
    length_regressor <- ncol(regressor_CIM)-ncol(z);
    Coefficients_CIM <- qr.coef(qr(regressor_CIM), y)[1:length_regressor];
    res_CIM <- qr.resid(qr(regressor_CIM), y);
  } else {
    # At least one of the IVs is selected as valid.

    if(intercept){
      regressor_CIM_temp <- cbind(Covariates, 1, z[, -selVec]);
      regressor_CIM <- cbind(fitted(lm(Covariates[,1] ~ Exogenous)), Covariates[,-1], 1, z[, -selVec]);
    }else{
      regressor_CIM_temp <- cbind(Covariates, z[, -selVec]);
      regressor_CIM <- cbind(fitted(lm(Covariates[,1] ~ Exogenous)), Covariates[,-1], z[, -selVec])
    }
    length_regressor <- ncol(regressor_CIM_temp) - ncol(z[, -selVec]);
    iv_CIM <- AER::ivreg(y ~ regressor_CIM_temp - 1 | Exogenous);
    Coefficients_CIM <- coef(iv_CIM)[1:length_regressor];
    if(robust){
      Coefficients_CIM_GMM <- (CIM.HansenJTest(y, regressor_CIM_temp, Exogenous)[[2]])[1:length_regressor];
      }
    res_CIM <- resid(iv_CIM);
  }

  # Standard Error and Confidence Interval
  if (robust) {
    sd_CIM <- sqrt(diag(
      solve(crossprod(regressor_CIM)) %*%
      crossprod(res_CIM * regressor_CIM) %*%
      solve(crossprod(regressor_CIM))
    ))[1:length_regressor];
    sd_CIM_GMM = (CIM.HansenJTest(y, regressor_CIM_temp, Exogenous)[[3]])[1:length_regressor];
    ci_CIM <- c(
      Coefficients_CIM[1] - qnorm(1-alpha/2) * sd_CIM[1],
      Coefficients_CIM[1] + qnorm(1-alpha/2) * sd_CIM[1]
    );
    ci_CIM_GMM <- c(
      Coefficients_CIM_GMM[1] - qnorm(1-alpha/2) * sd_CIM_GMM[1],
      Coefficients_CIM_GMM[1] + qnorm(1-alpha/2) * sd_CIM_GMM[1]
    );
  } else {
    sd_CIM <- sqrt(diag(
      mean(res_CIM^2) * solve(crossprod(regressor_CIM))
    ))[1:length_regressor];
    ci_CIM <- c(
      Coefficients_CIM[1] - qnorm(1-alpha/2) * sd_CIM[1],
      Coefficients_CIM[1] + qnorm(1-alpha/2) * sd_CIM[1]
    );
  }

  # Results
  if(robust){
    object <- list(
      robust = robust,
      if(is.null(selVec)){Invalid_Instruments = InstrumentNames[1:pz]}else{Invalid_Instruments = InstrumentNames[c(1:pz)[-selVec]]},
      Invalid_Instruments = Invalid_Instruments,
      Nr_invalid = Nr_invalid,
      Covariate_Names = CovariatesNames,
      Coefficients_CIM = Coefficients_CIM,
      sd_CIM = sd_CIM,
      Coefficients_CIM_GMM = Coefficients_CIM_GMM,
      sd_CIM_GMM = sd_CIM_GMM,
      ci_CIM = ci_CIM,
      ci_CIM_GMM = ci_CIM_GMM,
      HansenJ_CIM = Psar
    )
  }else{
    object <- list(
    robust = robust,
    if(is.null(selVec)){Invalid_Instruments = InstrumentNames[1:pz]}else{Invalid_Instruments = InstrumentNames[c(1:pz)[-selVec]]},
    Invalid_Instruments = Invalid_Instruments,
    Nr_invalid = Nr_invalid,
    Covariate_Names = CovariatesNames,
    Coefficients_CIM = Coefficients_CIM,
    sd_CIM = sd_CIM,
    ci_CIM = ci_CIM,
    Sargan_CIM = Psar
    )
  }

  class(object) <- "CIIV"

  object
}

#'Internal CIIV Functions
#'
#'These are not to be called by the user.
print.CIIV <- function(object,robust = object$robust,...){
  cat("\nInvalid Instruments:\n", object$Invalid_Instruments, "\n","\nNumber of Invalid Instruments:\n", object$Nr_invalid,"\n");
  cat("\nCoefficients:\n")

  names(object$Coefficients_CIM) = names(object$sd_CIM) = names(object$ci_CIM) = NULL;
  if(robust){
    names(object$Coefficients_CIM_GMM) = names(object$sd_CIM_GMM) = names(object$ci_CIM_GMM) = NULL;
    coef_cim <- cbind(object$Coefficients_CIM, object$sd_CIM, object$Coefficients_CIM_GMM,object$sd_CIM_GMM);
    colnames(coef_cim) = c("2SLS Estimate", "2SLS Std. Error", "GMM Estimate", "GMM Std. Error");
    rownames(coef_cim) = object$Covariate_Names;
    print(coef_cim, quote = FALSE);

    cat("\nConfidence Interval 2SLS: [", object$ci_CIM[1], ",", object$ci_CIM[2], "]", "\n", sep = '');
    cat("\nConfidence Interval GMM: [", object$ci_CIM_GMM[1], ",", object$ci_CIM_GMM[2], "]", "\n", sep = '');

    cat("\np-value of Hansen-J: ", object$HansenJ_CIM, sep = '');
  }else{
    coef_cim <- cbind(object$Coefficients_CIM, object$sd_CIM);
    colnames(coef_cim) = c("2SLS Estimate", "Std. Error");
    rownames(coef_cim) = object$Covariate_Names;
    print(coef_cim, quote = FALSE);

    cat("\nConfidence Interval: [", object$ci_CIM[1], ",", object$ci_CIM[2], "]", "\n", sep = '');

    cat("\np-value of Sargan: ", object$Sargan_CIM, sep = '');
  }
}
