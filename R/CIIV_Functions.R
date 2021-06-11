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
#'@param robust Logical. If robust = TRUE, the linear model is robust to
#'  heteroskedasticity (default = TRUE).
#'@param firststage Logical. If firststage = TRUE, a first-stage thresholding is 
#'  implemented to select the relevant instrument variables (default = FALSE).
#'@param firsttuning A numeric scalar specifiying the threshold critical value for the
#'  first stage t-test (default = sqrt(2.01*log(pz))).
#'@return Valid instruments:
#'    Identities of the valid instrumental variables
#'  selected by the algorithm.
#'@return Number of Valid Instruments:
#'    The number of the selected valid instrumental variables.
#'@return Relevant instruments:
#'    Identities of the relevant instrumental variables
#'  selected by the first stage thresholding if firststage = TRUE.
#'@return Number of Relevant Instruments:
#'    The number of the selected relevant instrumental variables 
#'  if firststage = TRUE.
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
#'library(AER); library(sandwich)
#'# the MASS package is only needed to
#'run the working example
#'library(MASS)
#'#Generate data
#'n = 2000; L = 10; s = 3
#'pi = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(0, rep(1,(L-2)), 0)
#'epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
#'epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
#'Z = matrix(rnorm(n*L),n,L)
#'D = 0.5 + Z %*% gamma + epsilon[,1]
#'Y = -0.5 + Z %*% pi + D * beta + epsilon[,2]
#'result = CIIV(Y,D,Z,robust=FALSE, firststage=TRUE)
#'result
#' @export
CIIV <- function(Y, D, Z, X, alpha = 0.05, tuning = 0.1/log(length(Y)), robust = TRUE, firststage = FALSE, firsttuning = sqrt(2.01*log(ncol(Z)))) {
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
  stopifnot(is.logical(firststage), is.logical(robust));
  stopifnot(is.numeric(alpha), length(alpha) == 1,alpha <= 1,alpha >= 0);
  stopifnot(is.numeric(tuning), length(tuning) == 1);

  # Define Constants
  n <- length(Y);
  pz <- ncol(Z);

  # Preserve Data
  d = D;
  y = Y;
  z = Z;
  
  if (!missing(X)){
    X <- cbind(1, X);
    if(is.null(colnames(X))){colnames(X) = c("intercept", paste0('X', 1:ncol(X)))};
  }else{
    X <- matrix(1, n, 1);
    colnames(X) <- "intercept";
  };
  Covariates = cbind(D,X);
  Exogenous = cbind(Z,X);
  
  # Variable Names
  colnames(Covariates)[1] <- Dname;
  CovariatesNames <- colnames(Covariates);
  if(!is.null(colnames(Z))){InstrumentNames = colnames(Z)}else{InstrumentNames = paste0('Z', 1:ncol(Z))};
  
  # Centralization
  D <- qr.resid(qr(X), D);
  Z <- scale(qr.resid(qr(X), Z), center = FALSE, scale = TRUE);
  
  # First Stage
  if (firststage){
    
    lm_first <- lm(D ~ Z - 1);
    coef_first <- coef(lm_first);
    sd_first <- sqrt(diag(if(robust) vcovHC(lm_first, "HC0") else vcov(lm_first)));
    
    first_index <- as.numeric(abs(coef_first/sd_first) >= firsttuning);
    relevant_instruments <- which(first_index == 1);
    Nr_relevant <- length(relevant_instruments);
    
    if (sum(first_index) <= 1){
      warning("Less than two IVs are individually relevant, treat all IVs as strong");
      relevant_instruments_one <- which(first_index == 1);
      relevant_instruments <- c(1:pz);
    }
    
    if (length(relevant_instruments) < pz){
      X <- cbind(X, z[, -relevant_instruments]);
      Z <- as.matrix(z[, relevant_instruments], nrow = n, ncol = Nr_relevant);
      pz <- ncol(Z);
      
      D <- qr.resid(qr(X), d);
      Z <- scale(qr.resid(qr(X), Z), center = FALSE, scale = TRUE);
    }
    
  }else{
    relevant_instruments <- c(1:pz);
    Nr_relevant <- pz;
  };
  
  Y <- qr.resid(qr(X), Y); 

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
  Covmat_gamma <- if(robust) vcovHC(lm_Reduced, "HC0") else vcov(lm_Reduced);
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
    robust = robust,
    gamma_D = gamma_D,
    Covmat_gamma = Covmat_gamma,
    firststage = firststage,
    relevant_instruments = relevant_instruments,
    Nr_relevant = Nr_relevant
  );
}


#Selection by the Sargan/Hansen Downward Testing Procedure

#'Internal CIIV Functions
#'
#'These are not to be called by the user.
CIIV.TestingSelectionEstimation <- function(
  Y, D, Z, U, Covariates, Exogenous, y, z, CovariatesNames, InstrumentNames, betaIV, sdIV, alpha = 0.05, tuning = 0.1/log(length(Y)), robust = TRUE,
gamma_D, Covmat_gamma, firststage = FALSE, relevant_instruments, Nr_relevant) {

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
  
  colnames(Z) = paste0(relevant_instruments);
  colnames(z) = paste0(c(1:ncol(z)));
  
  # Break Points
  crit <- matrix(0, pz, pz);
  rownames(crit) = colnames(crit) <- paste0(relevant_instruments);
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
    ciIV[, 3] <- relevant_instruments;

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
      
      Zv <- Z[,!colnames(Z) %in% paste0(selVec)];
      if(robust) {
        sar_CI <- CIM.HansenJTest(Y, cbind(D, Zv), Z)[[1]];
      } else {
        res_CI <- resid(AER::ivreg(Y ~ cbind(D, Zv) - 1 | Z));
        sar_CI <- non_robust_sar_CI(res_CI, Z, U, n);
      }
      Psi <- max(crit[paste0(selVec), paste0(selVec)]) - epsi; # updated psi
      Psar <- pchisq(sar_CI, maxs - 1, lower.tail = FALSE);
      
    } else if(maxs >= 2) {  # Tie
      selCI <- cbind(ciIV[, 3], t(selCI));
      selVec <- matrix(0, maxs, ncol(selCI) - 1);
      Psar <- Psi <- rep(0, ncol(selCI) - 1);

      for (i in 1:(ncol(selCI) - 1)) {
        selVec[, i] <- sort(selCI[selCI[, i+1] == 1, 1]);

        if(robust){
          sar_CI <- CIM.HansenJTest(Y, cbind(D, Z[,!colnames(Z) %in% paste0(selVec[, i])]), Z)[[1]];
          #sar_CI <- CIM.HansenJTest(Y, cbind(D, Z[, -selVec[, i]]), Z)[[1]];
        } else{
          res_CI <- resid(AER::ivreg(Y ~ D + Z[,!colnames(Z) %in% paste0(selVec[, i])] - 1 | Z));
          sar_CI <- non_robust_sar_CI(res_CI, Z, U, n);
        }

        Psar[i] <- pchisq(sar_CI, maxs - 1, lower.tail = FALSE);
        Psi[i] <- max(crit[paste0(selVec[, i]), paste0(selVec[, i])]);
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

  # Post-Selection Estimation
  if (Nr_valid == 0) {
    print("None of the instruments is selected as valid, do OLS.")

    regressor_CIM_temp = regressor_CIM <- cbind(Covariates, z);
    length_regressor <- ncol(regressor_CIM)-ncol(z);
    Coefficients_CIM <- qr.coef(qr(regressor_CIM), y)[1:length_regressor];
    res_CIM <- qr.resid(qr(regressor_CIM), y);
  } else {
    # At least one of the IVs is selected as valid.
    
    z_invalid <- matrix(z[,!colnames(z) %in% paste0(selVec)], ncol = (ncol(z) - Nr_valid), nrow = n);
    
    regressor_CIM_temp <- cbind(Covariates, z_invalid);
    regressor_CIM <- cbind(fitted(lm(Covariates[,1] ~ Exogenous)), Covariates[,-1], z_invalid)
    length_regressor <- ncol(regressor_CIM_temp) - ncol(z_invalid);
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
    ci_CIM <- c(
      Coefficients_CIM[1] - qnorm(1-alpha/2) * sd_CIM[1],
      Coefficients_CIM[1] + qnorm(1-alpha/2) * sd_CIM[1]
    );
    if (Nr_valid == 0){
      Coefficients_CIM_GMM <- NA;
      sd_CIM_GMM <- NA;
      ci_CIM_GMM <- NA;
    }else{
      sd_CIM_GMM <- (CIM.HansenJTest(y, regressor_CIM_temp, Exogenous)[[3]])[1:length_regressor];
      ci_CIM_GMM <- c(
        Coefficients_CIM_GMM[1] - qnorm(1-alpha/2) * sd_CIM_GMM[1],
        Coefficients_CIM_GMM[1] + qnorm(1-alpha/2) * sd_CIM_GMM[1]
      );
    }
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
      if(is.null(selVec)){
        Valid_Instruments = paste0("None instruments selected as valid.");
        }else{
          Valid_Instruments = InstrumentNames[c(1:length(z))[selVec]];
      },
      Valid_Instruments = Valid_Instruments,
      Nr_valid = Nr_valid,
      if(firststage){
        if(Nr_relevant == 0){
          Relevant_Instruments = paste0("None instruments selected as relevant.")
        }else{
          if (Nr_relevant == 1){
            Relevant_Instruments = InstrumentNames[c(1:length(z))[relevant_instruments_one]];
          }else{
            Relevant_Instruments = InstrumentNames[c(1:length(z))[relevant_instruments]];
          }
        };
      }else{
        Relevant_Instruments = paste0("No first stage selection.");
        Nr_relevant = paste0("No first stage selection.");
      },
      Relevant_Instruments = Relevant_Instruments,
      Nr_relevant = Nr_relevant,
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
    if(is.null(selVec)){
      Valid_Instruments = paste0("None instruments selected as valid.");
    }else{
      Valid_Instruments = InstrumentNames[c(1:ncol(z))[selVec]];
    },
    Valid_Instruments = Valid_Instruments,
    Nr_valid = Nr_valid,
    if(firststage){
      if(Nr_relevant == 0){
        Relevant_Instruments = paste0("None instruments selected as relevant.")
      }else{
        if (Nr_relevant == 1){
          Relevant_Instruments = InstrumentNames[c(1:length(z))[relevant_instruments_one]];
        }else{
          Relevant_Instruments = InstrumentNames[c(1:length(z))[relevant_instruments]];
        }
      };
    }else{
      Relevant_Instruments = paste0("No first stage selection.");
      Nr_relevant = paste0("No first stage selection.");
    },
    Relevant_Instruments = Relevant_Instruments,
    Nr_relevant = Nr_relevant,
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
  cat("\nValid Instruments:\n", object$Valid_Instruments, "\n","\nNumber of Valid Instruments:\n", object$Nr_valid,"\n");
  cat("\nRelevant Instruments:\n", object$Relevant_Instruments, "\n","\nNumber of Relevant Instruments:\n", object$Nr_relevant,"\n");
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
