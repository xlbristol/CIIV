# CIIV
CIIV is a set of functions to implement the Confidence Interval method, which is to select the valid instrumental variables from a set of potential instruments that may contain invalid ones in linear models. It also provides post-selection IV estimation for the causal effect of the treatment/exposure on the outcome. The method in these R functions is based on the paper *The Confidence Interval Method for Selecting Valid Instrumental Variables (Windmeijer, Liang, Hartwig, Bowden, 2020)*.

## Installation
source("https://raw.githubusercontent.com/xlbristol/CIIV/main/R/CIIV_Functions.R")

## Input and Output of CIIV Functions
* **Input**
  * **Y:** A numeric vector of outcomes.
  * **D:** A numeric vector of exposures/treatments.
  * **Z:** A numeric matrix of instrumental variables, with each column referring to one instrument.
  * **X:** An optional numeric matrix of exogenous explanatory variables, with each column referring to one variable.
  * **alpha:** A numeric scalar between 0 and 1 specifying the significance level for the confidence interval for the causal effect estimate (default = 0.05).
  * **tuning:** A numeric scalar specifiying the threshold p-value for the Saran/Hansen-J test (default = 0.1/log(n)).
  * **intercept:** Logical. If intercept = TRUE, an intercept term is included in the linear model (default = TRUE).
  * **robust:** Logical. If robust = TRUE, the linear model is robust to heteroskedasticity.
  * **firststage:** Logical. If firststage = TRUE, a first-stage thresholding is implemented to select the relevant instrument variables (default = FALSE).
* **Output**
  * **Valid Instruments:** Identities of the valid instrumental variables selected by the algorithm.
  * **Number of Valid Instruments:** The number of the selected valid instrumental variables.
  * **Relevant Instruments:** Identities of the relevant instrumental variables selected by the first stage thresholding if firststage = TRUE.
  * **Number of Relevant Instruments:** The number of the selected relevant instrumental variables if firststage = TRUE.
  * **Coefficients:** The matrix for the post-selection IV estimation results for the coefficients of the exposure/treatment variable and exogenous
  explanatory variables using the selected valid instruments as IV and controlling for the selected invalid instruments. The first two columns are
  2SLS estimates and their standard errors. If robust = TRUE, the second column is heteroskedasticity-robust stand errors. The two-step GMM estimates
  and their standard errors are also reported. If intercept = TRUE, the estimates for the intercept are reported in the last row.
  * **Confidence Interval:** The confidence interval for the 2SLS estimates for the coefficient of the exposure/treatment variable with significance level
  specified by alpha (default = 0.05). If robust = TRUE, the confidence interval for the two-step GMM estimate is also reported.
  * **p-value of Sargan:** The p-value for the Sargan overidentifying test for the selected valid instruments. If robust = TRUE, the p-value for the
  Hansen-J test is reported.
  
  ## Working Example
```
# Load CIIV Functions
source("https://raw.githubusercontent.com/xlbristol/CIIV/main/R/CIIV_Functions.R")

# Load Dependent R Packages
library(AER)
library(sandwich)
library(MASS) ## The MASS package is only needed to run the working example.

# Generate Data
n = 2000; L = 10; s = 3
pi = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(0, rep(1,(L-2)), 0)
epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
Z = matrix(rnorm(n*L),n,L)
D = 0.5 + Z %*% gamma + epsilon[,1]
Y = -0.5 + Z %*% pi + D * beta + epsilon[,2]

result = CIIV(Y,D,Z,robust=FALSE, firststage=TRUE)
result
```
