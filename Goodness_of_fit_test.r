# This version: 22 May 2018
# This code implements several tests studied in 
# K. Hayakawa (2018?) "Corrected Goodness of Fit Test in Covariance Structure Analysis"
# forthcoming in Psychological Methods
# It also provides fit indices based on T_RLS.  
# If you find any bugs, please let me know.

GFT <- function(fit1, fit2, fit3, X) 
{
S_N <- var(X)
N   <- nrow(X)
p   <- ncol(X)

# Conventinal goodness of fit test (T_ML)
T_ML  <- lavInspect(fit1, "test")[[1]]$stat
df_ML <- lavInspect(fit1, "test")[[1]]$df
p_ML  <- lavInspect(fit1, "test")[[1]]$pvalue


# Corrected goodness of fit test (T_RLS)
Sigma_hat0 <- fitted(fit1)[1]
Sigma_hat  <- matrix(unlist(Sigma_hat0), ncol = p, byrow = TRUE)
R <- (S_N-Sigma_hat)%*%solve(Sigma_hat)
T_RLS  <- (N-1)/2*sum(diag( R%*%R ))
df_RLS <- lavInspect(fit1, "test")[[1]]$df
p_RLS  <- 1 - pchisq(T_RLS,df_RLS)


# Sattora-Bentler scaled tests based on T_ML
T_SB1_ML  <- (N-1)/N*lavInspect(fit2, "test")[[2]]$stat
df_SB1_ML <- lavInspect(fit2, "test")[[2]]$df
p_SB1_ML  <- lavInspect(fit2, "test")[[2]]$pvalue


# Sattora-Bentler scaled tests based on T_RLS
scale_factor_SB1 <- N/(N-1)*lavInspect(fit2, "test")[[2]]$scaling.factor
T_SB1_RLS  <- T_RLS/scale_factor_SB1
df_SB1_RLS <- df_SB1_ML 
p_SB1_RLS  <- 1 - pchisq(T_SB1_RLS,df_SB1_RLS)


# Sattora-Bentler adjusted tests based on T_ML
T_SB2_ML  <- (N-1)/N*lavInspect(fit3, "test")[[2]]$stat
df_SB2_ML <- lavInspect(fit3, "test")[[2]]$df
p_SB2_ML  <- lavInspect(fit3, "test")[[2]]$pvalue


# Sattora-Bentler adjusted tests based on T_RLS
scale_factor_SB2 <- N/(N-1)*lavInspect(fit3, "test")[[2]]$scaling.factor
T_SB2_RLS  <- T_RLS/scale_factor_SB2
df_SB2_RLS <- df_SB2_ML 
p_SB2_RLS  <- 1 - pchisq(T_SB2_RLS, df_SB2_RLS)


# Sattora-Benlter adjusted test with improved scaling factor
Xbar <- colMeans(X)
C_N  <- NULL
for (i in 1:N) {
   S_i <- (X[i,]-Xbar)%*%t(X[i,]-Xbar)
   C_i <- lav_matrix_vech( S_i )
   C_N <- rbind(C_N, C_i)        
 }  

U  <- lavInspect(fit1, "UfromUGamma")
Q1 <- svd(U, LINPACK = FALSE)$u
W1 <- svd(U, LINPACK = FALSE)$v
D1 <- diag(sqrt(svd(U, LINPACK = FALSE)$d))
Usqrt <- Q1 %*% D1 %*% t(W1)

UX <- (Usqrt)%*%t(C_N)    
Y  <- UX - rowMeans(UX)%*%matrix(1,1,N)

V <- Y%*%t(Y)
M <- t(Y)%*%Y
D <- diag(diag(M))

trUG  <- sum(diag(V))/(N-1)
trUG2 <- sum(diag(V%*%V))/(N-1)^2

if (nrow(V) < nrow(M)) {
    trUG2_corr <- ( (N-2)*(N-1)*sum(diag(V%*%V)) - N*(N-1)*sum(diag(D%*%D)) + ( sum(diag(V)) )^2  )/(N*(N-1)*(N-2)*(N-3))
 } else {
    trUG2_corr <- ( (N-2)*(N-1)*sum(diag(M%*%M)) - N*(N-1)*sum(diag(D%*%D)) + ( sum(diag(D)) )^2  )/(N*(N-1)*(N-2)*(N-3))
 }


scale_factor_SB2_c <- trUG2_corr/trUG
T_SB2_ML_c  <- T_ML/scale_factor_SB2_c
T_SB2_RLS_c <- T_RLS/scale_factor_SB2_c

df_SB2_ML_c  <- trUG^2/trUG2_corr
df_SB2_RLS_c <- trUG^2/trUG2_corr 
p_SB2_ML_c   <- 1 - pchisq(T_SB2_ML_c, df_SB2_ML_c)
p_SB2_RLS_c  <- 1 - pchisq(T_SB2_RLS_c,df_SB2_RLS_c)



cat(" T_ML : ", "classical goodness of fit test", "(option: estimator = ML)", "\n",
"Test Statistic       = ", T_ML,  "\n", 
"Degrees of freedom   = ", df_ML, "\n", 
"P-value (Chi-square) = ", p_ML,  "\n",  "\n" )

cat(" T_RLS : ", "corrected goodness of fit test", "(option: estimator = ML)", "\n",
"Test Statistic       = ", T_RLS,  "\n", 
"Degrees of freedom   = ", df_RLS, "\n", 
"P-value (Chi-square) = ", p_RLS,  "\n",  "\n" )

cat("T_SB1_ML : ", "Sattora-Bentler scale test based on T_ML", "(option: estimator = MLM)", "\n", 
"Test Statistic       = ", T_SB1_ML,  "\n", 
"Degrees of freedom   = ", df_SB1_ML, "\n", 
"P-value (Chi-square) = ", p_SB1_ML,  "\n",
"Scaling  factor      = ", scale_factor_SB1,  "\n",  "\n"  )


cat("T_SB1_RLS : ", "Sattora-Bentler scale test based on T_RLS", "(option: estimator = MLM)", "\n",
"Test Statistic       = ", T_SB1_RLS,  "\n", 
"Degrees of freedom   = ", df_SB1_RLS, "\n", 
"P-value (Chi-square) = ", p_SB1_RLS,  "\n",
"Scaling  factor      = ", scale_factor_SB1,  "\n",  "\n"  )


cat("T_SB2_ML : ", "Sattora-Bentler adjusted test based on T_ML", "(option: estimator = MLMVS)", "\n", 
"Test Statistic       = ", T_SB2_ML,  "\n", 
"Degrees of freedom   = ", df_SB2_ML, "\n", 
"P-value (Chi-square) = ", p_SB2_ML,  "\n",
"Scaling  factor      = ", scale_factor_SB2,  "\n",  "\n"  )


cat("T_SB2_RLS : ", "Sattora-Bentler adjusted test based on T_RLS", "(option: estimator = MLMVS)", "\n", 
"Test Statistic       = ", T_SB2_RLS,  "\n", 
"Degrees of freedom   = ", df_SB2_RLS, "\n", 
"P-value (Chi-square) = ", p_SB2_RLS,  "\n",
"Scaling  factor      = ", scale_factor_SB2,  "\n",  "\n"  )


cat("T_SB2_ML_c : ", "Sattora-Bentler adjusted test based on T_ML with improved scaling factor", "(option: estimator = MLMVS)", "\n", 
"Test Statistic       = ", T_SB2_ML_c,  "\n", 
"Degrees of freedom   = ", df_SB2_ML_c, "\n", 
"P-value (Chi-square) = ", p_SB2_ML_c,  "\n",
"Scaling  factor      = ", scale_factor_SB2_c,  "\n" ,  "\n" )


cat("T_SB2_RLS_c : ", "Sattora-Bentler adjusted test based on T_RLS with improved scaling factor", "(option: estimator = MLMVS)", "\n",
"Test Statistic       = ", T_SB2_RLS_c,  "\n", 
"Degrees of freedom   = ", df_SB2_RLS_c, "\n", 
"P-value (Chi-square) = ", p_SB2_RLS_c,  "\n",
"Scaling  factor      = ", scale_factor_SB2_c,  "\n" )

###############################################################################
# Fit indices
T_base  <- fitMeasures(fit1, "baseline.chisq")
df_base <- p*(p+1)/2 - p

NFI_ML <- fitMeasures(fit1, "nfi")
TLI_ML <- fitMeasures(fit1, "tli")
CFI_ML <- fitMeasures(fit1, "cfi")
IFI_ML <- fitMeasures(fit1, "ifi")
MFI_ML <- fitMeasures(fit1, "mfi")
RMSEA_ML <- fitMeasures(fit1, "rmsea")

NFI_RLS <- (T_base - T_RLS)/T_base
TLI_RLS <- ((T_base/df_base) - (T_RLS/df_RLS))/(T_base/df_base-1)
CFI_RLS <- 1 - max(c(T_RLS-df_RLS, 0))/max(c(T_RLS-df_RLS, T_base-df_base, 0))
IFI_RLS <- (T_base - T_RLS)/(T_base - df_RLS)
MFI_RLS <- exp(-0.5*((T_RLS - df_RLS)/N))
RMSEA_RLS <- sqrt( max(c((T_RLS - df_RLS)/(N-1),0))/df_RLS )

cat("**************************************", "\n")
cat("Fit indices: ", "\n")
cat(
"\t    ", "ML",    "\t  ",  " RLS ",    "\n",
"NFI   = ", NFI_ML,   "\t",  NFI_RLS,   "\n",
"TLI   = ", TLI_ML,   "\t",  TLI_RLS,   "\n",
"CFI   = ", CFI_ML,   "\t",  CFI_RLS,   "\n",
"IFI   = ", IFI_ML,   "\t",  IFI_RLS,   "\n",
"MFI   = ", MFI_ML,   "\t",  MFI_RLS,   "\n",
"RMSEA = ", RMSEA_ML, "\t",  RMSEA_RLS, "\n" )




 }