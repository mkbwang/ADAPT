library(nloptr)


#### first part: try with different outcomes
eval_f <- function(x)
{
  return ( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 )
}


eval_grad_f <- function(x) {
  return( c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
             200 * (x[2] - x[1] * x[1]) ) )
}

x0 <- c( -1.2, 1 )

opts <- list("algorithm"="NLOPT_LD_LBFGS",
             "xtol_rel"=1.0e-8)

# solve Rosenbrock Banana function
res <- nloptr( x0=x0,
               eval_f=eval_f,
               eval_grad_f=eval_grad_f,
               opts=opts)


print(res)


## Rosenbrock Banana function and gradient in one function
eval_f_list <- function(x) {
  common_term <- x[2] - x[1] * x[1]
  return( list( "objective" = 100 * common_term^2 + (1 - x[1])^2,
                "gradient"  = c( -400 * x[1] * common_term - 2 * (1 - x[1]),
                                 200 * common_term) ) )
}


res <- nloptr( x0=x0,
               eval_f=eval_f_list,
               opts=opts)
print( res )


##### minimization with inequality constraints

# objective function
eval_f0 <- function( x, a, b ){
  return( sqrt(x[2]) )
}

# gradient of objective function
eval_grad_f0 <- function( x, a, b ){
  return( c( 0, .5/sqrt(x[2]) ) )
}


# constraint function
eval_g0 <- function( x, a, b ) {
  return( (a*x[1] + b)^3 - x[2] )
}

# jacobian of constraint
eval_jac_g0 <- function( x, a, b ) {
  return( rbind( c( 3*a[1]*(a[1]*x[1] + b[1])^2, -1.0 ),
                 c( 3*a[2]*(a[2]*x[1] + b[2])^2, -1.0 ) ) )
}


# define parameters
a <- c(2,-1)
b <- c(0, 1)

# Solve using NLOPT_LD_MMA with gradient information supplied in separate function
res0 <- nloptr( x0=c(1.234,5.678),
                eval_f=eval_f0,
                eval_grad_f=eval_grad_f0,
                lb = c(-Inf,0),
                ub = c(Inf,Inf),
                eval_g_ineq = eval_g0,
                eval_jac_g_ineq = eval_jac_g0,
                opts = list("algorithm" = "NLOPT_LD_MMA",
                            "xtol_rel"=1.0e-8,
                            "print_level" = 2,
                            "check_derivatives" = TRUE,
                            "check_derivatives_print" = "all"),
                a = a,
                b = b )


# Solve using NLOPT_LN_COBYLA without gradient information
res1 <- nloptr( x0=c(1.234,5.678),
                eval_f=eval_f0,
                lb = c(-Inf,0),
                ub = c(Inf,Inf),
                eval_g_ineq = eval_g0,# smaller than zero
                opts = list("algorithm"="NLOPT_LN_COBYLA",
                            "xtol_rel"=1.0e-8),
                a = a,
                b = b ) # minimize function
print( res1 )


# prepare X matrix and V diagonal matrix, maximize B

proflik <- function(B, V, uu, K, X){
  R <- B+V
  invmatR <- diag(1/R)

  invXRX <- function(xvecs, invR){
    XRX <- t(xvecs) %*% invR %*% xvecs
    eigen_result <- eigen(XRX)
    return (eigen_result$vectors %*% diag(1/eigen_result$values) %*% t(eigen_result$vectors))
  }
  URX <- t(uu) %*% invmatR %*% X

  URU <- t(uu) %*% invmatR %*% uu

  result <- sum(log(R)) + K * log(URU - URX %*% invXRX(X, invmatR) %*% t(URX))
  return(as.vector(result))
}

UU0 <- c(7,3,6,5)
V0 <- c(2,5,3,6)
X0 <- matrix(c(1,1,1,1,6,7,12,4), nrow=4)
K0 <- 4

res1 <- nloptr( x0=10,
                eval_f=proflik,
                lb = 0,
                ub = 300,
                opts = list("algorithm"="NLOPT_LN_COBYLA",
                            "xtol_rel"=1.0e-8),
                V = V0,
                uu = UU0,
                K = K0,
                X = X0) # minimize function



