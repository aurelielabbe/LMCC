#' @title Singular Value Decomposition of X and Y
#' @description
#' Perform a SVD of the covariates matrix X and the incomplete matrix Y
#' @param Y_obs an incomplete matrix, data frame or data table to be imputed. Missing values are coded as NAs.
#' @param sites a numeric vector of the site locations with a length equal to the number of columns of Y_obs.
#' @param X a matrix, data frame or data table of covariates with the same number of rows as Y_obs.
#' @param tol_eig the threshold for the eigen values of the SVD decomposition.
#' @return
#' a list of objects needed for the LMCC model estimation and prediction:
#' \item{Y_obs}{incomplete data matrix}
#' \item{input}{observed methylation sites (scaled)}
#' \item{output}{observed data matrix}
#' \item{scale_sites}{methylation sites (scaled)}
#' \item{A}{matrices A and X stacked by columns (see reference below)}
#' \item{A_output_t}{matrices B and V stacked by rows (see reference below)}
#' \item{X}{covariates matrix}
#'
#' @author Melina Ribaud
#' @references Melina Ribaud, Aurélie Labbe and Karim Oualkacha.
#' Imputation in genetic methylation studies: LMCC, a Linear Model of Coregionalization with informative Covariates.
#' 2023. hal-00000000
#'
#' @examples
#'
#'library(FastGP)
#'N = 100
#'K = 10
#'D = 2
#'sites = sort(runif(N))
#'beta = c(runif(D,2,5),runif(K-D,10,1000))
#'nugget = c(rep(0,D),runif(K-D,0.001,0.05))
#'A = matrix(runif(K*(K-D),0,0.1),nrow = K, ncol = K-D )
#'X = matrix(runif(D*K),ncol = D, nrow= K)
#'Hx = matrix(solve(t(X)%*%X,t(X)),nrow = length(X)/K)
#'A = A-X%*%Hx%*%A
#'A = cbind(X,A)
#'V = matrix(NA, nrow = K, ncol = N)
#'R00 = abs(outer(sites, sites, '-'))
#'for (d in 1:K) {
#'  R = matern_5_2_kernel(R00, beta = beta[d])
#'  R_tilde = R + nugget[d] * diag(N)
#'  V[d, ] = rcpp_rmvnorm_stable(1, R, rep(0, N))
#'}
#'Y_obs = A %*% V
#'obj_lmcc = svd_lmcc(Y_obs,sites,X,tol_eig = 1e-6)
#' @export
#'



svd_lmcc = function(Y_obs,sites,X,tol_eig = 0){

  if(sum(sites != sort(sites))>0){
    stop("sites must be sorted")
  }

  if(length(ncol(X))<1){
    nrow_X = length(X)
  }else{
    nrow_X = nrow(X)
  }


  K = nrow(Y_obs)
  D = K - (length(X)/K)
  n_index = apply(is.na(Y_obs),2,sum)==0
  scale_sites = scale_01(sites)
  input = scale_sites[n_index]

  if (length(input) < K) {
    warning(
      "ncol(Y_obs) without NA is lower than K, the SVD of I-Px(Y) will not be of rank K"
    )
  }

  if (K != nrow_X) {
    stop(
      "nrow(X) must be equal to nrow(Y_obs)"
    )
  }

  output = Y_obs[,n_index]
  C=(max(input)-min(input))/ncol(output)
  rowMeans_t_output=rep(0,K)

  output_mean = output - rowMeans_t_output



  svd_X = svd(X)
  eig_X = svd_X$d
  select_eig = eig_X >= tol_eig
  eig_X[eig_X<=0] = 1e-30
  D_svd = diag(eig_X[select_eig],sum(select_eig),sum(select_eig))
  X_svd = svd_X$u[,1:sum(select_eig)]%*%D_svd

  inv_D_svd_2 = diag(1/(eig_X[select_eig]^2),sum(select_eig),sum(select_eig))
  B_svd = inv_D_svd_2%*%t(X_svd)%*%output_mean




  Y_k = X_svd%*%B_svd
  X = X_svd
  B = B_svd




  svd_A = svd(output_mean-Y_k)
  eig_A = svd_A$d
  select_eig = eig_A > tol_eig
  eig_A[eig_A<=0] = 1e-30
  D_svd = diag(eig_A[select_eig],sum(select_eig),sum(select_eig))
  A = -(svd_A$u[,1:sum(select_eig)]%*%D_svd/sqrt(ncol(output)))[,1:sum(select_eig)]
  A_output_t = t(-svd_A$v[,1:sum(select_eig)]*sqrt(ncol(output)))


  A = cbind(X,A)
  A_output_t = rbind(B,A_output_t)

  obj = list(
    Y_obs = Y_obs,
    input = input,
    output = output,
    C = C,
    rowMeans_t_output = rowMeans_t_output,
    scale_sites = scale_sites,
    A = A,
    A_output_t = A_output_t,
    X = X
  )

  return(obj)
}
