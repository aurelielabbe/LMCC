#' @title LMCC model: fit
#' @description
#' Estimate the parameters of the LMCC model
#' @param obj a list of objects created by the svd_lmcc function.
#' @return
#' a list of objects needed for estimation and prediction, among others
#'
#' \item{Y_obs}{incomplete data matrix}
#' \item{input}{observed methylation sites (scaled)}
#' \item{output}{observed data matrix}
#' \item{scale_sites}{methylation sites (scaled)}
#' \item{A}{matrices A and X stacked by columns (see reference below)}
#' \item{A_output_t}{matrices B and V stacked by rows (see reference below)}
#' \item{beta_record}{estimated range parameters}
#' \item{eta_record}{estimated nugget parameters}
#' \item{sigma2_record}{estimated variance parameters}
#' \item{X}{covariates matrix}
#' @author Melina Ribaud
#' @references Melina Ribaud, Aurélie Labbe and Karim Oualkacha.
#' Imputation in genetic methylation studies: LMCC, a Linear Model of Coregionalization with informative Covariates.
#' 2023. hal-00000000
#'
#' Gu, M., & Xu, Y. (2020).
#' Fast nonseparable Gaussian stochastic process with application to methylation level interpolation.
#' Journal of Computational and Graphical Statistics, 29(2), 250-260.
#' \doi{10.1080/10618600.2019.1665534}
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
#'obj_lmcc = fit_lmcc(obj_lmcc)
#' @export
#'


fit_lmcc = function(obj){
  Y_obs = input = output = C = rowMeans_t_output = scale_sites = A = A_output_t = X = NULL

  for (ind in 1:length(obj))
  {
    assign(names(obj)[ind], obj[[ind]])
  }

  D = ncol(A)

  beta_record=rep(0,D)
  eta_record=rep(0,D)
  val_record=rep(0,D)
  sigma2_record=rep(0,D)



  for(d in 1:D){
    #print(d)
    fgasp.model=fgasp(input,A_output_t[d,], have_noise = T)
    test = try({
      tt_all <-
        optim(c(log(1 / C), 1), function(par, object)
          return(log_post(par, object, C)), object = fgasp.model, lower = c(-100,-15),
          method = "L-BFGS-B",
          control = list(fnscale = -1, maxit = 30))
    })

    count = 0
    while (class(test)[1] == "try-error" & count<=100) {
      count = count + 5
      test = try({
        tt_all <-
          optim(c(log(1 / C), 1), function(par, object)
            return(log_post(par, object, C)), object = fgasp.model, lower = c(-100/count,-15),
            method = "L-BFGS-B",
            control = list(fnscale = -1, maxit = 30))
      })
    }

    #tt_all
    val_record[d]=tt_all$value

    beta_record[d]=exp(tt_all$par)[1]
    eta_record[d]=exp(tt_all$par)[2]

    sigma2_record[d]=Get_log_det_S2(param=log(c(beta_record[d],eta_record[d])),fgasp.model@have_noise,fgasp.model@delta_x,
                                    fgasp.model@output,fgasp.model@kernel_type)[[2]]/fgasp.model@num_obs

  }


  obj = list(
    Y_obs = Y_obs,
    input = input,
    output = output,
    rowMeans_t_output = rowMeans_t_output,
    scale_sites = scale_sites,
    A = A,
    A_output_t = A_output_t,
    beta_record = beta_record,
    eta_record = eta_record,
    val_record = val_record,
    sigma2_record = sigma2_record,
    X = X)

  return(obj)
}
