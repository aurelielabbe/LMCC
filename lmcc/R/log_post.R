#' @title Log posterior of Matern5_2 kernel
#' @description
#' Compute the log posterior of the Matern5_2 kernel.
#' @param param hyperparameters.
#' @param object a FastGaSP object.
#' @param C parameter.
#' @references Gu, M., & Xu, Y. (2020).
#' Fast nonseparable Gaussian stochastic process with application to methylation level interpolation.
#' Journal of Computational and Graphical Statistics, 29(2), 250-260.
#' \doi{10.1080/10618600.2019.1665534}
#' @examples
#'
#'sites = sort(runif(100))
#'Y = sin(2*pi*sites)
#'C = diff(range(sites))/length(sites)
#'fgasp.model = fgasp(sites,Y, have_noise = TRUE)
#'log_post(c(log(1 / C), 1),fgasp.model,C)
#'
#' @export
#'

log_post<-function(param,object,C){
  log_lik_val=log_lik(param,object)
  beta=exp(param[1])
  eta=exp(param[2])  ####constraint
  a=1/2
  b=1
  logprior=approx_ref_matern_5_2(a,b,C,beta,  eta)
  log_lik_val+ logprior
}
