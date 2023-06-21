#' @title Approximate Mater5_2 kernel
#' @description
#' Compute the approximate Matern5_2 kernel.
#' @param a parameter.
#' @param b parameter.
#' @param C parameter.
#' @param beta_i the inverse of range parameter.
#' @param eta_i the nugget parameter.
#' @references Gu, M., & Xu, Y. (2020). Fast nonseparable Gaussian stochastic process with application to methylation level interpolation. Journal of Computational and Graphical Statistics, 29(2), 250-260.
#' @examples
#'
#'approx_ref_matern_5_2(1,1/2,0.1,0.5, 2)
#'
#' @export
#'

approx_ref_matern_5_2<-function(a,b,C,beta_i,eta_i ){
  t=C*beta_i+eta_i
  a*log(t)-b*t
}
