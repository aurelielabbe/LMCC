#' @title Mater5_2 kernel
#' @description
#' Compute the Matern5_2 kernel.
#' @param d the distance.
#' @param beta the inverse of range parameter.
#' @references Gu, M., & Xu, Y. (2020).
#' Fast nonseparable Gaussian stochastic process with application to methylation level interpolation.
#' Journal of Computational and Graphical Statistics, 29(2), 250-260.
#' \doi{10.1080/10618600.2019.1665534}
#' @examples
#'
#'matern_5_2_kernel(0.5, beta = 2)
#'
#' @export
#'

matern_5_2_kernel<-function(d,beta){
  res=(1+sqrt(5)*beta*d + 5*beta^2*d^2/3 )*exp(-sqrt(5)*beta*d)
  res
}
