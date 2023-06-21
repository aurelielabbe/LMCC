#' @title LMCC model: prediction of missing values
#' @description
#' Predict missing data values using the LMCC model
#' @param obj a list of objects created by the function fit_lmcc.
#' @return a data matrix completed with missing value predictions
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
#'Ypred = pred_lmcc(obj_lmcc)
#' @export
#'


pred_lmcc = function(obj){

  Y_obs = input = output = C = rowMeans_t_output = scale_sites = A = A_output_t = X = beta_record = eta_record = val_record = sigma2_record = NULL

  for (ind in 1:length(obj))
  {
    assign(names(obj)[ind], obj[[ind]])
  }


  K = nrow(Y_obs)
  D = ncol(A)

  ncol_X = 0
  if(length(X)>0){
    ncol_X = 1:(length(X)/K)
  }

  predict_all_dlm = sd_all_dlm = matrix(0,D, length(scale_sites))

  for(d in 1:D){


    fgasp.model=fgasp(input,A_output_t[d,],have_noise = T)
    pred_fast=predict(param=log(c(beta_record[d],eta_record[d])),object=fgasp.model,
                      testing_input=scale_sites)

    predict_all_dlm[d,]=pred_fast@mean
    pred_var = pred_fast@var
    pred_var[pred_var<0] = 0
    sd_all_dlm[d,]=sqrt(pred_var)
  }

  pred_all_record=A%*%(predict_all_dlm)

  pred_testing_record_cond = var_testing_record_cond = matrix(NA,K,length(scale_sites))


  for(i in 1:length(scale_sites)){
    if(i%%10000==0){
      print(i)
    }
    testing_ppl_i = is.na(Y_obs[,i])

    if(sum(testing_ppl_i) == nrow(Y_obs)){
      pred_testing_record_cond[testing_ppl_i,i]=pred_all_record[testing_ppl_i,i]
      var_predict=Sigma_hat[testing_ppl_i,testing_ppl_i]
      var_testing_record_cond[testing_ppl_i,i]=diag(var_predict)

    }else if(sum(testing_ppl_i)>0){

      Sigma_hat=A%*%diag(sd_all_dlm[,i]^2)%*%t(A)
      inv_sig = Sigma_hat[!testing_ppl_i,!testing_ppl_i]

      test = try({inv_sig = solve(Sigma_hat[!testing_ppl_i,!testing_ppl_i])},silent = T)
      if (class(test)[1] == "try-error") {
        inv_sig = solve(Sigma_hat[!testing_ppl_i,!testing_ppl_i] +
                          diag(rep(1e-6, sum(!testing_ppl_i))))
      }
      useful_block=Sigma_hat[testing_ppl_i,!testing_ppl_i]%*%inv_sig

      pred_testing_record_cond[testing_ppl_i,i]=pred_all_record[testing_ppl_i,i]+
        useful_block%*%(Y_obs[!testing_ppl_i,i]-rowMeans_t_output[!testing_ppl_i] -
                          pred_all_record[!testing_ppl_i,i])+ rowMeans_t_output[testing_ppl_i]


      var_predict=Sigma_hat[testing_ppl_i,testing_ppl_i]-useful_block%*%Sigma_hat[!testing_ppl_i,testing_ppl_i]
      var_testing_record_cond[testing_ppl_i,i]=diag(var_predict)

    }

    if(sum(!testing_ppl_i)>0){
      pred_testing_record_cond[!testing_ppl_i,i]=Y_obs[!testing_ppl_i,i]
      var_testing_record_cond[!testing_ppl_i,i] = 0
    }


  }
  return(pred_testing_record_cond)
}
