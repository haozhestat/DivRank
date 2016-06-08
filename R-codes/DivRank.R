# =============================================================================
# File Name:     DivRank.R
# Author:        Haozhe Zhang
# Contact:       haozhe@iastate.edu
# Creation Date: 2016-06-07
# Last Modified: 2016-06-08 13:28:18 CDT
# =============================================================================

DivRank = function(weight_mat, lambda, alpha, prior_score=NULL,converg_error=10^(-3),max_iter=1000){
  n = nrow(weight_mat)
  if(is.null(prior_score))
    prior_score = rep(1/n,n)
  num_iter = 0
  diff_error = 10^6
  marg_prob = rep(1/n,n)
  prob_int_mat = as.matrix(alpha*weight_mat/matrix(rep(apply(weight_mat, 1, sum),n),n,n))
  diag(prob_int_mat) = (1-alpha)
  while((num_iter<max_iter)&(diff_error>converg_error)){
    tmp = prob_int_mat*t(matrix(rep(marg_prob,n),n,n))
    prob_mat = (1-lambda)*t(matrix(rep(prior_score,n),n,n))+lambda*tmp/matrix(rep(rowSums(tmp),n),n,n)
    marg_prob_new = marg_prob%*%prob_mat
    num_iter = num_iter + 1
    diff_error = sum(abs(marg_prob_new-marg_prob))/sum(abs(marg_prob))
    marg_prob = marg_prob_new
  }
  return(list(divrank=order(marg_prob,decreasing=TRUE),marg_prob=as.vector(marg_prob),num_iter=num_iter))
}
