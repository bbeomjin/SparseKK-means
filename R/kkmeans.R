# #' tune.skkm
# #' 
# #' tuning kkmeans parameter
# #' @importFrom parallel mclapply
# #' @importFrom kernlab kkmeans
# #' @param x input data
# #' @return out results of kkmeans
# #' @export 
# tune.kkmeans = function(x, nCluster, nPerms = 20, nStart = 10, weights = NULL,
#                         kernel = "linear", kparam = 1, opt = TRUE, nCores = 1, ...)
# {
#   out = list()
#   call = match.call()
#   kernel = match.arg(kernel, c("linear", "gaussian"))
#   
#   x = as.matrix(x)
#   n = nrow(x)
# 
#   perm_list = vector("list", nPerms)
#   for (i in 1:nPerms) {
#     perm_list[[i]] = sapply(1:p, function(j) sample(x[, j]))
#   }
# 
#   org_bcd = unlist(parallel::mclapply(1:length(kparam), FUN = function(j) {
#     org_fit = kkmeans(x = x, nCluster = nCluster, nStart = nStart, weights = weights,
#                       kernel = kernel, kparam = kparam[j], opt = TRUE, ...)
#     return(org_fit$maxBcd)
#   }, mc.cores = nCores))
# 
# 
#   perm_bcd_list = matrix(0, nrow = nPerms, ncol = length(kparam))
#   for (b in 1:nPerms) {
#     perm_bcd = unlist(parallel::mclapply(1:length(kparam), FUN = function(j) {
#       perm_fit = kkmeans(x = perm_list[[b]], nCluster = nCluster, nStart = nStart, weights = weights,
#                          kernel = kernel, kparam = kparam[j], opt = TRUE, ...)
#       return(perm_fit$maxBcd)
#     }, mc.cores = nCores))
#     
#     perm_bcd_list[b, ] = perm_bcd
#   }
#   out$orgBcd = org_bcd
#   out$permBcd = perm_bcd_list
#   out$gaps = log(org_bcd) - colMeans(log(perm_bcd_list))
#   out$optInd = min(which(out$gaps == max(out$gaps)))
#   out$opt_kparam = kparam[out$optInd]
#   
#   if (opt) {
#     opt_fit = kkmeans(x = x, nCluster = nCluster, nStart = nStart, weights = weights,
#                       kernel = kernel, kparam = out$opt_kparam, opt = TRUE, ...)  
#     out$optModel = opt_fit
#   }
#   out$call = call
#   return(out)
# }

#' kkmeans
#' 
#' fitting kkmeans 
#' @importFrom kernlab kkmeans
#' @param x input data
#' @return out list of clusters
#' @export 
kkmeans = function(x, nCluster, nStart = 10, weights = NULL,
                   kernel = "linear", kparam = 1, opt = TRUE, ...) 
{
  out = list()
  call = match.call()
  kernel = match.arg(kernel, c("linear", "gaussian"))
  
  # if (kernel == "linear") {
  #   kernel_fun = kernlab::vanilladot()
  #   kpar = list()
  # } else if (kernel == "gaussian") {
  #   kernel_fun = kernlab::rbfdot(sigma = kparam)
  #   kpar = list(sigma = kparam)
  # }

  x = as.matrix(x)
  n = nrow(x)
  
  if (is.null(weights)) {
    weights = rep(1, n)
  }
  
  Kmat = list()
  # Kmat$K = list(kernlab::kernelMatrix(kernel_fun, x, x))
  Kmat$K = list(kernelMatrix(x, x, kernel = kernel, kparam = kparam))
  Kmat$numK = 1
  td = GetWCD(Kmat, clusters = rep(1, n), weights = weights)  

  res = vector("list", length = nStart)
  wcd = numeric(length(nStart))
  seeds = seq(1, nStart, by = 1)
  for (j in 1:length(seeds)) {
    
    try_error = try({
      # res[[j]] = kernlab::kkmeans(x = x, centers = nCluster, 
      #                            kernel = kernel, kpar = kpar, ...)
      # res[[j]] = kernlab::kkmeans(Kmat$K[[1]], centers = nCluster, ...)
      clusters0 = sample(1:nCluster, size = n, replace = TRUE)
      res[[j]] = updateCs(anovaKernel = Kmat, theta = 1, 
                          clusters = clusters0, weights = weights)
    })
    if (inherits(try_error, "try-error")) {
      wcd[j] = Inf
    } else {
      # wcd[j] = GetWCD(Kmat, clusters = res[[j]]@.Data, weights = weights)  
      wcd[j] = GetWCD(Kmat, clusters = res[[j]]$clusters, weights = weights)  
    }
  }
  bcd = td - wcd
  optInd = which.max(bcd)
  if (opt) {
    out$optRes = res[[optInd]]
    out$optClusters = out$optRes$clusters
    out$maxBcd = bcd[optInd]
  }
  out$td = td
  out$wcd = wcd
  out$bcd = bcd
  out$optInd = optInd
  out$res = res
  return(out)
}