tune.kkmeans = function(x, nCluster, nPerms = 20, nStart = 10, weights = NULL,
                        kernel = "linear", kparam = 1, opt = TRUE, nCores = 1, ...)
{
  out = list()
  call = match.call()
  kernel = match.arg(kernel, c("linear", "gaussian"))
  
  x = as.matrix(x)
  n = nrow(x)

  perm_list = vector("list", nPerms)
  for (i in 1:nPerms) {
    perm_list[[i]] = sapply(1:p, function(j) sample(x[, j]))
  }

  org_bcd = unlist(parallel::mclapply(1:length(kparam), FUN = function(j) {
    org_fit = kkmeans(x, nCluster = nCluster, nStart = nStart, weights = weights,
                      kernel = kernel, kparam = kparam[j], opt = TRUE, ...)
    return(org_fit$maxBcd)
  }, mc.cores = nCores))


  perm_bcd_list = matrix(0, nrow = nPerms, ncol = length(kparam))
  for (b in 1:nPerms) {
    perm_bcd = unlist(parallel::mclapply(1:length(kparam), FUN = function(j) {
      perm_fit = kkmeans(x = perm_list[[b]], nCluster = nCluster, nStart = nStart, weights = weights,
                         kernel = kernel, kparam = kparam[j], opt = TRUE, ...)
      return(perm_fit$maxBcd)
    }, mc.cores = nCores))
    
    perm_bcd_list[b, ] = perm_bcd
  }
  out$orgBcd = org_bcd
  out$permBcd = perm_bcd_list
  out$gaps = log(org_bcd) - colMeans(log(perm_bcd_list))
  out$optInd = min(which(out$gaps == max(out$gaps)))
  out$opt_kparam = kparam[out$optInd]
  
  if (opt) {
    opt_fit = kkmeans(x = x, nCluster = nCluster, nStart = nStart, weights = weights,
                  kernel = kernel, kparam = out$opt_kparam, opt = TRUE, ...)  
    out$optModel = opt_fit
  }
  out$call = call
  return(out)
}


kkmeans = function(x, nCluster, nStart = 10, weights = NULL,
                   kernel = "linear", kparam = 1, opt = TRUE, ...) 
{
  out = list()
  call = match.call()
  kernel = match.arg(kernel, c("linear", "gaussian"))
  kernel = switch(kernel, 
                  "linear" = "vanilladot",
                  "gaussian" = "rbfdot")
  
  if (kernel == "linear") {
    kernel_fun = kernlab::vanilladot()
    kpar = list()
  } else if (kernel == "rbfdot") {
    kernel_fun = kernlab::rbfdot(sigma = kparam)
    kpar = list(sigma = kparam)
  }

  x = as.matrix(x)
  n = nrow(x)
  
  if (is.null(weights)) {
    weights = rep(1, n)
  }
  
  Kmat = list()
  Kmat$K = list(kernlab::kernelMatrix(kernel_fun, x, x))
  Kmat$numK = 1
  td = GetWCD(Kmat, clusters = rep(1, n), weights = weights)  

  res = vector("list", length = nStart)
  wcd = numeric(length(nStart))
  seeds = seq(1, nStart, by = 1)
  for (j in 1:length(seeds)) {
    
    try_error = try({
      res[[j]] = kernlab::kkmeans(x = x, centers = nCluster, 
                                 kernel = kernel, kpar = kpar, ...)
    })
    if (inherits(try_error, "try-error")) {
      wcd[j] = Inf
    } else {
      wcd[j] = GetWCD(Kmat, clusters = res[[j]]@.Data, weights = weights)  
    }
  }
  bcd = td - wcd
  optInd = which.max(bcd)
  if (opt) {
    out$optRes = res[[optInd]]
    out$optClusters = out$optRes@.Data
    out$maxBcd = bcd[optInd]
  }
  out$td = td
  out$wcd = wcd
  out$bcd = bcd
  out$res = res
  return(out)
}
