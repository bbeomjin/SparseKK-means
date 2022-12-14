#' tune.skkm
#' 
#' tuning skkm parameter
#' @importFrom parallel mclapply
#' @param x input data
#' @return s tuning parameter
#' @export 
tune.skkm = function(x, nCluster, nPerms = 20, s = NULL, ns = 100, nStart = 10, weights = NULL, 
                     kernel = "linear", kparam = 1, search = "exact", opt = TRUE, nCores = 1, ...)
{
  out = list()
  call = match.call()
  kernel = match.arg(kernel, c("linear", "gaussian", "spline-t",
                               "gaussian-2way", "spline-t-2way"))
  search = match.arg(search, c("exact", "binary"))
  if (!is.matrix(x)) {
    x = as.matrix(x)
  }
  
  p = ncol(x)
  
  if (is.null(s)) {
    if (kernel %in% c("gaussian-2way", "spline-t-2way")) {
      nv = p + p * (p - 1) / 2
    } else {
      nv = p
    }
    s = exp(seq(log(1), log(sqrt(nv)), length.out = ns))
  } else {
    s = sort(s)
  }
  
  params = expand.grid(s = s, kparam = kparam)

  perm_list = vector("list", nPerms)
  for (i in 1:nPerms) {
    perm_list[[i]] = sapply(1:p, function(j) sample(x[, j]))
  }
    
  org_bcd = unlist(parallel::mclapply(1:nrow(params), FUN = function(j) {
    org_fit = skkm(x = x, nCluster = nCluster, nStart = nStart, s = params$s[j], weights = weights,
                   kernel = kernel, kparam = params$kparam[j], opt = TRUE, ...)
    return(org_fit$maxBcd)
  }, mc.cores = nCores))
  
  # org_bcd = unlist(parallel::mclapply(s, FUN = function(ss) {
  #   org_fit = skkm(x = x, nCluster = nCluster, nStart = nStart, s = ss, weights = weights,
  #                  kernel = kernel, kparam = kparam, opt = TRUE, ...)
  #   return(org_fit$maxBcd)
  # }, mc.cores = nCores))


  perm_bcd_list = matrix(0, nrow = nPerms, ncol = nrow(params))
  for (b in 1:nPerms) {
    perm_bcd = unlist(parallel::mclapply(1:nrow(params), FUN = function(j) {
      perm_fit = skkm(x = perm_list[[b]], nCluster = nCluster, nStart = nStart, s = params$s[j], weights = weights,
                      kernel = kernel, kparam = params$kparam[j], opt = TRUE, ...)
      return(perm_fit$maxBcd)
    }, mc.cores = nCores))
   
    perm_bcd_list[b, ] = perm_bcd
  }

  # perm_bcd_list = matrix(0, nrow = nPerms, ncol = ns)
  # for (b in 1:nPerms) {
  #   perm_bcd = unlist(parallel::mclapply(s, FUN = function(ss) {
  #     perm_fit = skkm(x = perm_list[[b]], nCluster = nCluster, nStart = nStart, s = ss, weights = weights,
  #                     kernel = kernel, kparam = kparam, opt = TRUE, ...)
  #     return(perm_fit$maxBcd)
  #   }, mc.cores = nCores))
    
  #   perm_bcd_list[b, ] = perm_bcd
  # } 

  out$orgBcd = org_bcd
  out$permBcd = perm_bcd_list
  out$gaps = log(org_bcd) - colMeans(log(perm_bcd_list))
  out$optInd = min(which(out$gaps == max(out$gaps)))
  # out$opt_s = s[out$optInd]
  out$opt_s = params[out$optInd, "s"]
  out$opt_kparam = params[out$optInd, "kparam"]
    
  if (opt) {
    opt_fit = skkm(x = x, nCluster = nCluster, nStart = nStart, s = out$opt_s, weights = weights,
                  kernel = kernel, kparam = out$opt_kparam, opt = TRUE, ...)  
    out$optModel = opt_fit
  }
  
  # if (opt) {
  #   opt_fit = skkm(x = x, nCluster = nCluster, nStart = nStart, s = out$opt_s, weights = weights,
  #                 kernel = kernel, kparam = kparam, opt = TRUE, ...)  
  #   out$optModel = opt_fit
  # }

  out$call = call
  return(out)
}

#' skkm
#' 
#' fitting skkm 
#' @param x input data
#' @return out list of clusters
#' @export 
skkm = function(x, nCluster, nStart = 10, s = 1.5, weights = NULL,
               kernel = "linear", kparam = 1, search = "exact", opt = TRUE, ...) 
{
  out = list()
  call = match.call()
  kernel = match.arg(kernel, c("linear", "gaussian", "spline-t",
                               "gaussian-2way", "spline-t-2way"))
  search = match.arg(search, c("exact", "binary"))
  x = as.matrix(x)
  n = nrow(x)
  # p = ncol(x)
  
  if (is.null(weights)) {
    weights = rep(1, n)
    # attr(weights, "type") = "auto"
  }
  
  res = vector("list", length = nStart)
  seeds = seq(1, nStart, by = 1)
  for (j in 1:length(seeds)) {
    # initialization
    # set.seed(seeds[j])
    # clusters0 = sample(1:nCluster, size = n, replace = TRUE)
    # aa = make_anovaKernel(x, x, kernel = kernel, kparam = sigma)
    # theta = rep(1 / sqrt(3), 3)
    # K = combine_kernel(aa, theta)
    # fit = kkmeans2(K, centers = nCluster)
    # clusters0 = fit@.Data
    
    res[[j]] = skkm_core(x = x, clusters = nCluster, theta = NULL, s = s, weights = weights,
                         kernel = kernel, kparam = kparam, search = search, ...)
  }
  if (opt) {
    bcd_list = sapply(res, function(x) {
      bcd = max(x$bcd)
    })
    
    optInd = which(bcd_list == max(bcd_list))
    if (length(optInd) > 1) {
      warning("")
      optInd = optInd[1]
    }
    
    out$optClusters = res[[optInd]]$clusters
    out$optTheta = res[[optInd]]$theta
    out$maxBcd = bcd_list[optInd]
  }
  out$res = res
  return(out)
}

skkm_core = function(x, clusters = NULL, nInit = 20, theta = NULL, s = 1.5, weights = NULL,
               kernel = "linear", kparam = 1, search = "exact", maxiter = 100, eps = 1e-8) 
{
  out = list()
  call = match.call()
  n = nrow(x)
  # p = ncol(x)
  
  # initialization
  anovaKernel = make_anovaKernel(x = x, y = x, kernel = kernel, kparam = kparam)
  theta0 = theta
  
  if (is.null(theta0)) {
    theta0 = rep(1 / sqrt(anovaKernel$numK), anovaKernel$numK)
  }
  
  if (is.null(weights)) {
    weights = rep(1, n)
  }
  
  if (length(clusters) == 1) {
    nCluster = clusters
    init_wcd_vec = numeric(nInit)
    init_clusters_list = vector("list", nInit)
    for (i in 1:nInit) {
      clusters0 = sample(1:nCluster, size = n, replace = TRUE)
      clusters = updateCs(anovaKernel = anovaKernel, theta = theta0, 
                          clusters = clusters0, weights = weights)$clusters
      init_clusters_list[[i]] = clusters
      wcd = GetWCD(anovaKernel, clusters = clusters, weights = weights)
      init_wcd_vec[i] = sum(theta0 * wcd)
    }
    init_clusters = clusters0 = init_clusters_list[[which.min(init_wcd_vec)]]
  } else {
    init_clusters = clusters0 = clusters
  }
  
  td_vec = wcd_vec = bcd_vec = c()
  
  for (iter in 1:maxiter) {
    
    # Update clusters
    clusters = updateCs(anovaKernel = anovaKernel, theta = theta0, 
                        clusters = clusters0, weights = weights)$clusters  
    # plot(dat$x[, 1:2], col = clusters)
    
    # Update theta
    wcd = GetWCD(anovaKernel, clusters = clusters, weights = weights)
    td = GetWCD(anovaKernel, rep(1, length(clusters)), weights = weights)
    bcd = td - wcd
    
    if (search == "exact") {
      terror = try({
        suppressWarnings({delta = ExactSearch(coefs = bcd, s = s)})
      }, silent = TRUE)
      if (inherits(terror, "try-error")) {
        delta = BinarySearch(coefs = bcd, s = s)  
      }  
    } else {
      delta = BinarySearch(coefs = bcd, s = s)
    }
    
    theta_tmp = soft_threshold(bcd, delta = delta)
    theta = normalization(theta_tmp)
    
    # td_new = GetWCD(anovaKernel, rep(1, length(clusters)), weights = weights)
    # wcd_new = GetWCD(anovaKernel, clusters = clusters, weights = weights)
    # sum(theta0 * td) - sum(theta * td)
    # theta0 * td
    # theta * td
    # sum(theta0 * wcd) - sum(theta * wcd)
    
    td_vec[iter] = sum(theta * td)
    wcd_vec[iter] = sum(theta * wcd)
    bcd_vec[iter] = sum(theta * bcd)
    
    # print((sum(abs(theta - theta0)) / sum(theta0)))
    if ((sum(abs(theta - theta0)) / (sum(theta0) + 1e-12)) < eps) {
      break
    } else {
      theta0 = theta
      clusters0 = clusters
    }
  }
  out$clusters = clusters
  out$theta = theta
  out$weights = weights
  out$td = td_vec
  out$wcd = wcd_vec
  out$bcd = bcd_vec
  out$init_clusters = init_clusters
  out$iteration = iter
  return(out)
}
