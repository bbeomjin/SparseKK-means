kernelMatrix = function(x, y, kernel = "gaussian", kparam = 1.0) {
  
  x = as.matrix(x)
  y = as.matrix(y)
  p = ncol(x)
  
  if (NCOL(x) == 0) {
    x = matrix(0, nrow = nrow(x), ncol = 1)
  }
  
  if (NCOL(y) == 0) {
    y = matrix(0, nrow = nrow(y), ncol = 1)
  }
  
  if (kernel == "poly") {
    K = (x %*% t(y) + 1.0)^kparam
  } else if(kernel == "gaussian" | kernel == "gaussian-2way") {
    normx = rowSums(x^2)
    normy = rowSums(y^2)
    temp = x %*% t(y)
    temp = (-2.0 * temp) + outer(normx, rep(1.0, nrow(y)), "*") + outer(rep(1.0, nrow(x)), normy, "*")
    K = exp(-temp * kparam)
    # K = kernlab:::kernelMatrix(rbfdot(sigma = kparam), x, y)
  } else if (kernel == "linear") {
    K = tcrossprod(x, y)
  } else if (kernel == "anova_gaussian") {
    K = 0
    for (d in 1:p) {
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = kernelMatrix(A, B, kernel = "gaussian", kparam = kparam)
      K = K + K_temp
    }
  } else {
    K = NULL
  }
  return(K)
}

make_anovaKernel = function(x, y, kernel, kparam)
{
  x = as.matrix(x)
  y = as.matrix(y)
  dimx = ncol(x)
  
  # calculate anova kernels for main effects
  if (kernel == "spline") {
    # assign the number of anova kernels
    numK = 2 * dimx
    # list of kernel matrices
    anova_kernel = vector(mode = "list", numK)
    # list of kernel coordinate indices
    kernelCoord = vector(mode = "list", numK)
    index = 0
    
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep="")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep="")
    }
    
  } else if (kernel == 'spline-2way') {
    numK = (2 * dimx) + (2 * dimx * (2 * dimx - 1) / 2 - dimx)
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    # main effects
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep = "")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep = "")
    }
    # two-way interactions
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A_linear = as.matrix(anova_kernel[[2 * i - 1]])
        A_smooth = as.matrix(anova_kernel[[2 * i]])
        B_linear = as.matrix(anova_kernel[[2 * j - 1]])
        B_smooth = as.matrix(anova_kernel[[2 * j]])
        anova_kernel[[index]] = A_linear * B_linear
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " linear", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_linear * B_smooth
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " smooth", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_smooth * B_linear
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " linear", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_smooth * B_smooth
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " smooth", sep = "")
      }
    }
  } else if (kernel == "spline-t") {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
  } else if (kernel == 'spline-t-2way') {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
    
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep = "")
      }
    }
  } else if (kernel == "gaussian-2way") {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      anova_kernel[[index]] = kernelMatrix(A, B, kernel, kparam)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
    
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep = "")
      }
    }
  } else {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    for (d in 1:dimx) {
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      anova_kernel[[d]] = kernelMatrix(A, B, kernel, kparam)
      kernelCoord[[d]] = paste("x", d, sep = "")
    }
  }
  return(list(K = anova_kernel, coord = kernelCoord, numK = numK, kernel = kernel, kparam = kparam))
}

spline_kernel = function(x, u)
{
  x = as.matrix(x)
  u = as.matrix(u)
  K1x = (x - 1 / 2)
  K1u = (u - 1 / 2)
  K2x = (K1x^2 - 1 / 12) / 2
  K2u = (K1u^2 - 1 / 12) / 2
  ax = x %x% matrix(1, 1, nrow(u))
  au = u %x% matrix(1, 1, nrow(x))
  b = abs(ax - t(au))
  K1 = K1x %x% t(K1u)
  K2 = K2x %x% t(K2u) - ((b - 1 / 2)^4 - (b - 1 / 2)^2 / 2 + 7 / 240) / 24
  return(list(K1 = K1, K2 = K2))
}

combine_kernel = function(anovaKernel, 
                          theta = rep(1 / sqrt(anovaKernel$numK),
                                      anovaKernel$numK))
{
  K = 0
  for (v in 1:length(theta)) {
    K = (K + theta[v] * anovaKernel$K[[v]])
  }
  return(K)
}

# make_a_vec = function(anovaKernel, clusters) {
#   uc = sort(unique(clusters))
#   d = anovaKernel$numK
#   n = length(clusters)
#   K1 = 2 * n * sapply(anovaKernel$K, function(x) sum(diag(x)))
#   K2 = 2 * sapply(anovaKernel$K, sum)
#   
#   K3 = K4 = numeric(d)
#   for (j in 1:d) {
#     r_part = lapply(uc, function(i) {
#       gind = clusters == i
#       temp_K = anovaKernel$K[[j]][gind, gind]
#       temp_K3 = sum(gind) * sum(diag(temp_K))
#       temp_K4 = sum(temp_K)
#       return(list(temp_K3, temp_K4))
#     })
#     K3[j] = sum(sapply(r_part, "[[", 1))
#     K4[j] = sum(sapply(r_part, "[[", 2))
#   }
#   return(K1 - K2 - K3 + K4)
# }

GetWCD = function(anovaKernel, clusters, weights)
{
  uc = unique(clusters)
  d = anovaKernel$numK
  wcd = numeric(d)
  for (v in 1:d) {
    K = anovaKernel$K[[v]]
    
    for (g in 1:length(uc)) {
      ind = clusters == uc[g]
      swt = weights[ind]
      subK = K[ind, ind]
      # wcd[v] = wcd[v] + sum(diag(subK)) - (1 / sum(swt)) * sum((subK * tcrossprod(swt)))
      wcd[v] = wcd[v] + sum(swt * diag(subK)) - (1 / sum(swt)) * sum((subK * tcrossprod(swt)))
    }
  }
  return(wcd)
  # return(list(td = td, wcd = wcd))
}


updateCs = function(anovaKernel, theta, clusters, weights, maxiter = 100) {
  
  # Initialization
  clusters0 = clusters
  
  for (i in 1:maxiter) {
    uc = unique(clusters0)
    RKHS_dist = sapply(uc, function(g) {
      
      ind = clusters0 == g
      swt = weights[ind]
      
      K = combine_kernel(anovaKernel, theta = theta)
      
      Kxx = diag(K)
      
      # Kxy = list()
      # Kxy$K = lapply(anovaKernel$K, function(x) x[ind, , drop = FALSE] * swt)
      # Kxy = colSums(combine_kernel(Kxy, theta))
      Kxy = colSums(K[ind, , drop = FALSE] * swt)
      # Kxy = (K[, ind, drop = FALSE] %*% swt) / sum(swt)
      
      # Kyy = list()
      # Kyy$K = lapply(anovaKernel$K, function(x) x[ind, ind, drop = FALSE] * tcrossprod(swt))
      # Kyy = sum(combine_kernel(Kyy, theta))
      Kyy = sum(K[ind, ind, drop = FALSE] * tcrossprod(swt))
      # Kyy = sum(drop(crossprod(K[ind, ind], swt)) * swt) / sum(swt)^2
      
      return(Kxx - (2 * Kxy / sum(swt)) + (Kyy / sum(swt)^2))
      # return(Kxx - 2 * Kxy + Kyy)
      # return(list(Kxx, Kxy, Kyy))
    })
    clusters = uc[apply(RKHS_dist, 1, which.min)]
    if (sum(clusters != clusters0) == 0) {
      break
    } else {
      clusters0 = clusters
    }
  }
  return(list(clusters = clusters, iteration = i, RKHS_dist = RKHS_dist))
}


# RKHS_dist2 = function(x, theta, g, kernel, kparam) {
#   new_x = x
#   G = unique(g)
#   split_X = lapply(G, function(i) x[g == i, , drop = FALSE])
#   g_index_temp = sapply(split_X, function(x) {
#     n = nrow(x)
#     anova_tmp = make_anovaKernel(new_x, new_x, kernel = kernel, kparam = kparam)
#     anova_K = combine_kernel(anova_tmp, theta = theta)
#     
#     anova_temp2 = make_anovaKernel(new_x, x, kernel = kernel, kparam = kparam)
#     anova_temp3 = make_anovaKernel(x, x, kernel = kernel, kparam = kparam)
#     
#     K1 = diag(anova_K)
#     K2 = 2 * rowSums(combine_kernel(anova_temp2, theta)) / n
#     # K2 = rowSums(combine_kernel(anova_temp2, theta))
#     
#     K3 = sum(combine_kernel(anova_temp3, theta)) / n^2
#     # K3 = sum(combine_kernel(anova_temp3, theta))
#     return(K1 + K3 - K2)
#     # return(list(K1, K2, K3))
#   })
#   return(g_index_temp)
# }


soft_threshold = function(x, delta) {
  w = sign(x) * pmax(abs(x) - delta, 0)
  return(w)
}

normalization = function(x) {
  return(x / sqrt(sum(x^2)))
}

# the function from sparcl package
BinarySearch = function(coefs, s) 
{
  if((sum(coefs^2) == 0) | (sum(abs(normalization(coefs))) <= s)) return(0)
  lamb1 = 0
  lamb2 = max(abs(coefs)) - 1e-6
  iter = 0
  while ((iter <= 30) & ((lamb2 - lamb1) > 1e-5)) {
    iter = iter + 1
    w_tmp = soft_threshold(coefs, (lamb1 + lamb2) / 2)
    w = normalization(w_tmp)
    if (sum(abs(w)) < s) {
      lamb2 = (lamb1 + lamb2) / 2
    } else {
      lamb1 = (lamb1 + lamb2) / 2
    }
  }
  return((lamb1 + lamb2) / 2)
}

#
ExactSearch = function(coefs, s)
{
  p = length(coefs)
  sorted_coefs = sort(coefs, decreasing = TRUE)
  lambda_vec = numeric(p)
  for (i in 1:p) {
    sub_coefs = coefs[1:i]
    ft = sum(sub_coefs) / i
    st_tmp = sum(sub_coefs)^2 / i - (sum(sub_coefs)^2 - sum(sub_coefs^2) * s^2) / (i - s^2)
    st_tmp[abs(st_tmp) < 1e-12] = 0
    if (is.nan(st_tmp) | (abs(st_tmp) == Inf) | (st_tmp < 0)) {
      lambda_vec[i] = Inf
    } else {
      st = (1 / sqrt(i)) * sqrt(st_tmp)
      lambda_vec[i] = min(ft + st, ft - st)
    }
  }
  return(min(lambda_vec, na.rm = TRUE))
}


