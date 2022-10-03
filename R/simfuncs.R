#' @export 
generateThreeNormal = function(n, p = 2, seed = 1, with_noise = TRUE, noise_p = 1, noise_sd = 1)
{
  set.seed(seed)
  nclass = 3
  each_n = floor(n / nclass)
  n = each_n * nclass
  x1 = matrix(rnorm(each_n * p, c(0, -2), 1),
            nrow = each_n, ncol = p, byrow = TRUE)
  x2 = matrix(rnorm(each_n * p, c(2, 2), 1),
              nrow = each_n, ncol = p, byrow = TRUE)
  x3 = matrix(rnorm(each_n * p, c(-2, 2), 1),
              nrow = each_n, ncol = p, byrow = TRUE)

  X = rbind(x1, x2, x3)
  y = rep(c(1, 2, 3), each = each_n)
  
  if (with_noise) {
    noise_dat = matrix(rnorm(n * noise_p, sd = noise_sd), n, noise_p)
    X = cbind(X, noise_dat)
  }
  return(list(x = X, y = y))
}



#' @export 
generateSmiley = function(n, p = 2, seed = 1, with_noise = TRUE, noise_p = 1, noise_sd = 2)
{
  set.seed(seed)
  nclass = 3
  each_n = floor(n / nclass)
  n = each_n * nclass
  sigma = 0.1
  x1 = matrix(rnorm(each_n * p, c(-1, -0.5), 0.2),
            nrow = each_n, ncol = p, byrow = TRUE)
  x2 = matrix(rnorm(each_n * p, c(1, -0.5), 0.2),
              nrow = each_n, ncol = p, byrow = TRUE)
  seq_x = runif(each_n, pi, 2 * pi)
  x3 = cbind(3 * cos(seq_x) + rnorm(each_n) * sigma, 2 * sin(seq_x) + 0.0 + rnorm(each_n) * sigma)

  X = rbind(x1, x2, x3)
  y = rep(c(1, 2, 3), each = each_n)
  if (with_noise) {
    noise_dat = matrix(rnorm(n * noise_p, sd = noise_sd), n, noise_p)
    X = cbind(X, noise_dat)
  }
  return(list(x = X, y = y))
}


#' @export 
generateTaegeuk = function(n, p = 2, seed = 1, with_noise = TRUE, noise_p = 1, noise_sd = 2)
{
  set.seed(seed)
  nclass = numeric(2)
  each_n = floor(n / 2)
  n = each_n * 2
  X = matrix(nrow = n, ncol = p)
  y = numeric(n)
  k = 1

  while (k <= n) {
    x = rnorm(p, sd = 2)
    sx = sum(x^2)
    if (sx < 4 & x[2] > (2 / 3 * sin(x[1] * 1 / 2 * pi)) + 1 / 3) {
      nclass[1] = nclass[1] + 1
      y[k] = 1
      X[k, ] = x
      k = k + 1
    } else if (sx < 4 & x[2] < (2 / 3 * sin(x[1] * 1 / 2 * pi)) - 1 / 3) {
      nclass[2] = nclass[2] + 1
      y[k] = 2
      X[k, ] = x
      k = k + 1
    }
  }

  if (with_noise) {
    noise_dat = matrix(rnorm(n * noise_p, sd = noise_sd), n, noise_p)
    X = cbind(X, noise_dat)
  }
  return(list(x = X, y = y))
}


#' @export 
generateMultiorange = function(n, p = 2, seed = 1, with_noise = TRUE, noise_p = 1, noise_sd = 2)
{
  set.seed(seed)
  nclass = numeric(3)
  each_n = floor(n / 3)
  n = each_n * 3
  X = matrix(nrow = n, ncol = p)
  y = numeric(n)
  k = 1
  while (k <= n) {
    x = rnorm(p, sd = 2)
    sx = sum(x^2)
    if (sx <= 0.5 & nclass[1] < each_n) {
      nclass[1] = nclass[1] + 1
      y[k] = 1
      X[k, ] = x
      k = k + 1
    }
    else if (2.0 < sx & sx <= 3 & nclass[2] < each_n) {
      nclass[2] = nclass[2] + 1
      y[k] = 2
      X[k, ] = x
      k = k + 1
    }
    else if (6 < sx & sx <= 7 & nclass[3] < each_n) {
      nclass[3] = nclass[3] + 1
      y[k] = 3
      X[k, ] = x
      k = k + 1
    }
  }
  if (with_noise) {
    noise_dat = matrix(rnorm(n * noise_p, sd = noise_sd), n, noise_p)
    X = cbind(X, noise_dat)
  }
  return(list(x = X, y = y))
}

#' @export 
generateTwoorange = function(n, p = 2, seed = 1, with_noise = TRUE, noise_p = 1, noise_sd = 2)
{
  set.seed(seed)
  nclass = numeric(2)
  each_n = floor(n / 2)
  n = each_n * 2
  X = matrix(nrow = n, ncol = p)
  y = numeric(n)
  k = 1
  while (k <= n) {
    x = rnorm(p, sd = 2)
    sx = sum(x^2)
    if (sx <= 0.5 & nclass[1] < each_n) {
      nclass[1] = nclass[1] + 1
      y[k] = 1
      X[k, ] = x
      k = k + 1
    }else if (5 < sx & sx <= 7 & nclass[2] < each_n) {
      nclass[2] = nclass[2] + 1
      y[k] = 2
      X[k, ] = x
      k = k + 1
    }
  }
  if (with_noise) {
    noise_dat = matrix(rnorm(n * noise_p, sd = noise_sd), n, noise_p)
    X = cbind(X, noise_dat)
  }
  return(list(x = X, y = y))
}

#' @export 
generateMultiMoon = function(each_n = 100, sigma = 1, noise_p = 4, noise_sd = 3, seed = NULL)
{
  set.seed(seed)
  x = runif(each_n, 0, pi)
  # c1 = cbind(5 * cos(x) - 3.5 + rnorm(each_n) * sigma, 10 * sin(x) -
  #              2.5 + rnorm(each_n) * sigma)
  # x = runif(each_n, pi, 2 * pi)
  # c2 = cbind(5 * cos(x) + 3.5 + rnorm(each_n) * sigma, 10 * sin(x) +
  #              0.5 + rnorm(each_n) * sigma)
  # x = runif(each_n, 0, pi)
  # c3 = cbind(5 * cos(x) + 10.5 + rnorm(each_n) * sigma, 10 * sin(x) -
  #              2.5 + rnorm(each_n) * sigma)
  c1 = cbind(7.5 * cos(x) - 5.5 + rnorm(each_n) * sigma, 10 * sin(x) -
               2.5 + rnorm(each_n) * sigma)
  x = runif(each_n, pi, 2 * pi)
  c2 = cbind(7.5 * cos(x) + 3.5 + rnorm(each_n) * sigma, 10 * sin(x) +
               3.5 + rnorm(each_n) * sigma)
  x = runif(each_n, 0, pi)
  c3 = cbind(7.5 * cos(x) + 12.5 + rnorm(each_n) * sigma, 10 * sin(x) -
               2.5 + rnorm(each_n) * sigma)
  X = rbind(c1, c2, c3)
  noise_X = matrix(rnorm(3 * each_n * noise_p, 0, noise_sd), nrow = 3 * each_n, ncol = noise_p)
  X = cbind(X, noise_X)
  y = rep(c(1, 2, 3), each = each_n)
  return(list(x = X, y = y))
}

#' @export 
generateTwoMoon = function(each_n = 100, sigma = 1, noise_p = 4, noise_sd = 3, seed = NULL)
{
  set.seed(seed)
  x = runif(each_n, 0, pi)
  c1 = cbind(7.5 * cos(x) - 3.5 + rnorm(each_n) * sigma, 10 * sin(x) -
               2.5 + rnorm(each_n) * sigma)
  x = runif(each_n, pi, 2 * pi)
  c2 = cbind(7.5 * cos(x) + 3.5 + rnorm(each_n) * sigma, 10 * sin(x) +
               2.5 + rnorm(each_n) * sigma)
  
  X = rbind(c1, c2)
  noise_X = matrix(rnorm(2 * each_n * noise_p, 0, noise_sd), nrow = 2 * each_n, ncol = noise_p)
  X = cbind(X, noise_X)
  y = rep(c(1, 2), each = each_n)
  return(list(x = X, y = y))
}