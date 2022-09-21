rm(list = ls())
gc()

# setwd(r"(C:\Users\Beom\GitHub\SparseKK-means)")
setwd(r"(C:\Users\user\GitHub\SparseKK-means)")

require(sparcl)
require(kernlab)
require(caret)
require(fossil)
source(r"(.\R\subfuncs.R)")
source(r"(.\R\main.R)")
source(r"(.\R\simfuncs.R)")

n = 100
p = 2
noise_p = c(0, 5, 25, 50)
# noise_p = c(0)
iter = 100

# save ARI results
# ARI_mat = matrix(NA, nrow = iter, ncol = 3)
# colnames(ARI_mat) = c("skkm", "skm", "kkm")
ari_list = list()

# save fitting results
skkm_list = skm_list = kkm_list = list()
skkm_res_list = skm_res_list = kkm_res_list = list()

# save time results
# time_mat = matrix(NA, nrow = iter, ncol = 3)
# colnames(time_mat) = c("skkm", "skm", "kkm")
time_list = list() 

i = 1
j = 4

for (j in j:length(noise_p)) {
  i = 1
  
  ARI_mat = matrix(NA, nrow = iter, ncol = 3)
  colnames(ARI_mat) = c("skkm", "skm", "kkm")
  
  time_mat = matrix(NA, nrow = iter, ncol = 3)
  colnames(time_mat) = c("skkm", "skm", "kkm")
  
  for (i in i:iter) {
    cat(j, "th setting", i, "th iteration \n")
    # dat = generateMultiorange(n = n, p = p, seed = 2, with_noise = TRUE, noise_p = 5)
    dat = generateTwoorange(n = n, p = p, seed = i, with_noise = TRUE, noise_p = noise_p[j])
    # dat = generateMultiMoon(each_n = n, sigma = 0.5, seed = 1, noise_p = 5, noise_sd = 3)
    # dat = generateTwoMoon(each_n = n, sigma = 0.5, seed = 1, noise_p = 5, noise_sd = 3)
    
    # sigma = kernlab::sigest(scale(dat$x), scale = FALSE)[3]
    # sigma = 1.5
    sigma = 1.0
    
    # Sparse kernel k-means algorithm
    skkm_t = system.time({
      tuned_skkm = tune.skkm(x = dat$x, nCluster = 2, s = NULL, ns = 10, nPerms = 25,
                             nStart = 1, kernel = "gaussian-2way", kparam = sigma, opt = TRUE,
                             nInit = 20)
    })
    skkm_clusters = tuned_skkm$optModel$opt_clusters
    ari_skkm = adj.rand.index(dat$y, skkm_clusters)
    ARI_mat[i, "skkm"] = ari_skkm
    skkm_list[[i]] = tuned_skkm
    time_mat[i, "skkm"] = skkm_t[3]
    # plot(dat$x[, 1:2], col = skkm_clusters,
    #      pch = 16, cex = 1.5,
    #      xlab = "x1", ylab = "y1", main = "Proposed method")
    
    
    # Sparse k-means algorithm
    skm_t = system.time({
      tuned_skm = KMeansSparseCluster.permute(x = dat$x, K = 2, nvals = 10, nperms = 25, silent = TRUE)
      opt_skm = KMeansSparseCluster(x = dat$x, K = 2, wbounds = tuned_skm$bestw, silent = TRUE)
    })
    skm_clusters = opt_skm[[1]]$Cs
    ari_skm = adj.rand.index(dat$y, skm_clusters)
    ARI_mat[i, "skm"] = ari_skm
    skm_list[[i]] = opt_skm
    time_mat[i, "skm"] = skm_t[3]
    
    # Kernel k-means algorithm
    kkm_t = system.time({
      kkm_fit_list = list()
      kkm_fit_wcd = c()
      for (kk in 1:20) {
        KKK = list()
        kkm_res_temp = kkmeans(dat$x, centers = 2, kernel = "rbfdot", kpar = list(sigma = sigma))
        kkm_fit_list[[kk]] = kkm_res_temp
        
        # computing within-cluster distance of kkmeans
        KKK$K = list(kernlab::kernelMatrix(rbfdot(sigma = sigma), dat$x, dat$x))
        KKK$numK = 1
        kkm_fit_wcd[kk] = GetWCD(KKK, clusters = kkm_res_temp@.Data, weights = rep(1, nrow(dat$x)))
      }
      kkm_res = kkm_fit_list[[which.min(kkm_fit_wcd)]]
      # plot(dat$x[, 1:2], col = kkm_res@.Data,
      #      pch = 16, cex = 1.5,
      #      xlab = "x1", ylab = "y1", main = "kkmeans")
    })
    kkm_clusters = kkm_res@.Data
    ari_kkm = adj.rand.index(dat$y, kkm_clusters)
    ARI_mat[i, "kkm"] = ari_kkm
    kkm_list[[i]] = kkm_res
    time_mat[i, "kkm"] = kkm_t[3]
    
    save.image("./orange_simulation_n=100_20220919_scale=FALSE_sigma=1_p50.Rdata")
  }
  skkm_res_list[[j]] = skkm_list
  skm_res_list[[j]] = skm_list
  kkm_res_list[[j]] = kkm_list
  
  ari_list[[j]] = ARI_mat
  time_list[[j]] = time_mat
  
  
  save.image("./orange_simulation_n=100_20220919_scale=FALSE_sigma=1_p50.Rdata")
}

sapply(ari_list, colMeans, na.rm = TRUE)


