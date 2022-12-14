rm(list = ls())
gc()

# setwd(r"(C:\Users\Beom\GitHub\SparseKK-means)")
setwd(r"(C:\Users\user\GitHub\SparseKK-means)")

# remove.packages("skkm")
# devtools:::install_github("bbeomjin/skkm")
source("./R/main.R")
source("./R/kkmeans.R")
source("./R/subfuncs.R")
source("./R/simfuncs.R")

require(sparcl)
require(kernlab)
require(caret)
require(fossil)
# require(skkm)

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
j = 1

##### just test ########
dat = generateTaegeuk(n = n, seed = i, noise_p = noise_p[j])
nclusters = 2
# dat = generateSmiley(n = n, p = p, seed = i, with_noise = TRUE, noise_p = noise_p[j], noise_sd = 2)
x = dat$x
nCluster = nclusters
s = NULL
ns = 20
nPerms = 25
nStart = 1
kernel = "gaussian-2way"
kparam = c(0.25, 0.5, 0.75, 1)
opt = TRUE
nInit = 20
weights = NULL
nCores = 1
p = ncol(dat$x)
nv = p + p * (p - 1) / 2
s = sqrt(nv)
search = "exact"

sigmas = seq(0.1, 1, by = 0.1)
# sigma = 1
# sigma = 0.5
# sigma = 0.25
# sigma = 0.1

ari = c()
tttt = c()
cccc = c()
for (sig in sigmas) {
  cat(sig)
   skkm_t = system.time({
      tuned_skkm = tune.skkm(x = x, nCluster = 2, s = NULL, ns = 10, nPerms = 25,
                             nStart = 1, kernel = "gaussian-2way", kparam = sig, opt = TRUE,
                             nInit = 20)
    })
   aaa = adj.rand.index(dat$y, tuned_skkm$optModel$optClusters)    
   ari = append(ari, aaa)
   anovaK = make_anovaKernel(x, x, kernel = "gaussian-2way", kparam = sig)
   K = combine_kernel(anovaK, tuned_skkm$optModel$optTheta)
   e = eigen(K, symmetric = TRUE)
   e$values[e$values < 0] = 0
   kk = diag(sqrt(e$values)) %*% e$vectors
   ttt = 0
   for (i in unique(tuned_skkm$optModel$optClusters)) {
     Ksub = kk[tuned_skkm$optModel$optClusters == i, , drop = FALSE]
     ttt = ttt + sum(rowSums((Ksub - rowMeans(Ksub))^2))
   }
  tttt = append(tttt, ttt) 
  ccc = GetWCD(anovaK, clusters = tuned_skkm$optModel$optClusters, weights = rep(1, n))
  cccc = append(cccc, sum(tuned_skkm$optModel$optTheta * ccc))
}


ari
tttt
which.min(tttt)

adj.rand.index(dat$y, tuned_skkm$optModel$optClusters)

aa = tune.kkmeans(x, 2, nPerms = 25, nStart = 20, kernel = "gaussian", kparam = sigmas)
aa$gaps
sigmas[which.max(aa$gaps)]

# tuned_skkm$optModel$optTheta

##########################

tt = tune.kkmeans(x = x, nCluster = 2, nPerms = 25, nSpart = 20, kernel = "gaussian", kparam = c(0.5, 1))
tt$optModel$optClusters


for (j in j:length(noise_p)) {
  i = 1
  
  ARI_mat = matrix(NA, nrow = iter, ncol = 3)
  colnames(ARI_mat) = c("skkm", "skm", "kkm")
  
  time_mat = matrix(NA, nrow = iter, ncol = 3)
  colnames(time_mat) = c("skkm", "skm", "kkm")
  
  for (i in i:iter) {
    cat(j, "th setting", i, "th iteration \n")
    # dat = generateMultiorange(n = n, p = p, seed = 2, with_noise = TRUE, noise_p = 5)
    # dat = generateTwoorange(n = n, p = p, seed = i, with_noise = TRUE, noise_p = noise_p[j])
    dat = generateSmiley(n = n, p = p, seed = i, with_noise = TRUE, noise_p = noise_p[j], noise_sd = 2)
    # dat = generateMultiMoon(each_n = n, sigma = 0.5, seed = 1, noise_p = 5, noise_sd = 3)
    # dat = generateTwoMoon(each_n = n, sigma = 0.5, seed = 1, noise_p = 5, noise_sd = 3)
    
    # sigma = kernlab::sigest(scale(dat$x), scale = FALSE)[3]
    # sigma = 1.5
    sigma = 1.0
    
    # Sparse kernel k-means algorithm
    skkm_t = system.time({
      tuned_skkm = tune.skkm(x = dat$x, nCluster = 3, s = NULL, ns = 10, nPerms = 25,
                             nStart = 1, kernel = "gaussian-2way", kparam = sigma, opt = TRUE,
                             nInit = 20)
    })
    skkm_clusters = tuned_skkm$optModel$opt_clusters
    ari_skkm = adj.rand.index(dat$y, skkm_clusters)
    ARI_mat[i, "skkm"] = ari_skkm
    skkm_list[[i]] = tuned_skkm
    time_mat[i, "skkm"] = skkm_t[3]
    # tuned_skkm$opt_s
    # tuned_skkm$optModel$opt_theta
    plot(dat$x[, 1:2], col = skkm_clusters,
          pch = 16, cex = 1.5,
          xlab = "x1", ylab = "y1", main = "Proposed method")
    
    
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


