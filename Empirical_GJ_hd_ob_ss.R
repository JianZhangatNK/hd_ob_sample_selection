# Empirical Estimation: OB decomposition in a High-dimensional Sample Selection Framework
#  Features: Outcome variable Y is binary;
# 
#
# variables:
# y1 - hourly wage, y2 - monthly wage, y3 - yearly wage, y4 - health, y5 - good job indicator
# D1 - hukou, D2 - male, D3 - formal, D4-minority
# x - other control variables
#
#
# Clear all and load packages:
rm(list = ls())
library('MASS')
library('Matrix')
library('VGAM')
library('foreach')
library('glmnet')
library('bigmemory')
library('pracma')
library('parallel')
library('doSNOW')
library('snow')
library('bigmemory')

# source file:
# Our Functions:
source('source_hdobss.R')

# Auto-ML
source('primitives.R')
source('stage1v3.R')

###### Loading data ####
load("OBhimdimcfps2018v2GJ.RData")


##### plug in values ######
n.total.sample = length(y1)
ste.sample.index = sample(1:n.total.sample, n.total.sample, replace = FALSE)   #19480 c(1:n.total.sample) 46397

### change y1~y4 and D1~D4 to consider different dependent variables and treatments
y = y5[ste.sample.index]
y.label = "good job"
D = as.logical(D1[ste.sample.index]) # Reverse baseline group: use 1-D1 or 1-D2 instead of D1 or D2
D.label = "Hukou:Urban" # change label accordingly

w = w.cross[ste.sample.index,] # z.cross = (x.c#x.d, x.d#x.d)
w.d = w.d.cross[ste.sample.index,] # z.d.cross = x.d#x.d
x = x.zl[ste.sample.index,] # x.zl = (x.c, x.d)
x.c.ep = x.c[ste.sample.index,]
n = length(y) # sample size

# z.data = cbind(z1, z2, z1^2, z2^2) 
z.data = cbind(z1, z1^2)
z = z.data[ste.sample.index,]

#######################################
############ Sample selection indicator
S = (y1[ste.sample.index]>0)
######################################
######################################

##### grid points #####
ngrid = 1
ylim.list = 1 #quantile(y, c(1:ngrid)/(1+ngrid))

# y.gap = (y.max-y.min)/(ngrid+1)
# ngrid = 4
# ylim.list = c(1,2,3,4)

##### Create Folds####
##### random splitting, f-fold: #####
##### first (f-1) folds: size = floor(n/f) 
##### The f-th fold: size = n - floor(n/f)*(f-1)
ndata.rd = sample(1:n, n, replace = FALSE)
f = 10
n_int = floor(n/f)*f
datafolds = matrix(ndata.rd[1:n_int], floor(n/f), f)
newcol = ndata.rd[(floor(n/f)*(f-1)+1):n]
datafolds.list1 = as.list(as.data.frame(datafolds[,1:(f-1)]))
datafolds.list2 = list(newcol)
names(datafolds.list2) = paste0("v",f)
datafolds.list = c(datafolds.list1, datafolds.list2)

### Weighted bootstrap
B = 500
b.xi = matrix(rexp(B*n, rate = 1), B, n)

##########################
##########################
### Estimation of rho ###
# Go through some bugs
s=S
xpoly = cbind(x, x.c.ep^2, w)
xzpoly = cbind(xpoly, z) # If z is binary, delete z^2 term
p=dim(xpoly)[2]

# Set up grid points
tgrid = seq(0.01, 0.99, 0.01) # t=F*(y|x)=G(x*beta(y))
n_tgrid = length(tgrid)
n_rho = 19 # rho^d  copula coefficients
rho_grid = seq(from = -0.95, to = 0.95, length.out = n_rho)


# Find rho_d first
# Step 1
s0 = s[D==0]
xz0poly = xzpoly[D==0,]
s0.cv.out = cv.glmnet(xz0poly, s0, type.measure = 'mse', family = "binomial")
s0coeff = coefficients(s0.cv.out, s = s0.cv.out$lambda.min)
s0_index0 = s0coeff@i[-1] + 1

s1 = s[D==1]
xz1poly = xzpoly[D==1,]
s1.cv.out = cv.glmnet(xz1poly, s1, type.measure = 'mse', family = "binomial")
s1coeff = coefficients(s1.cv.out, s = s1.cv.out$lambda.min)
s1_index0 = s1coeff@i[-1] + 1

# Step 2
xz01poly = xzpoly[(D==0)&(s==1),]
xz11poly = xzpoly[(D==1)&(s==1),]

FS_g01_coef = vector("list", ngrid)
FS_g11_coef = vector("list", ngrid)
FS_g01_index = vector("list", ngrid)
FS_g11_index = vector("list", ngrid)

for (jj in c(1:ngrid)) {
  cat('jj =', jj, '\n')
  y_grid = ylim.list[jj]
  # d=0
  yvar_g = (y==y_grid) 
  yvar_g01 = yvar_g[(D==0)&(s==1)]
  y_g01.cv.out = cv.glmnet(xz01poly, yvar_g01, type.measure = 'mse')
  y_g01coeff = coefficients(y_g01.cv.out, s = y_g01.cv.out$lambda.min)
  y_g01_index0 = y_g01coeff@i[-1] + 1
  
  FS_g01_coef[[jj]] = y_g01coeff
  FS_g01_index[[jj]] = y_g01_index0
  
  # d=1
  yvar_g11 = yvar_g[(D==1)&(s==1)]
  y_g11.cv.out = cv.glmnet(xz11poly, yvar_g11, type.measure = 'mse')
  y_g11coeff = coefficients(y_g11.cv.out, s = y_g11.cv.out$lambda.min)
  y_g11_index0 = y_g11coeff@i[-1] + 1
  
  FS_g11_coef[[jj]] = y_g11coeff
  FS_g11_index[[jj]] = y_g11_index0
  
}

PI0hat_ii = as.big.matrix(matrix(NA, n, n))
PI1hat_ii = as.big.matrix(matrix(NA, n, n))
for (ii in c(1:n)) {
  cat('ii =', ii, '\n')
  xi_vec = matrix(rep(x[ii,], each = n), n, dim(x)[2])
  x.c.ep_ii = matrix(rep(x.c.ep[ii,], each = n), n, dim(x.c.ep)[2])
  w_ii = matrix(rep(w[ii,], each = n), n, dim(w)[2])
  xi_zpoly = cbind(xi_vec, x.c.ep_ii^2, w_ii, z)
  
  # d = 0
  PI0_ii.index= cbind(0, xi_zpoly[,s0_index0-1]) %*%
    c(0, s0coeff[s0_index0]) + s0coeff[1]
  PI0hat_ii[ii,] = as.vector(exp(PI0_ii.index)/(1+exp(PI0_ii.index))) # Pi_d(x,Z)
  
  # d = 1
  PI1_ii.index= cbind(0, xi_zpoly[,s1_index0-1]) %*%
    c(0, s1coeff[s1_index0]) + s1coeff[1]
  PI1hat_ii[ii,] = as.vector(exp(PI1_ii.index)/(1+exp(PI1_ii.index))) # Pi_d(x,Z)
}

save.image(file = "OBhimdimcfps2018Mainv2GJ_temp1.RData")
PI0hat_ii = as.matrix(PI0hat_ii)
PI1hat_ii = as.matrix(PI1hat_ii)

# Step 3 ## change to parallel 
### some variables to save time
d_z = dim(z)[2]
d_x = dim(xpoly)[2]


### Parallel Setup
cl <- makeCluster(10) ##### change according to the number of your PC cores
clusterEvalQ(cl, {
  library('VGAM')
  library('MASS')
  library('Matrix') 
  library('pracma')
  library('glmnet')
  library('bigmemory')
  library('primes')
})
varlist = c("n.total.sample", "ste.sample.index",  
            "w",  "z", "x", "x.c.ep",  "ylim.list", "xpoly", "d_z", "d_x",
            "n", "n_rho", "rho_grid", "ngrid", "n_tgrid", "tgrid", 'S', 's', 'PI0hat_ii', 'PI1hat_ii',
            'FS_g01_coef', 'FS_g11_coef', 'FS_g01_index', 'FS_g11_index',
            'Copula.n')
clusterExport(cl, varlist, envir = .GlobalEnv)

par_find_rho <- function(kk){
  #cat('loop_rho2_kk =', kk, '\n')
  obj_0_rho = 0
  obj_1_rho = 0
  that_0_xyr = matrix(NA, n, ngrid)
  that_1_xyr = matrix(NA, n, ngrid)
  
  rho_d = rho_grid[kk]
  for (ii in c(1:n)) {
   # cat('loop_rho2_ii =', ii, '\n')
    xi_poly = c(xpoly[ii,], rep(0, d_z))
    
    for (jj in c(1:ngrid)) {
      #cat('loop_rho2_jj =', jj, '\n')
      # d = 0
      y_g01coeff = FS_g01_coef[[jj]]
      y_g01_index0 = FS_g01_index[[jj]]
      FS_g01.index1= sum(c(0, xi_poly[y_g01_index0-1]) *
        c(0, y_g01coeff[y_g01_index0])) + y_g01coeff[1]
      y_g01_index0_z = y_g01_index0[y_g01_index0>(d_x+1)]
      FS_g01.index2= c(cbind(0, z[,y_g01_index0_z-1-d_x]) %*%
        c(0, y_g01coeff[y_g01_index0_z]))
      FS_g01hat = FS_g01.index1 + FS_g01.index2 # F^s(y|x,z)
      
      
      # d = 1
      y_g11coeff = FS_g11_coef[[jj]]
      y_g11_index0 = FS_g11_index[[jj]]
      FS_g11.index1= sum(c(0, xi_poly[y_g11_index0-1]) *
        c(0, y_g11coeff[y_g11_index0])) + y_g11coeff[1]
      y_g11_index0_z = y_g11_index0[y_g11_index0>(d_x+1)]
      FS_g11.index2= c(cbind(0, z[,y_g11_index0_z-1-d_x]) %*%
        c(0, y_g11coeff[y_g11_index0_z]))
      FS_g11hat = FS_g11.index1 + FS_g11.index2 # F^s(y|x,z)
      
      obj_0_t = rep(NA, n_tgrid)
      obj_1_t = rep(NA, n_tgrid)
      for (nt in c(1:n_tgrid)) {
        t_vec = rep(tgrid[nt], n)
        obj_0t_vec = Copula.n(t_vec, PI0hat_ii[ii,], rho_d) - 
          PI0hat_ii[ii,]*FS_g01hat
        obj_1t_vec = Copula.n(t_vec, PI1hat_ii[ii,], rho_d) - 
          PI0hat_ii[ii,]*FS_g11hat
        obj_0_t[nt] = mean(obj_0t_vec^2)
        obj_1_t[nt] = mean(obj_1t_vec^2)
      }
      that_0.index = which.min(obj_0_t)
      that_0_xyr[ii,jj] = tgrid[that_0.index] ####
      that_1.index = which.min(obj_1_t)
      that_1_xyr[ii,jj] = tgrid[that_1.index] ####
      
      ## Step 4 for selection
      # i. calculate objective functions
      obj_0_rho = obj_0_rho + obj_0_t[that_0.index]
      obj_1_rho = obj_1_rho + obj_1_t[that_1.index]
      
    }
    
  }
  
  return(list(obj_0_rho, obj_1_rho, that_0_xyr, that_1_xyr))
}

## Simulation Starting...
kk.list = c(1:n_rho)
par_find_rho_Output = parLapply(cl, kk.list, par_find_rho)
stopCluster(cl)

save.image(file = "OBhimdimcfps2018Mainv2GJ_temp2.RData")

obj_0_rho_par = rep(NA, n_rho)
obj_1_rho_par = rep(NA, n_rho)
for (kk in kk.list) {
  par1output = par_find_rho_Output[[kk]]
  obj_0_rho_par[kk] = par1output[[1]]
  obj_1_rho_par[kk] = par1output[[2]]
}

# Optimal rho_d:
rho_0_index = which.min(obj_0_rho_par)
rho_0_min = rho_grid[rho_0_index]
rho_1_index = which.min(obj_1_rho_par)
rho_1_min = rho_grid[rho_1_index]

# Estimate Fstar
F.star.hat_y0 = par_find_rho_Output[[rho_0_index]][[3]]
F.star.hat_y1 = par_find_rho_Output[[rho_1_index]][[4]]

rm(par_find_rho_Output, par1output, PI0hat_ii, PI1hat_ii)


###################################
###### Auto - Debias ##############
###################################


# Auto-Debias Parameters #
D_LB=0 #each diagonal entry of \hat{D} lower bounded by D_LB
D_add=.2 #each diagonal entry of \hat{D} increased by D_add. 0.1 for 0, 0,.2 otw
max_iter=10 #max number iterations in Dantzig selector iteration over estimation and weights

#p0=dim(X0) used in low-dim dictionary in the stage 1 tuning procedure
p0=ceiling(p/4) #p/2 works for low p
if (p>60){
  p0=ceiling(p/40)
}

# specifications
delta = 1  ## x.up = D*b(x), x.down = 0, x = (1-D)*b(x)

# Auto-Debias Riesz representor #
# without selection
alpha.aml = rep(NA, n)
ratio.aml = rep(NA, n)
beta.xz = rep(NA, n)
beta.aml = rep(NA, n)

# with selection
ratio1.S.aml = rep(NA, n)
alpha1.S.aml = rep(NA, n)
ratio0.S.aml = rep(NA, n)
alpha0.S.aml = rep(NA, n)
ratio00.S.aml = rep(NA, n)
alpha00.S.aml = rep(NA, n)

# sample splitting
ylim.ad = median(y)
for(k in 1:f){
  cat('fold_auto =', k, '\n')
  foldind = datafolds.list[[k]]
  D.nfold = D[-foldind]
  s.nfold = s[-foldind]
  z.nfold = z[-foldind]
  xpoly.nfold = xpoly[-foldind,]
  xzpoly.nfold = xzpoly[-foldind,]
  y.nfold = y[-foldind]
  yvar.ad.nfold = (y.nfold==ylim.list[3])
  
  ## Auto-ML method: DML without Selection (Chernozhukov,et al., 2022)
  xinput.nfold = (1-D.nfold)*xpoly.nfold
  x.up.nfold = D.nfold*xpoly.nfold
  x.down.nfold = matrix(0, dim(x.up.nfold)[1], dim(x.up.nfold)[2])
  rho_hat=RMD_stable(yvar.ad.nfold,xinput.nfold,x.up.nfold,x.down.nfold,
                     delta,p0,D_LB,D_add,max_iter,1,1) # lasso
  ratio.aml[foldind] = xpoly[foldind,]%*%rho_hat
  alpha.aml[foldind] = ratio.aml[foldind] * 
    (1-mean(D.nfold))/(mean(D.nfold))
  
  ## Auto-ML method: DML with selection
  # CDF when D=1 
  xinput.nfold = s.nfold*(D.nfold)*xpoly.nfold
  x.up.nfold = D.nfold*xpoly.nfold
  x.down.nfold = matrix(0, dim(x.up.nfold)[1], dim(x.up.nfold)[2])
  rho2_hat=RMD_stable(yvar.ad.nfold,xinput.nfold,x.up.nfold,x.down.nfold,
                      delta,p0,D_LB,D_add,max_iter,1,1) # lasso
  ratio1.S.aml[foldind] = xpoly[foldind,]%*%rho2_hat
  alpha1.S.aml[foldind] = ratio1.S.aml[foldind] * 
    (mean((D.nfold==1)&(s.nfold==1)))/(mean(D.nfold))
  
  # Counterfactual CDF with D=0
  xinput.nfold = s.nfold*(1-D.nfold)*xpoly.nfold
  x.up.nfold = D.nfold*xpoly.nfold
  x.down.nfold = matrix(0, dim(x.up.nfold)[1], dim(x.up.nfold)[2])
  rho3_hat=RMD_stable(yvar.ad.nfold,xinput.nfold,x.up.nfold,x.down.nfold,
                      delta,p0,D_LB,D_add,max_iter,1,1) # lasso
  ratio0.S.aml[foldind] = xpoly[foldind,]%*%rho3_hat
  alpha0.S.aml[foldind] = ratio0.S.aml[foldind] * 
    (mean((D.nfold==0)&(s.nfold==1)))/(mean(D.nfold))
  
  # CDF when D=0
  xinput.nfold = s.nfold*(1-D.nfold)*xpoly.nfold
  x.up.nfold = (1-D.nfold)*xpoly.nfold
  x.down.nfold = matrix(0, dim(x.up.nfold)[1], dim(x.up.nfold)[2])
  rho4_hat=RMD_stable(yvar.ad.nfold,xinput.nfold,x.up.nfold,x.down.nfold,
                      delta,p0,D_LB,D_add,max_iter,1,1) # lasso
  ratio00.S.aml[foldind] = xpoly[foldind,]%*%rho2_hat
  alpha00.S.aml[foldind] = ratio00.S.aml[foldind] * 
    (mean((D.nfold==1)&(s.nfold==1)))/(mean(D.nfold))
  
}


###################################
###################################
ncount = 1
ylim = ylim.list[ncount]
######
## DML estimation with or without selection
# bases functions: polynomials
xpoly = cbind(x, x.c.ep^2, w)
xzpoly = cbind(xpoly, z)

# without selection:
Fyx = rep(NA, n)


# Find rho_d first
# Step 1
s0 = s[D==0]
xz0poly = xzpoly[D==0,]

s1 = s[D==1]
xz1poly = xzpoly[D==1,]

# Step 2
xz01poly = xzpoly[(D==0)&(s==1),]
xz11poly = xzpoly[(D==1)&(s==1),]

# sample splitting
for(k in 1:f){
  cat('fold =', k, '\n')
  foldind = datafolds.list[[k]]
  D.nfold = D[-foldind]
  s.nfold = s[-foldind]
  z.nfold = z[-foldind]
  xpoly.nfold = xpoly[-foldind,]
  xzpoly.nfold = xzpoly[-foldind,]
  y.nfold = y[-foldind]
  yvar.nfold = (y.nfold==ylim)
  
  # conditional CDF
  x0poly = xpoly.nfold[as.logical(1-D.nfold),]
  x1poly = xpoly.nfold[D.nfold,]
  y0.nfold = y.nfold[as.logical(1-D.nfold)]
  yvar0.nfold = as.numeric(y0.nfold==ylim)
  cv.out = cv.glmnet(x0poly, yvar0.nfold, type.measure = 'mse',
                     family = "binomial")
  y0coeff = coefficients(cv.out, s = cv.out$lambda.min)
  nz_index0 = y0coeff@i[-1] + 1
  Fyx.index= cbind(0, xpoly[foldind,nz_index0-1]) %*%
    c(0, y0coeff[nz_index0]) + y0coeff[1]
  Fyx[foldind] = exp(Fyx.index)/(1+exp(Fyx.index))
  
}

## Estimate Counterfactual CDF (unconditional)
# without selection:
# 1. Y1 cdf:
yvar = as.numeric(y==ylim)
F1yhat = mean(yvar[D==1])

# MBootstrap:
F1yhat.B = (b.xi[,D==1]%*%yvar[D==1])/sum(D)
sd.y1.aml = sd(F1yhat.B)

# 2. Counterfactual Yc cdf:
yvar0 = yvar[D==0]
Fyx.c = Fyx[D]
Fyx0 = Fyx[as.logical(1-D)]
Fhatyx.c.aml = mean(Fyx.c) + mean(alpha.aml[as.logical(1-D)] * (yvar0-Fyx0))

# MBootstrap:
Fhatyx.c.aml.B = (b.xi[,D==1]%*%Fyx.c)/sum(D) + (b.xi[,D==0]%*%(alpha.aml[as.logical(1-D)] * 
                                                                  (yvar0-Fyx0)))/(n-sum(D))
sd.yc.aml = sd(Fhatyx.c.aml.B)

# 3. Y0 cdf
F0yhat = mean(yvar[D==0])

# MBootstrap:
F0yhat.B = (b.xi[,D==0]%*%yvar[D==0])/sum(D==0)
sd.y0.aml = sd(F0yhat.B)

# 4. Stru. Effect (CDF)
StEffectyx.c.aml = F1yhat - Fhatyx.c.aml
CompEffectyx.c.aml = Fhatyx.c.aml - F0yhat

# MBootstrap:
StEffectyx.c.aml.B = F1yhat.B - Fhatyx.c.aml.B
CompEffectyx.c.aml.B = Fhatyx.c.aml.B - F0yhat.B

# with selection:
# 1. Y1 cdf:
F1yhat_star = mean(as.vector(F.star.hat_y1[D==1,ncount])) + 
  mean(alpha1.S.aml[(D==1)&(s==1)]*(
    yvar[(D==1)&(s==1)] - as.vector(F.star.hat_y1[(D==1)&(s==1),ncount])))

# MBootstrap:
F1yhat_star.B = (b.xi[,D==1]%*%as.vector(F.star.hat_y1[D==1,ncount]))/sum(D) + 
  (b.xi[,(D==1)&(s==1)]%*%(alpha1.S.aml[(D==1)&(s==1)]*(
    yvar[(D==1)&(s==1)] - as.vector(F.star.hat_y1[(D==1)&(s==1),ncount]))))/sum(((D==1)&(s==1)))
sd.y1.aml_star = sd(F1yhat_star.B)

# 2. Counterfactual Yc cdf:
Fhatyx.c.aml_star = mean(as.vector(F.star.hat_y0[D==1,ncount])) + 
  mean(alpha0.S.aml[(D==0)&(s==1)]*(yvar[(D==0)&(s==1)] - 
                                      as.vector(F.star.hat_y0[(D==0)&(s==1),ncount])))

# MBootstrap:
Fhatyx.c.aml_star.B = (b.xi[,D==1]%*%as.vector(F.star.hat_y0[D==1,ncount]))/sum(D) + 
  (b.xi[,(D==0)&(s==1)]%*%(alpha0.S.aml[(D==0)&(s==1)]*(yvar[(D==0)&(s==1)] - 
                                                          as.vector(F.star.hat_y0[(D==0)&(s==1),ncount]))))/sum(((D==0)&(s==1)))
sd.yc.aml_star = sd(Fhatyx.c.aml_star.B)


# 3.Y0 cdf
F0yhat_star = mean(as.vector(F.star.hat_y0[D==0,ncount])) + 
  mean(alpha00.S.aml[(D==0)&(s==1)]*(
    yvar[(D==0)&(s==1)] - as.vector(F.star.hat_y0[(D==0)&(s==1),ncount])))

# MBootstrap
F0yhat_star.B = (b.xi[,D==0]%*%as.vector(F.star.hat_y0[D==0,ncount]))/sum(D==0) + 
  (b.xi[,(D==0)&(s==1)]%*%(alpha00.S.aml[(D==0)&(s==1)]*(
    yvar[(D==0)&(s==1)] - as.vector(F.star.hat_y0[(D==0)&(s==1),ncount]))))/sum(((D==0)&(s==1)))
sd.y0.aml_star = sd(F0yhat_star.B)

# 4. Stru. Effect (CDF)
StEffectyx.c.aml_star = F1yhat_star - Fhatyx.c.aml_star
CompEffectyx.c.aml_star = Fhatyx.c.aml_star - F0yhat_star
StEffectyx.c.aml_star.B = F1yhat_star.B - Fhatyx.c.aml_star.B
CompEffectyx.c.aml_star.B = Fhatyx.c.aml_star.B - F0yhat_star.B

## B-standard deviation
sd.ste.aml = sd(as.vector(StEffectyx.c.aml.B)) # without selection
sd.cpe.aml = sd(as.vector(CompEffectyx.c.aml.B))
sd.ste.aml_star = sd(as.vector(StEffectyx.c.aml_star.B)) # with selection
sd.cpe.aml_star = sd(as.vector(CompEffectyx.c.aml_star.B))
# save.image(file = "OBhimdimcfps2018Mainv2_z1z2_GJ.RData")
save.image(file = "OBhimdimcfps2018Mainv2_z1_GJ.RData")


## Output
# without selection
F1yhat
sd.y1.aml
F0yhat
sd.y0.aml
Fhatyx.c.aml
sd.yc.aml
StEffectyx.c.aml
sd.ste.aml
CompEffectyx.c.aml
sd.cpe.aml

# with selection
F1yhat_star
sd.y1.aml_star
F0yhat_star
sd.y0.aml_star
Fhatyx.c.aml_star
sd.yc.aml_star
StEffectyx.c.aml_star
sd.ste.aml_star
CompEffectyx.c.aml_star
sd.cpe.aml_star