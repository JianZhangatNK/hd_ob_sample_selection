##### Estimate Structure Effects (simulation version) ###
##### Model: High dimensional with Sample Selection ####
#####        OB decomposition                      #####
##### Ver. 4.1 ########################################
#### Features #######################################
#### Auto ML for DML and PSDB ############################
#####################################################

rm(list = ls())
library('Matrix')
library('foreach')
library('VGAM')


## Auto ML(Chernozhukov,et al., 2022) packages ##
library("foreign")
library("dplyr")
library("ggplot2")
library("quantreg")

library("MASS")
library("glmnet")	#lasso, group lasso, and ridge, for outcome models, LPM, mlogit. Also, unpenalized mlogit via glmnet is far faster than the mlogit package.
library("grplasso")	#glmnet group lasso requires the same n per group (ie multi-task learning), which is perfect for mlogit but wrong for the outcome model.
# library("mlogit")	#slower than unpenalized estimation in glmnet, but glmnet won't fit only an intercept
library("nnet")	#quicker multinomial logit
library("randomForest")
library("gglasso")
library("plotrix")
library("gridExtra")

## Parallel computation
library('parallel')
library('doSNOW')
library('snow')
##################################
# Our Functions:
source('source_hdobss.R')

# Auto Debias ML (chernozhukov et al., 2022)
source('primitives.R')
source('stage1v3.R')

###### Parameters ######
n = 500 
p = 50 # p = dim(x)
SIM = 100
ylim = 1 #0 0.5 1

# Treatment Effect params
cx = 1
cz = 0.5
cy = 0.5
Sy = 0.2

# Copula param
rho_0 = 0.5
rho_1 = 0.5

# Grid points:
ngrid = 3 # y
tgrid = seq(0.01, 0.99, 0.01) # t=F*(y|x)=G(x*beta(y))
n_tgrid = length(tgrid)
n_rho = 19 # rho^d  copula coefficients
rho_grid = seq(from = -0.95, to = 0.95, length.out = n_rho)

#### Auto ML params(Chernozhukov,et al.,2022) #####
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

# ~ functions
source('primitives.R')
source('stage1.R')

################################################

# cov matrix
rho = 0.6
x.var.mat = rho.mat(p, rho)

# coefficients
b.coeff = 0.8
b.vec = b.coeff^(c(1:(p-3)))
bsq = b.coeff^2
sd.b = bsq*(1-bsq^(p-3))/(1-bsq)

# Create Folds:
f = 10
datafolds = matrix(c(1:n), f, n/f)

###### Output/Variables ####
F1yhat = rep(NA, SIM)
F1yhat_star = rep(NA, SIM)

Fhatyx.c.aml = rep(NA, SIM)
Fhatyx.c.aml_star = rep(NA, SIM)

StEffectyx.c.aml = rep(NA, SIM)
StEffectyx.c.aml_star = rep(NA, SIM)

sd.ste.aml = rep(NA, SIM)
sd.ste.aml_star = rep(NA, SIM)


set.seed(5201204)
true_value = true_value_hdob(ylim, cx, cz, cy, sy, x.var.mat, rho_0, rho_1, b.vec)
for (sim in 1:SIM) {
  cat('sim =', sim, '\n')
  ## Data generating
  x0 = mvrnorm(n, rep(0,p), x.var.mat) #z = replicate(p, rnorm(n, 0, 1))
  xprime = rnorm(n, cz, 1)
  D = (runif(n, 0, 1)>0.5)
  # Treatment on x:
  x = x0
  x1 = x0
  x1[,1] = xprime
  x[,1] = (1-D)*x[,1] + D*xprime
  u0 = rnorm(n, 0, 1)
  u1 = rnorm(n, 0, 1)
  u = (1-D)*u0 + D*u1
  # latent outcome ystar:
  ystar = as.vector(x[,1]+x[,-(1:3)]%*%b.vec) + cy*D + u 
  
  # selection:
  z = rnorm(n, 0, 1)
  v0 = sqrt(1-rho_0^2)*rnorm(n, 0, 1) + rho_0*u0
  v1 = sqrt(1-rho_1^2)*rnorm(n, 0, 1) + rho_1*u1
  s0 = (1+0.5*x[,2] + z>v0)
  s1 = (2+0.25*x[,2] + z>v1)
  s = (1-D)*s0 + D*s1
  
  # observed outcome y:
  y = ystar * s
  
  ## End of Data generating
  
  ## DML estimation with or without selection
  # bases functions: polynomials
  xpoly = cbind(x, x^2)
  xzpoly = cbind(xpoly, z, z^2)
  xzdpoly = cbind(xzpoly, D*xzpoly)
  
  # without selection:
  alpha.aml = rep(NA, n)
  ratio.aml = rep(NA, n)
  beta.xz = rep(NA, n)
  beta.aml = rep(NA, n)
  Fyx = rep(NA, n)
  
  
  # with selection
  Fstar_hat_yx0 = rep(NA, n)
  Fstar_hat_yx1 = rep(NA, n)
  ratio1.S.aml = rep(NA, n)
  alpha1.S.aml = rep(NA, n)
  ratio0.S.aml = rep(NA, n)
  alpha0.S.aml = rep(NA, n)
  
  # y evaluating points:
  ylim.list = quantile(y, c(1:ngrid)/(1+ngrid))
  
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
  
  that_0_xyr = array(rep(NA, n*n_rho*ngrid), dim = c(n, ngrid, n_rho)) 
  that_1_xyr = array(rep(NA, n*n_rho*ngrid), dim = c(n, ngrid, n_rho)) 
  obj_0_rho = rep(0, n_rho)
  obj_1_rho = rep(0, n_rho)
  for (kk in c(1:n_rho)) {
    cat('loop1 =', kk, '\n')
    rho_d = rho_grid[kk]
    
    for (jj in c(1:ngrid)) {
      y_grid = ylim.list[jj]
      # d=0
      yvar_g = (y<=y_grid) 
      yvar_g01 = yvar_g[(D==0)&(s==1)]
      y_g01.cv.out = cv.glmnet(xz01poly, yvar_g01, type.measure = 'mse')
      y_g01coeff = coefficients(y_g01.cv.out, s = y_g01.cv.out$lambda.min)
      y_g01_index0 = y_g01coeff@i[-1] + 1
      
      
      # d=1
      yvar_g11 = yvar_g[(D==1)&(s==1)]
      y_g11.cv.out = cv.glmnet(xz11poly, yvar_g11, type.measure = 'mse')
      y_g11coeff = coefficients(y_g11.cv.out, s = y_g11.cv.out$lambda.min)
      y_g11_index0 = y_g11coeff@i[-1] + 1
      
      
      for (ii in c(1:n)) {
        xi_vec = matrix(rep(x[ii,], each = n), n, p)
        xi_zpoly = cbind(xi_vec, xi_vec^2, z, z^2)
        
        # d = 0
        PI0_ii.index= cbind(0, xi_zpoly[,s0_index0-1]) %*%
          c(0, s0coeff[s0_index0]) + s0coeff[1]
        PI0hat_ii = as.vector(exp(PI0_ii.index)/(1+exp(PI0_ii.index))) # Pi_d(x,Z)
        
        FS_g01.index= cbind(0, xi_zpoly[,y_g01_index0-1]) %*%
          c(0, y_g01coeff[y_g01_index0]) + y_g01coeff[1]
        FS_g01hat = FS_g01.index # F^s(y|x,z)
        
        
        # d = 1
        PI1_ii.index= cbind(0, xi_zpoly[,s1_index0-1]) %*%
          c(0, s1coeff[s1_index0]) + s1coeff[1]
        PI1hat_ii = as.vector(exp(PI1_ii.index)/(1+exp(PI1_ii.index))) # Pi_d(x,Z)
        
        FS_g11.index= cbind(0, xi_zpoly[,y_g11_index0-1]) %*%
          c(0, y_g11coeff[y_g11_index0]) + y_g11coeff[1]
        FS_g11hat = FS_g11.index # F^s(y|x,z)
        
        # Step 3
        obj_0_t = rep(NA, n_tgrid)
        obj_1_t = rep(NA, n_tgrid)
        for (nt in c(1:n_tgrid)) {
          t_vec = rep(tgrid[nt], n)
          obj_0t_vec = Copula.n(t_vec, PI0hat_ii, rho_d) - PI0hat_ii*FS_g01hat
          obj_1t_vec = Copula.n(t_vec, PI1hat_ii, rho_d) - PI0hat_ii*FS_g11hat
          obj_0_t[nt] = mean(obj_0t_vec^2)
          obj_1_t[nt] = mean(obj_1t_vec^2)
        }
        that_0.index = which.min(obj_0_t)
        that_0_xyr[ii,jj,kk] = tgrid[that_0.index] ####
        that_1.index = which.min(obj_1_t)
        that_1_xyr[ii,jj,kk] = tgrid[that_1.index] ####
        
        ## Step 4 for selection
        # i. calculate objective functions
        obj_0_rho[kk] = obj_0_rho[kk] + obj_0_t[that_0.index]
        obj_1_rho[kk] = obj_1_rho[kk] + obj_1_t[that_1.index]
        
      }
      
    }
  }
  
  # Optimal rho_d:
  rho_0_index = which.min(obj_0_rho)
  rho_0_min = rho_grid[rho_0_index]
  rho_1_index = which.min(obj_1_rho)
  rho_1_min = rho_grid[rho_1_index]
  
  # sample splitting
  for(k in 1:f){
    cat('fold =', k, '\n')
    foldind = datafolds[k,]
    D.nfold = D[-foldind]
    s.nfold = s[-foldind]
    z.nfold = z[-foldind]
    xpoly.nfold = xpoly[-foldind,]
    xzpoly.nfold = xzpoly[-foldind,]
    y.nfold = y[-foldind]
    yvar.nfold = (y.nfold<=ylim)
    
    ## Step 1 for Selection
    # d=0 
    s0.nfold = s.nfold[D.nfold==0]
    xz0poly.nfold = xzpoly.nfold[D.nfold==0,]
    s0.cv.out = cv.glmnet(xz0poly.nfold, s0.nfold, type.measure = 'mse', family = "binomial")
    s0coeff = coefficients(s0.cv.out, s = s0.cv.out$lambda.min)
    s0_index0 = s0coeff@i[-1] + 1
    PI0.index= cbind(0, xzpoly.nfold[,s0_index0-1]) %*%
      c(0, s0coeff[s0_index0]) + s0coeff[1]
    PI0hat = exp(PI0.index)/(1+exp(PI0.index))
    
    # d=1
    s1.nfold = s.nfold[D.nfold==1]
    xz1poly.nfold = xzpoly.nfold[D.nfold==1,]
    s1.cv.out = cv.glmnet(xz1poly.nfold, s1.nfold, type.measure = 'mse', family = "binomial")
    s1coeff = coefficients(s1.cv.out, s = s1.cv.out$lambda.min)
    s1_index0 = s1coeff@i[-1] + 1
    PI1.index= cbind(0, xzpoly.nfold[,s1_index0-1]) %*%
      c(0, s1coeff[s1_index0]) + s1coeff[1]
    PI1hat = exp(PI1.index)/(1+exp(PI1.index))
    
    ## Step 2 for Selection (script: ds)
    # d=0
    yvar01.nfold = yvar.nfold[(D.nfold==0)&(s.nfold==1)]
    xz01poly.nfold = xzpoly.nfold[(D.nfold==0)&(s.nfold==1),]
    y01.cv.out = cv.glmnet(xz01poly.nfold, yvar01.nfold, type.measure = 'mse', family = "binomial")
    y01coeff = coefficients(y01.cv.out, s = y01.cv.out$lambda.min)
    y01_index0 = y01coeff@i[-1] + 1
    FS01.index= cbind(0, xzpoly.nfold[,y01_index0-1]) %*%
      c(0, y01coeff[y01_index0]) + y01coeff[1]
    FS01hat = exp(FS01.index)/(1+exp(FS01.index))
    
    # d=1
    yvar11.nfold = yvar.nfold[(D.nfold==1)&(s.nfold==1)]
    xz11poly.nfold = xzpoly.nfold[(D.nfold==1)&(s.nfold==1),]
    y11.cv.out = cv.glmnet(xz11poly.nfold, yvar11.nfold, type.measure = 'mse', family = "binomial")
    y11coeff = coefficients(y11.cv.out, s = y11.cv.out$lambda.min)
    y11_index0 = y11coeff@i[-1] + 1
    FS11.index= cbind(0, xzpoly.nfold[,y11_index0-1]) %*%
      c(0, y11coeff[y11_index0]) + y11coeff[1]
    FS11hat = exp(FS11.index)/(1+exp(FS11.index))
    
    ## Step 3 for Selection
    l.fold = length(foldind)
    that_0_xyr = array(rep(NA, l.fold*n_rho*ngrid), dim = c(l.fold, ngrid, n_rho)) 
    that_1_xyr = array(rep(NA, l.fold*n_rho*ngrid), dim = c(l.fold, ngrid, n_rho)) 
    that2_0_xr = matrix(data = NA, l.fold, n_rho)
    that2_1_xr = matrix(data = NA, l.fold, n_rho)
    obj_0_rho = matrix(0, l.fold, n_rho)
    obj_1_rho = matrix(0, l.fold, n_rho)
    for (kk in c(1:n_rho)) {
      cat('loop2 =', kk, '\n')
      rho_d = rho_grid[kk]
      for (ii in c(1:l.fold)) {
        xi_vec = matrix(rep(x[ii,], each = n - l.fold), n - l.fold, p)
        xi_zpoly = cbind(xi_vec, xi_vec^2, z.nfold, z.nfold^2)
        
        # d = 0
        PI0_ii.index= cbind(0, xi_zpoly[,s0_index0-1]) %*%
          c(0, s0coeff[s0_index0]) + s0coeff[1]
        PI0hat_ii = as.vector(exp(PI0_ii.index)/(1+exp(PI0_ii.index))) # Pi_d(x,Z)
        FS01_ii.index= cbind(0, xi_zpoly[,y01_index0-1]) %*%
          c(0, y01coeff[y01_index0]) + y01coeff[1]
        FS01hat_ii = exp(FS01_ii.index)/(1+exp(FS01_ii.index))
        obj2_0_t = rep(NA, n_tgrid)
        for (nt in c(1:n_tgrid)) {
          t_vec = rep(tgrid[nt], n - l.fold)
          obj2_t_vec = Copula.n(t_vec, PI0hat_ii, rho_d) - PI0hat_ii*FS01hat_ii
          obj2_0_t[nt] = mean(obj2_t_vec^2)
        }
        that2_0.index = which.min(obj2_0_t)
        that2_0_xr[ii,kk] = tgrid[that2_0.index] ####
        
        # d = 1
        PI1_ii.index= cbind(0, xi_zpoly[,s1_index0-1]) %*%
          c(0, s1coeff[s1_index0]) + s1coeff[1]
        PI1hat_ii = as.vector(exp(PI1_ii.index)/(1+exp(PI1_ii.index))) # Pi_d(x,Z)
        FS11_ii.index= cbind(0, xi_zpoly[,y11_index0-1]) %*%
          c(0, y11coeff[y11_index0]) + y11coeff[1]
        FS11hat_ii = exp(FS11_ii.index)/(1+exp(FS11_ii.index))
        obj2_1_t = rep(NA, n_tgrid)
        for (nt in c(1:n_tgrid)) {
          t_vec = rep(tgrid[nt], n - l.fold)
          obj2_t_vec = Copula.n(t_vec, PI1hat_ii, rho_d) - PI1hat_ii*FS11hat_ii
          obj2_1_t[nt] = mean(obj2_t_vec^2)
        }
        that2_1.index = which.min(obj2_1_t)
        that2_1_xr[ii,kk] = tgrid[that2_1.index] ####
        
      }
    }
    

    
    # Estimate F*(y|x) for d=0 or 1
    Fstar_hat_yx0_loop = rep(NA, l.fold)
    Fstar_hat_yx1_loop = rep(NA, l.fold)
    for (ii in c(1:l.fold)) {
      
      # F*(y|x):
      Fstar_hat_yx0_loop[ii] = that2_0_xr[ii, rho_0_index ]
      Fstar_hat_yx1_loop[ii] = that2_1_xr[ii, rho_1_index ]
    }
    Fstar_hat_yx0[foldind] = Fstar_hat_yx0_loop
    Fstar_hat_yx1[foldind] = Fstar_hat_yx1_loop
    
    
    ## Auto-ML method: DML without Selection (Chernozhukov,et al., 2022)
    xinput.nfold = (1-D.nfold)*xpoly.nfold
    x.up.nfold = D.nfold*xpoly.nfold
    x.down.nfold = matrix(0, dim(x.up.nfold)[1], dim(x.up.nfold)[2])
    rho_hat=RMD_stable(yvar.nfold,xinput.nfold,x.up.nfold,x.down.nfold,
                       delta,p0,D_LB,D_add,max_iter,1,1) # lasso
    ratio.aml[foldind] = xpoly[foldind,]%*%rho_hat
    alpha.aml[foldind] = ratio.aml[foldind] * 
      (1-mean(D.nfold))/(mean(D.nfold))
    
    # conditional CDF
    x0poly = xpoly.nfold[as.logical(1-D.nfold),]
    x1poly = xpoly.nfold[D.nfold,]
    y0.nfold = y.nfold[as.logical(1-D.nfold)]
    yvar0.nfold = as.numeric(y0.nfold<=ylim)
    cv.out = cv.glmnet(x0poly, yvar0.nfold, type.measure = 'mse',
                       family = "binomial")
    y0coeff = coefficients(cv.out, s = cv.out$lambda.min)
    nz_index0 = y0coeff@i[-1] + 1
    Fyx.index= cbind(0, xpoly[foldind,nz_index0-1]) %*%
      c(0, y0coeff[nz_index0]) + y0coeff[1]
    Fyx[foldind] = exp(Fyx.index)/(1+exp(Fyx.index))
    
    ## Auto-ML method: DML with selection
    # CDF when D=1 
    xinput.nfold = s.nfold*(D.nfold)*xpoly.nfold
    x.up.nfold = D.nfold*xpoly.nfold
    x.down.nfold = matrix(0, dim(x.up.nfold)[1], dim(x.up.nfold)[2])
    rho2_hat=RMD_stable(yvar.nfold,xinput.nfold,x.up.nfold,x.down.nfold,
                        delta,p0,D_LB,D_add,max_iter,1,1) # lasso
    ratio1.S.aml[foldind] = xpoly[foldind,]%*%rho2_hat
    alpha1.S.aml[foldind] = ratio1.S.aml[foldind] * 
      (mean((D.nfold==1)&(s.nfold==1)))/(mean(D.nfold))
    
    # Counterfactual CDF with D=0
    xinput.nfold = s.nfold*(1-D.nfold)*xpoly.nfold
    x.up.nfold = D.nfold*xpoly.nfold
    x.down.nfold = matrix(0, dim(x.up.nfold)[1], dim(x.up.nfold)[2])
    rho3_hat=RMD_stable(yvar.nfold,xinput.nfold,x.up.nfold,x.down.nfold,
                        delta,p0,D_LB,D_add,max_iter,1,1) # lasso
    ratio0.S.aml[foldind] = xpoly[foldind,]%*%rho3_hat
    alpha0.S.aml[foldind] = ratio0.S.aml[foldind] * 
      (mean((D.nfold==0)&(s.nfold==1)))/(mean(D.nfold))
    
  }
  
  ## Estimate Counterfactual CDF (unconditional)
  # without selection:
  # 1. Y1 cdf:
  yvar = as.numeric(y<=ylim)
  F1yhat[sim] = mean(yvar[D==1])
  # 2. Counterfactual Yc cdf:
  yvar0 = yvar[D==0]
  Fyx.c = Fyx[D]
  Fyx0 = Fyx[as.logical(1-D)]
  Fhatyx.c.aml[sim] = mean(Fyx.c) + mean(alpha.aml[as.logical(1-D)] * (yvar0-Fyx0))
  # 3. Stru. Effect (CDF)
  StEffectyx.c.aml[sim] = F1yhat[sim] - Fhatyx.c.aml[sim]
  
  # with selection:
  # 1. Y1 cdf:
  F1yhat_star[sim] = mean(Fstar_hat_yx1[D==1]) + 
    mean(alpha1.S.aml[(D==1)&(s==1)]*(
      yvar[(D==1)&(s==1)] - Fstar_hat_yx1[(D==1)&(s==1)]))
  # 2. Counterfactual Yc cdf:
  Fhatyx.c.aml_star[sim] = mean(Fstar_hat_yx0[D==1]) + 
    mean(alpha0.S.aml[(D==0)&(s==1)]*(yvar[(D==0)&(s==1)] - 
                                        Fstar_hat_yx0[(D==0)&(s==1)] ))
  # 3. Stru. Effect (CDF)
  StEffectyx.c.aml_star[sim] = F1yhat_star[sim] - Fhatyx.c.aml_star[sim]
  
}

save.image(file = "HDOB_SS_v4_1.RData")


################################################################################
############# Generate pictures ################################################
## 图片
# Latex:
library(latex2exp)
##################################

par(
  cex.main = 2,      # 主标题
  cex.lab = 1.4,       # 坐标轴标签
  cex.axis = 1,      # 坐标轴刻度
  font.main = 1.5,       # 主标题字体
  font.lab = 1.5         # 坐标轴标签字体
)

## Counterfactural CDF
ccdf.true = true_value[2]
ccdf.den.nosel = den.fitted.plot(Fhatyx.c.aml)
ccdf.den.sel = den.fitted.plot(Fhatyx.c.aml_star)
xlim.max = quantile(c(Fhatyx.c.aml,Fhatyx.c.aml_star),0.975)
xlim.min = quantile(c(Fhatyx.c.aml,Fhatyx.c.aml_star),0.025)
ylim.max = max(c(ccdf.den.nosel, ccdf.den.sel))*1.05
ylim.min = min(c(ccdf.den.nosel, ccdf.den.sel))*0.95
plot(Fhatyx.c.aml[order(Fhatyx.c.aml)], ccdf.den.nosel[order(Fhatyx.c.aml)], 
     xlim=c(xlim.min, xlim.max), ylim=c(ylim.min, ylim.max),
     xlab = TeX('$\\hat{F}^*_{0|1}(1)$'), ylab ="拟合密度函数",
     main="反事实分布函数在y=1的估计: 样本选择vs无样本选择",
     type = 'l', pch = 4, col = "blue",lwd = 2.5)
lines(Fhatyx.c.aml_star[order(Fhatyx.c.aml_star)], ccdf.den.sel[order(Fhatyx.c.aml_star)], 
      type = 'l', pch = 4, col = "darkgreen", lwd = 2.5)
abline(v=ccdf.true, col = 'red', lwd = 2.5, lty = 'dashed')
legend("topleft", 
       legend=c("无样本选择", "样本选择", "真实值"),
       col=c( "blue", "darkgreen", "red"), 
       lty=c(1,1,2), lwd = c(2.5,2.5,2.5), cex=1.5)

## Structure Effect
StrE.true = true_value[3]
StrE.den.nosel = den.fitted.plot(StEffectyx.c.aml)
StrE.den.sel = den.fitted.plot(StEffectyx.c.aml_star)
xlim.max = 0.06 #quantile(c(StEffectyx.c.aml,StEffectyx.c.aml_star),0.975)
xlim.min = quantile(c(StEffectyx.c.aml,StEffectyx.c.aml_star),0.025)
ylim.max = max(c(StrE.den.nosel, StrE.den.sel))*1.05
ylim.min = min(c(StrE.den.nosel, StrE.den.sel))*0.95
plot(StEffectyx.c.aml[order(StEffectyx.c.aml)], StrE.den.nosel[order(StEffectyx.c.aml)], 
     xlim=c(xlim.min, xlim.max), ylim=c(ylim.min, ylim.max),
     xlab = "结构效应", ylab ="拟合密度函数",
     main="OB分解结构效应在y=1的估计: 样本选择vs无样本选择",
     type = 'l', pch = 4, col = "blue",lwd = 2.5)
lines(StEffectyx.c.aml_star[order(StEffectyx.c.aml_star)], StrE.den.sel[order(StEffectyx.c.aml_star)], 
      type = 'l', pch = 4, col = "darkgreen", lwd = 2.5)
abline(v=StrE.true, col = 'red', lwd = 2.5, lty = 'dashed')
legend("topright", 
       legend=c("无样本选择", "样本选择", "真实值"),
       col=c( "blue", "darkgreen", "red"), 
       lty=c(1,1,2), lwd = c(2.5,2.5,2.5), cex=1.5)




