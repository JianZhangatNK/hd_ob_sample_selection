# data loading and clearing
# Version 7: CFPS2018 data
# 9_2 Good Job case

rm(list=ls())
library("Matrix")
library(haven)
library(fastDummies)

data <- read.csv(file = "HDOBcfps2018_dataV2.csv", header = TRUE,  stringsAsFactors=FALSE)
#data<-read_dta("working_cfps2018.dta")
#data <- read.csv(file = "HDOBcfps2018_data.csv", header = TRUE,  stringsAsFactors=FALSE)
#data = as.numeric(data)
#data$hukou<-as.numeric(data$hukou)
#dummy_cols(data, select_columns =c("urban", "gender", "hukou", "spouse", "edu", "military1", "party1", "evsmoke1",  "agri_job"))

y1.raw <- data$lnwage
y2.raw <- data$week_whours
y3.raw <- data$employ
y4.raw <- data$health
y5.raw <- (data$stateown==1)*(data$wcollar==1) # good job indicator
D1.raw <- data$hukou
D2.raw <- data$gender
D3.raw <- data$formal
D4.raw <- data$minority
z1.raw <- log(data$gsubsidy_amount)
z1.raw[data$gsubsidy_amount==0] = 0
z2.raw <- data$lnsubsidy

# remember to delete treatment variable in x.d.raw
x.d.raw <- cbind( data$urban, data$familysize, data$spouse, data$military1, data$party1, data$evsmoke1, data$health, data$minority,
                  data$agri_job,  data$formal,
                  data$part_time,  data$service,  data$BJ, data$TJ, data$HeB, data$ShX, data$NMG, data$LN, data$JN, data$HLJ, data$SH, data$JS, data$ZJ, data$AnH, data$FJ, 
                  data$JX, data$ShD, data$HeN, data$HuB, data$HuN, data$GuD, data$GuX, data$HaN, data$ChQ, data$SiCh, data$GuZ, data$YunN, data$ShaX, data$GS, data$QinH, 
                  data$NinX, data$XinJ) # 
x.c.raw <- cbind(data$age, data$edu, data$experience) 
w_ini.raw <- x.d.raw


#To be continued
# delete NA rows
data.raw = cbind(y1.raw, y2.raw, y3.raw, y4.raw, y5.raw, D1.raw, D2.raw, x.c.raw, x.d.raw, z1.raw, z2.raw)
w_ini.na = is.na(data.raw)
row.na.count = rowSums(w_ini.na)
w_ini.row = which(row.na.count!=0)
cat(length(w_ini.row), 'NA rows are deleted.', '\n')
if(length(w_ini.row)==0){
  w_ini.row = 10^6
}

w_ini.nan = w_ini.raw[-w_ini.row,]
x.d = x.d.raw[-w_ini.row,]
x.c = x.c.raw[-w_ini.row,]
y1 = y1.raw[-w_ini.row]
y2 = y2.raw[-w_ini.row]
y3 = y3.raw[-w_ini.row]
y4 = y4.raw[-w_ini.row]
y5 = y5.raw[-w_ini.row]
D1 = D1.raw[-w_ini.row]
D2 = D2.raw[-w_ini.row]
D3 = D3.raw[-w_ini.row]
D4 = D4.raw[-w_ini.row]
z1 = z1.raw[-w_ini.row]
z2 = z2.raw[-w_ini.row]





# Sparse Matrix
w_ini <- Matrix(as.matrix(w_ini.nan), sparse = TRUE)
w_ini.label <- c('urban', 'familysize', 'spouse', 'military1', 'party1', 'evsmoke1', 
                 'health', 'minority', 'agri_job', 
                  'formal', 'part_time',  'service',  'BJ',
                 'TJ', 'HeB', 'ShX', 'NMG', 'LN', 'JN', 'HLJ', 'SH', 'JS', 'ZJ', 'AnH', '
                 FJ', 'JX', 'ShD', 'HeN', 'HuB', 'HuN', 'GuD', 'GuX', 'HaN', 'ChQ', 'SiCh',
                 'GuZ', 'YunN', 'ShaX', 'GS', 'QinH', 'NinX', 'XinJ')


# create cross terms
w.cross.a = 1
w.cross.a.lab = c('const')
w.len = dim(w_ini)[2]
for (i in 1:(w.len-1)) {
  for(j in (i+1):w.len){
    if(i!=j){
      d = w_ini[,i]*w_ini[,j]
      d.sum = sum(d)
      if(d.sum>0){
        w.cross.a = cbind(w.cross.a, d)
        d.lab = paste(w_ini.label[i], w_ini.label[j], sep = '_')
        w.cross.a.lab = cbind(w.cross.a.lab, d.lab)
      }
    }
  }
}
w.d.cross = w.cross.a[,-1]
w.d.cross.lab = w.cross.a.lab[-1]
x.zl = cbind(x.c, x.d)
w.cross = cbind(x.c[,1]*x.d, x.c[,2]*x.d, x.c[,3]*x.d, w.d.cross)

# continuous cross dummy labels
cd.label.a = c('delete')
cont.label = c('age', 'edu', 'expr')
cont2.label = c('age2', 'edu2', 'expr2')
xd.len = length(w_ini.label)
for (i in 1:3) {
  for (j in 1:xd.len) {
    cd.lab.app = paste(cont.label[i], w_ini.label[j], sep = '_')
    cd.label.a = cbind(cd.label.a, cd.lab.app)
  }
}
cd.label = cd.label.a[-1]
OBhim.label = c(cont.label, w_ini.label, cont2.label, cd.label, w.d.cross.lab)

save(y1, y2, y3,y4,y5, D1, D2, D3, D4, x.d, x.c, x.zl, z1, z2, w.d.cross, w.cross, file = "OBhimdimcfps2018v2GJ.RData")
save(OBhim.label, file = "OBhimdimcfps2018LBv2GJ.RData")
save.image(file="OBhimdimcfps2018v2GJ_all.RData")