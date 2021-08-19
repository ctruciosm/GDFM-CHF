## Read Condicional Covariance Matrizes
# all(double(1) == 0)


library(Rsolnp)
library(data.table)
source("optim.mvp.bounded.R")
retornos = read.table("retornos_APP3.txt",sep=" ",head=TRUE)
Ht = fread("H_GDFM_CHF_DCC_NL.txt",header=FALSE)
Ht = as.matrix(Ht)



WR = dim(retornos)[1]-750
WR_BACKUP = WR
WR1 = WR_BACKUP
WR = WR1 #-89
r_RM_c = rep(0, WR)
cons_mvp_H = matrix(NA,ncol=656,nrow=WR)

for (i in 1:WR){ 
  print(i)
  H_aux = matrix(Ht[i,],ncol=656)  
  # if comes from matlab matrix(Ht[i,],ncol=656)
  # if comes from R matrix(Ht[,i],ncol=656)
  H_aux = 0.5*(H_aux+t(H_aux))
  Y_rea = as.matrix(retornos[i+750,]) 
  cons_mvp_H[i,] = optim.mvp(H_aux)
  r_RM_c[i] = sum(cons_mvp_H[i,]*Y_rea)
}

INI =   1 
AV_RM = mean(r_RM_c[c(INI:WR)])*252
SD_RM = sd(r_RM_c[c(INI:WR)])*sqrt(252)
IR_RM = AV_RM/SD_RM
S_RM  = ifelse(r_RM_c[c(INI:WR)]<0,r_RM_c[c(INI:WR)]^2,0)
RO_RM  = AV_RM/sqrt(252*mean(S_RM))
ResultsAnn = matrix(c(AV_RM,SD_RM,IR_RM,RO_RM))
round(ResultsAnn,4)

write.csv(r_RM_c,"OoSp_GDFM_CHF_DCC_NL.csv")


GDFM_CHF_DCC_NL <- read.csv("OoSp_GDFM_CHF_DCC_NL.csv")[,-1]
PCA_MGARCH_DCC_N <- read.csv("OoSp_H_PCA_MGARCH_DCC_NL.csv")[,-1]

OoS_GDFM_PCA <- data.frame(GDFM = GDFM_CHF_DCC_NL, PCA = PCA_MGARCH_DCC_N)
write.csv(OoS_GDFM_PCA,"OoS_GDFM_PCA.csv")

GDFM_CHF_DCC_NL <- read.csv("OoSp_GDFM_CHF_DCC_NL.csv")[,-1]
DCC_NL <- read.csv("OoSp_DCC_NL.csv")[,-1]

OoS_GDFM_DCCNL <- data.frame(GDFM = GDFM_CHF_DCC_NL, DCCNL = DCC_NL)
write.csv(OoS_GDFM_DCCNL,"OoS_GDFM_DCCNL.csv")

GDFM_CHF_DCC_NL <- read.csv("OoSp_GDFM_CHF_DCC_NL.csv")[,-1]
AFM_DCC_NL <- read.csv("OoSp_AFM_DCC_NL.csv")[,-1]

OoS_GDFM_AFMDCCNL <- data.frame(GDFM = GDFM_CHF_DCC_NL, AFMDCCNL = AFM_DCC_NL)
write.csv(OoS_GDFM_AFMDCCNL,"OoS_GDFM_AFMDCCNL.csv")




