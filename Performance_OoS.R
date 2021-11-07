## Read Condicional Covariance Matrizes
# all(double(1) == 0)


library(Rsolnp)
library(data.table)
source("/home/ctruciosm/GDFM-CHF/aux_functions/optim.mvp.bounded.R")
retornos = read.table("/home/ctruciosm/GDFM-CHF/Data/retornos_APP3.txt",sep = " ",head = TRUE)
Ht = fread("H_GDFM_CHF_DCC_NL_normal_bic.txt", header = FALSE)
Ht = as.matrix(Ht)



WR = dim(retornos)[1] - 750
WR_BACKUP = WR
WR1 = WR_BACKUP
WR = WR1 
r_RM_c = rep(0, WR)
cons_mvp_H = matrix(NA, ncol = 656, nrow = WR)

for (i in 1:WR) { 
  print(i)
  H_aux = matrix(Ht[i,],ncol = 656)  
  # If comes from matlab use matrix(Ht[i,], ncol = 656)
  # If comes from R use matrix(Ht[,i], ncol = 656)
  H_aux = 0.5*(H_aux + t(H_aux))
  Y_rea = as.matrix(retornos[i + 750,]) 
  cons_mvp_H[i,] = optim.mvp(H_aux)
  r_RM_c[i] = sum(cons_mvp_H[i,]*Y_rea)
}

INI =   1 
AV_RM = mean(r_RM_c[c(INI:WR)])*252
SD_RM = sd(r_RM_c[c(INI:WR)])*sqrt(252)
IR_RM = AV_RM/SD_RM
S_RM  = ifelse(r_RM_c[c(INI:WR)] < 0,r_RM_c[c(INI:WR)]^2,0)
RO_RM  = AV_RM/sqrt(252*mean(S_RM))
ResultsAnn = matrix(c(AV_RM,SD_RM,IR_RM,RO_RM))
round(ResultsAnn,4)

write.csv(r_RM_c,"OoSp_H_GDFM_CHF_DCC_NL_normal_bic.csv")
