#########################################################
#####   POET #######
#########################################################
# Loading packages
library(POET)
library(data.table)

# Reading data
Y =  fread("/Users/ctruciosm/Desktop/GDFM-CHF/Data/retornos_APP3.txt")[,-1]
n = nrow(Y)
p = ncol(Y)

#Setting parameters for the RW scheme.
W =  750 
L =  n  
WR =  L-750

H = matrix(0, ncol = p, nrow = p)
Hvector =  matrix(NA,nrow = p*p, ncol = WR)

for (l in 1:WR){
  print(l)
  demean_Y =  t(scale(Y[l:(l+750-1),],scale = FALSE))
  H = POET(demean_Y, K = 8, C = 0.5)$SigmaY
  Hvector[,l] =  as.vector(H)
}
Sigma <- Hvector[,1:WR]
write.table(data.frame(Sigma),"H_POET.txt",sep=",",dec=".",col.names=FALSE,row.names = FALSE)

