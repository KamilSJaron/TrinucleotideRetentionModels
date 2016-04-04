# ---- Packages ----

library(MASS)
library(DAAG)

# ---- Data ----

retention <- read.csv('/Volumes/dump/data/chemical_properties_nucl/retention.csv',sep=',',dec=',')
chemizalize <- read.csv('/Volumes/dump/data/chemical_properties_nucl/chemicalize_propertes.csv',sep='\t')

retention$logk <- log(retention$k,10)

# inicialization and setting of dummy variables of presence of a nuclotide on a position, these will be used for testing of outliers
retention$A_1 = 0
retention$G_1 = 0
retention$T_1 = 0
retention$A_2 = 0
retention$G_2 = 0
retention$T_2 = 0
retention$A_3 = 0
retention$G_3 = 0
retention$T_3 = 0

retention$A_1[which(retention$st == 'A')] = 1
retention$T_1[which(retention$st == 'T')] = 1
retention$G_1[which(retention$st == 'G')] = 1
retention$A_2[which(retention$nd == 'A')] = 1
retention$T_2[which(retention$nd == 'T')] = 1
retention$G_2[which(retention$nd == 'G')] = 1
retention$A_3[which(retention$rd == 'A')] = 1
retention$T_3[which(retention$rd == 'T')] = 1
retention$G_3[which(retention$rd == 'G')] = 1

# iniciation and merging surface estimates of nucleorides and sugars (dR) using sasa for every nucleotide (Ns_index) and for sum values (sN)


# measurements of retention factor and sasa models are in different orders, this for loop reorganises surface data to retention table
# for(i in 1:27){
#   trink = paste(retention[i,2],retention[i,3],retention[i,4],sep = '')
#   j = which(surf2[,1] == trink)
#   retention[i,16:27] = as.numeric(as.matrix(surf2[j,c(-1)]))
#   retention$sA[i] = as.numeric(as.matrix(surf$Total.A.[j]))
#   retention$sG[i] = as.numeric(as.matrix(surf$Total.G.[j]))
#   retention$sT[i] = as.numeric(as.matrix(surf$Total.T.[j]))
#   retention$sdR[i] = as.numeric(as.matrix(surf$Total.dR.[j]))
# }

# computation of transformed varibles polsur = sum (polarity * surface) and poldivsur = sum (polarity / surface) and sum of polarities of all nucleotides
# retention$polsur <- 
#   chemizalize[as.character(retention$st),'logP'] * rowSums(retention[,c('As_1','Gs_1','Ts_1')]) +
#   chemizalize[as.character(retention$nd),'logP'] * rowSums(retention[,c('As_2','Gs_2','Ts_2')]) +
#   chemizalize[as.character(retention$rd),'logP'] * rowSums(retention[,c('As_3','Gs_3','Ts_3')])
# 
# retention$poldivsur <- 
#   chemizalize[as.character(retention$st),'logP'] / rowSums(retention[,c('As_1','Gs_1','Ts_1')]) +
#   chemizalize[as.character(retention$nd),'logP'] / rowSums(retention[,c('As_2','Gs_2','Ts_2')]) +
#   chemizalize[as.character(retention$rd),'logP'] / rowSums(retention[,c('As_3','Gs_3','Ts_3')])

retention$pol <- 
  (retention$A_1 + retention$A_2 + retention$A_3) * chemizalize['A','logP'] +
  (retention$G_1 + retention$G_2 + retention$G_3) * chemizalize['G','logP'] +
  (retention$T_1 + retention$T_2 + retention$T_3) * chemizalize['T','logP']

# retention$sur <- retention$sA + retention$sG + retention$sT

# cleaning of redundant variables (only table with polarities and retention factors are kept)
# rm(list = c('i','j','trink'))

sasaALL <- read.csv('/Volumes/dump/data/chemical_properties_nucl/sasa_all.csv',sep=',',dec='.', skip = 1)
sesaALL <- read.csv('/Volumes/dump/data/chemical_properties_nucl/sesa_all.csv',sep=',',dec='.', skip = 1)

sasaPSA <- read.csv('/Volumes/dump/data/chemical_properties_nucl/sasa_psa.csv',sep=',',dec='.', skip = 1)
sasaNSA <- read.csv('/Volumes/dump/data/chemical_properties_nucl/sasa_nsa.csv',sep=',',dec='.', skip = 1)
sesaPSA <- read.csv('/Volumes/dump/data/chemical_properties_nucl/sesa_psa.csv',sep=',',dec='.', skip = 1)
sesaNSA <- read.csv('/Volumes/dump/data/chemical_properties_nucl/sesa_nsa.csv',sep=',',dec='.', skip = 1)

egbFile <- read.csv('/Volumes/dump/data/chemical_properties_nucl/egb.csv',sep=',',dec='.', skip = 1)

retention$sasadR <- 0
retention$sesadR <- 0
retention$sasadR1 <- 0
retention$sesadR1 <- 0
retention$sasadR2 <- 0
retention$sesadR2 <- 0
retention$sasadR3 <- 0
retention$sesadR3 <- 0

retention$sasaPSAsur <- 0
retention$sasaNSAsur <- 0
retention$sasaPolSurChemi <- 0
retention$sasaPolSurSkype <- 0
logPG = -0.6
logPT = -0.4
logPA = -0.1

for(i in 1:27){
  j = which(sesaALL$seq == retention$trink[i])
  retention$sasadR[i] <- (sesaALL$R1 + sesaALL$R2 + sesaALL$R3)[j]
  retention$sesadR[i] <- (sasaALL$R1 + sasaALL$R2 + sasaALL$R3)[j]
  retention$sasadR1[i] <- (sesaALL$R1)[j]
  retention$sesadR1[i] <- (sasaALL$R1)[j]
  retention$sasadR2[i] <- (sesaALL$R2)[j]
  retention$sesadR2[i] <- (sasaALL$R2)[j]
  retention$sasadR3[i] <- (sesaALL$R3)[j]
  retention$sesadR3[i] <- (sasaALL$R3)[j]
  retention$sasaPSAsur[i] <- sasaPSA$Total[j]
  retention$sasaNSAsur[i] <- sasaNSA$Total[j]
  retention$sasaPolSurChemi[i] <-   sum(sasaPSA[j,c('B1.A','B2.A','B3.A')] * chemizalize['A','logP'] +
                                          sasaPSA[j,c('B1.G','B2.G','B3.G')] * chemizalize['G','logP'] +
                                          sasaPSA[j,c('B1.T','B2.T','B3.T')] * chemizalize['T','logP'])
  retention$sasaPolSurSkype[i] <-   sum(sasaPSA[j,c('B1.A','B2.A','B3.A')] * logPA +
                                          sasaPSA[j,c('B1.G','B2.G','B3.G')] * logPG +
                                          sasaPSA[j,c('B1.T','B2.T','B3.T')] * logPT)
  retention$egb[i] <- egbFile$Egb[j]
}

retention$polSkype <- 
  (retention$A_1 + retention$A_2 + retention$A_3) * logPA +
  (retention$G_1 + retention$G_2 + retention$G_3) * logPG +
  (retention$T_1 + retention$T_2 + retention$T_3) * logPT

retention$Pol1 <- (retention$A_1 * logPA) + (retention$T_1 * logPT) + (retention$G_1 * logPG)
retention$Pol2 <- (retention$A_2 * logPA) + (retention$T_2 * logPT) + (retention$G_2 * logPG)
retention$Pol3 <- (retention$A_3 * logPA) + (retention$T_3 * logPT) + (retention$G_3 * logPG)

colnames(sasaNSA) <- paste("NSA",colnames(sasaNSA),sep='')
colnames(sasaPSA) <- paste("PSA",colnames(sasaPSA),sep='')
sasa <- cbind(sasaNSA,sasaPSA)

colnames(sesaNSA) <- paste("NSA",colnames(sesaNSA),sep='')
colnames(sesaPSA) <- paste("PSA",colnames(sesaPSA),sep='')
sesa <- cbind(sesaNSA,sesaPSA)

retention$sasaALL <- retention$sasaPSAsur + retention$sasaNSAsur


