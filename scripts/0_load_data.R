# ---- Data ----

retention <- read.csv('./data/retention.csv',sep=',',dec=',')

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

#sasa and egb

sasaALL <- read.csv('./data/sasa_all.csv',sep=',',dec='.', skip = 1)
egbFile <- read.csv('./data/egb.csv',sep=',',dec='.', skip = 1)

# the retention times were stored in different order than sasa and egb, therefore following loop is pairing correct trinucleotides
for(i in 1:27){
   j = which(sasaALL$seq == retention$trink[i])
   retention$egb[i] <- egbFile$Egb[j]
   retention$sasaALL[i] <- sasaALL$Total[j]
}

