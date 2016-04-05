source('./scripts/test_trink.R')
# ======= NNN model =======
# log(k) = a*E_gb + b*SASA + c

# setting for crossvalidation: 27 = leave-one out, changing this to lower number will cause different type of cv
m_for_cv = 27 

model_formula <- "logk ~ egb + sasaALL"
model1 <- lm(model_formula, data=retention)
kfCV <- cv.lm(retention, model1, m = 27, plotit = F)

MSE_tab <- data.frame(nk = 'intercept', MSE = mean(model1$residuals^2), cvMSE = mean((kfCV$cvpred - kfCV$logk)^2), Adj.R2 = summary(model1)$adj.r.squared)

# test_trink is the function which tests all possible (those which are not included yet) nucleoties and neibouring nucleotide combinations taking fomula, source data.frame and parameter for crossvalidation as arguments.
res <- test_trink(model_formula, retention, cv_m = m_for_cv)

while(dim(res)[1] > 0){
  addedNk <- levels(res[1,]$component)[which(levels(res[1,]$component) == res[1,]$component)]
  model_formula <- paste(model_formula,addedNk, sep = ' + ')
  model1 <- lm(model_formula, data=retention)
  MSE_tab <- rbind(MSE_tab, data.frame(nk = addedNk, MSE = res$MSE[1], cvMSE = res$cvMSE[1], Adj.R2 = summary(model1)$adj.r.squared))
  res <- test_trink(model_formula, retention, cv_m = m_for_cv)
}

write.csv(MSE_tab, file = './output/3_deviant_nucleotides.csv')

capture.output(summary(model1), file = './output/3_NGA_model_summary.txt')
capture.output(print('ANOVA'), file = './output/3_NGA_model_summary.txt', append = T)
capture.output(anova(model1), file = './output/3_NGA_model_summary.txt', append = T)