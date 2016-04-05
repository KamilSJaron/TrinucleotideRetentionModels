# rm(list = ls())
options(digits=5)
setwd('/Volumes/dump/projects/Retence/160113_finalization/')
source('./160115_load_data.R')
source('./test_trink.R')

# ======= CV fyzikalne spravneho modelu =======
# log(k) = a*E_gb + b*SASA + c

model1 <- lm(logk ~ egb + sasaALL, data=retention)
model1_wi <- lm(logk ~ -1 + egb + sasaALL, data=retention)
kfCV <- cv.lm(retention, model1, m = 27, plotit = F)
kfCV_wi <- cv.lm(retention, model1_wi, m = 27, plotit = F)
anova(model1, model1_wi)
m_for_cv = 27

MSE_tab <- data.frame(nk = 'intercept', MSE = mean(model1$residuals^2), cvMSE = mean((kfCV$cvpred - kfCV$logk)^2), Adj.R2 = summary(model1)$adj.r.squared,
                      nk_wi = '-', MSE_wi = mean(model1_wi$residuals^2), MSEcv_vi = mean((kfCV_wi$cvpred - kfCV_wi$logk)^2), Adj.R2wi = summary(model1_wi)$adj.r.squared)

res <- test_trink(model1, retention, cv_m = m_for_cv)
res <- rbind(res,test_trink(model1_wi, retention, cv_m = m_for_cv))

addedNk <- levels(res[1,]$component)[which(levels(res[1,]$component) == res[1,]$component)]
addedNk_wi <- levels(res[2,]$component)[which(levels(res[2,]$component) == res[2,]$component)]

model1 <- lm(logk ~ egb + sasaALL + A_2, data=retention)
model1_wi <- lm(logk ~ -1 + egb + sasaALL + A_2, data=retention)
anova(model1, model1_wi)

MSE_tab <- rbind(MSE_tab, data.frame(nk = addedNk, MSE = res$MSE[1], cvMSE = res$cvMSE[1], Adj.R2 = summary(model1)$adj.r.squared,
                                     nk_wi = addedNk_wi, MSE_wi = res$MSE[2], MSEcv_vi = res$cvMSE[2], Adj.R2wi = summary(model1_wi)$adj.r.squared))

res <- rbind(test_trink(model1, retention, cv_m = m_for_cv),test_trink(model1_wi, retention, cv_m = m_for_cv))

addedNk <- levels(res[1,]$component)[which(levels(res[1,]$component) == res[1,]$component)]
addedNk_wi <- levels(res[2,]$component)[which(levels(res[2,]$component) == res[2,]$component)]

model1 <- lm(logk ~ egb + sasaALL + A_2 + T_1:G_2, data=retention)
model1_wi <- lm(logk ~ -1 + egb + sasaALL + A_2 + T_1:G_2, data=retention)
anova(model1, model1_wi)

MSE_tab <- rbind(MSE_tab, data.frame(nk = addedNk, MSE = res$MSE[1], cvMSE = res$cvMSE[1], Adj.R2 = summary(model1)$adj.r.squared,
                                     nk_wi = addedNk_wi, MSE_wi = res$MSE[2], MSEcv_vi = res$cvMSE[2], Adj.R2wi = summary(model1_wi)$adj.r.squared))

res <- rbind(test_trink(model1, retention, cv_m = m_for_cv),test_trink(model1_wi, retention, cv_m = m_for_cv))

addedNk <- levels(res[1,]$component)[which(levels(res[1,]$component) == res[1,]$component)]
addedNk_wi <- levels(res[2,]$component)[which(levels(res[2,]$component) == res[2,]$component)]

model1 <- lm(logk ~ egb + sasaALL + A_2 + T_1:G_2 + A_1, data=retention)
model1_wi <- lm(logk ~ -1 + egb + sasaALL + A_2 + T_1:G_2 + A_1, data=retention)
anova(model1, model1_wi)

MSE_tab <- rbind(MSE_tab, data.frame(nk = addedNk, MSE = res$MSE[1], cvMSE = res$cvMSE[1], Adj.R2 = summary(model1)$adj.r.squared,
                                     nk_wi = addedNk_wi, MSE_wi = res$MSE[2], MSEcv_vi = res$cvMSE[2], Adj.R2wi = summary(model1_wi)$adj.r.squared))

res <- rbind(test_trink(model1, retention, cv_m = m_for_cv),test_trink(model1_wi, retention, cv_m = m_for_cv))

addedNk <- levels(res[1,]$component)[which(levels(res[1,]$component) == res[1,]$component)]
addedNk_wi <- levels(res[2,]$component)[which(levels(res[2,]$component) == res[2,]$component)]

model1 <- lm(logk ~ egb + sasaALL + A_2 + T_1:G_2 + A_1 + T_2:A_3, data=retention)
model1_wi <- lm(logk ~ -1 + egb + sasaALL + A_2 + T_1:G_2 + A_1 + T_2:A_3, data=retention)
anova(model1, model1_wi)

MSE_tab <- rbind(MSE_tab, data.frame(nk = addedNk, MSE = res$MSE[1], cvMSE = res$cvMSE[1], Adj.R2 = summary(model1)$adj.r.squared,
                                     nk_wi = addedNk_wi, MSE_wi = res$MSE[2], MSEcv_vi = res$cvMSE[2], Adj.R2wi = summary(model1_wi)$adj.r.squared))

res <- rbind(test_trink(model1, retention, cv_m = m_for_cv),test_trink(model1_wi, retention, cv_m = m_for_cv))

addedNk <- levels(res[1,]$component)[which(levels(res[1,]$component) == res[1,]$component)]
addedNk_wi <- levels(res[2,]$component)[which(levels(res[2,]$component) == res[2,]$component)]

model1 <- lm(logk ~ egb + sasaALL + A_2 + T_1:G_2 + A_1 + T_2:A_3 + G_2:A_3, data=retention)
model1_wi <- lm(logk ~ -1 + egb + sasaALL + A_2 + T_1:G_2 + A_1 + T_2:A_3, data=retention)
anova(model1, model1_wi)

MSE_tab <- rbind(MSE_tab, data.frame(nk = addedNk, MSE = res$MSE[1], cvMSE = res$cvMSE[1], Adj.R2 = summary(model1)$adj.r.squared,
                                     nk_wi = addedNk_wi, MSE_wi = res$MSE[2], MSEcv_vi = res$cvMSE[2], Adj.R2wi = NA))

res <- rbind(test_trink(model1, retention, cv_m = m_for_cv),test_trink(model1_wi, retention, cv_m = m_for_cv))

addedNk <- levels(res[1,]$component)[which(levels(res[1,]$component) == res[1,]$component)]
addedNk_wi <- levels(res[2,]$component)[which(levels(res[2,]$component) == res[2,]$component)]

MSE_tab <- rbind(MSE_tab, data.frame(nk = addedNk, MSE = res$MSE[1], cvMSE = res$cvMSE[1], Adj.R2 = summary(model1)$adj.r.squared,
                                     nk_wi = addedNk_wi, MSE_wi = res$MSE[2], MSEcv_vi = res$cvMSE[2], Adj.R2wi = summary(model1_wi)$adj.r.squared))

source('/Volumes/dump/scripts/R/start_up/table_print.R')
MSE_tab <- MSE_tab[-7,]
MSE_tab

# print.wikitable(MSE_tab)
# 
# 
# pdf('EGB_SASA.pfd')
#   plot(MSE_tab$MSEcv_vi, ylim = c(0,MSE_tab$MSEcv_vi[1]), pch = 20, col = 'red', xaxt='n', ann=FALSE, ylab = "MSE")
#   points(MSE_tab$MSE_wi, col = 'red')
#   
#   points(MSE_tab$MSE)
#   points(MSE_tab$cvMSE, pch = 20)
#   
#   legend('topright', c('MSE', 'cvMSE','MSE_wi','cvMSE_wi'), pch = c(1, 20, 1, 20), col = c('black','black','red','red'), title = 'log(k) = a*E_gb + b*SASA + c', bty = 'n')
#   
#   axis(levels(MSE_tab$nk)[1:6], side = 1, at = 1:6)
# dev.off()


MSE_tab$nk_wi <- MSE_tab$cvMSE / MSE_tab$MSE
MSE_tab$cvMSEratio_wi <- MSE_tab$MSEcv_vi / MSE_tab$MSE_wi
colnames(MSE_tab)[5] <- 'cvMSEratio'

print.wikitable(MSE_tab)

model1 <- lm(logk ~ egb + sasaALL + A_2 + T_1:G_2 + A_1 + T_2:A_3 + G_2:A_3, data=retention)
plot(model1$residuals)
hist(model1$residuals)
shapiro.test(model1$residuals)
summary(model1)

model1 <- lm(logk ~ egb + sasaALL, data=retention)
plot(model1$residuals)
hist(model1$residuals)
shapiro.test(model1$residuals)
summary(model1)
