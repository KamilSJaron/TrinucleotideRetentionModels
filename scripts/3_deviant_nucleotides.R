# ======= CV fyzikalne spravneho modelu =======
# log(k) = a*E_gb + b*SASA + c

retention$sasaALL <- retention$sasaPSAsur + retention$sasaNSAsur

fcor_model <- lm(logk ~ egb + sasaALL, data = retention)
fcor_model_wo_intercept <- lm(logk ~ -1 + egb + sasaALL, data = retention)

summary(fcor_model)
summary(fcor_model_wo_intercept)

kfCV <- cv.lm(retention, fcor_model, m = 27, plotit = F)
kfCV_wo_intercept <- cv.lm(retention, fcor_model, m = 27, plotit = F)

anova(fcor_model)
anova(fcor_model_wo_intercept)

hist(fcor_model$residuals)
hist(fcor_model_wo_intercept$residuals)

mean(fcor_model$residuals^2)
mean((kfCV$cvpred - kfCV$logk)^2)

mean(fcor_model_wo_intercept$residuals^2)
mean((kfCV_wo_intercept$cvpred - kfCV_wo_intercept$logk)^2)