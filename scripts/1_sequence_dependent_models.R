# ------------SEQ MODELS----------------
# position and nk
BS_model <- lm(formula = logk ~ st + nd + rd, data = retention)
kfCV <- cv.lm(retention, BS_model, m = 27, plotit = F)

summary(BS_model)
mean(BS_model$residuals^2)
mean((kfCV$cvpred - kfCV$logk)^2)

# position, nk and neighbours
AS_model <- lm(logk ~ st + nd + rd + st:nd + nd:rd, data = retention)
kfCV2 <- cv.lm(retention, AS_model, m = 27, plotit = F)

summary(AS_model)
anova(AS_model)
shapiro.test(AS_model$residuals)
mean(AS_model$residuals^2)
mean((kfCV2$cvpred - kfCV2$logk)^2)

anova(BS_model,AS_model)
# at this point, we see a strong dependency of retention factor to sequence of trinucletide.