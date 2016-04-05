# ------------SEQ MODELS----------------
# position and nk
BS_model <- lm(formula = logk ~ st + nd + rd, data = retention)
kfCV <- cv.lm(retention, BS_model, m = 27, plotit = F)

capture.output(summary(BS_model), file = './output/1_BS_model_summary.txt')
capture.output(print('ANOVA'), file = './output/1_BS_model_summary.txt', append = T)
capture.output(anova(BS_model), file = './output/1_BS_model_summary.txt', append = T)
capture.output(print('Shapiro-wilk test of residual errors'), append = T, file = './output/1_BS_model_summary.txt')
capture.output(shapiro.test(BS_model$residuals), file = './output/1_BS_model_summary.txt', append = T)
capture.output(print('MSE'), append = T, file = './output/1_BS_model_summary.txt')
capture.output(mean(BS_model$residuals^2), file = './output/1_BS_model_summary.txt', append = T)
capture.output(print('MSE crossvalidated'), append = T, file = './output/1_BS_model_summary.txt')
capture.output(mean((kfCV$cvpred - kfCV$logk)^2), file = './output/1_BS_model_summary.txt', append = T)

# position, nk and neighbours
AS_model <- lm(logk ~ st + nd + rd + st:nd + nd:rd, data = retention)
kfCV2 <- cv.lm(retention, AS_model, m = 27, plotit = F)

capture.output(summary(AS_model), file = './output/1_AS_model_summary.txt')
capture.output(print('ANOVA'), file = './output/1_AS_model_summary.txt', append = T)
capture.output(anova(AS_model), file = './output/1_AS_model_summary.txt', append = T)
capture.output(print('Shapiro-wilk test of residual errors'), append = T, file = './output/1_AS_model_summary.txt')
capture.output(shapiro.test(AS_model$residuals), file = './output/1_AS_model_summary.txt', append = T)
capture.output(print('MSE'), append = T, file = './output/1_AS_model_summary.txt')
capture.output(mean(AS_model$residuals^2), file = './output/1_AS_model_summary.txt', append = T)
capture.output(print('MSE crossvalidated'), append = T, file = './output/1_AS_model_summary.txt')
capture.output(mean((kfCV2$cvpred - kfCV2$logk)^2), file = './output/1_AS_model_summary.txt', append = T)

capture.output(anova(BS_model,AS_model), file = './output/1_BS_AS_model_ANOVA.txt')
# at this point, we see a strong dependency of retention factor to sequence of trinucletide.