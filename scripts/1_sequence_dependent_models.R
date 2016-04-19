# ------------SEQ MODELS----------------
# BS_model
BS_model <- lm(formula = logk ~ st + nd + rd, data = retention)

# cv.lm computes a crossvalidation. m is equal to number of mesurements, therefore it is leave-one crossvalidation
kfCV <- cv.lm(retention, BS_model, m = 27, plotit = F)

BS_summary_file <- './output/1_BS_model_summary.txt'

# BS model statistics
capture.output(summary(BS_model), file = BS_summary_file)
capture.output(print('ANOVA'), file = BS_summary_file, append = T)
capture.output(anova(BS_model), file = BS_summary_file, append = T)
capture.output(print('Shapiro-wilk test of residual errors'), append = T, file = BS_summary_file)
capture.output(shapiro.test(BS_model$residuals), file = BS_summary_file, append = T)
capture.output(print('MSE'), append = T, file = BS_summary_file)
capture.output(mean(BS_model$residuals^2), file = BS_summary_file, append = T)
capture.output(print('MSE crossvalidated'), append = T, file = BS_summary_file)
capture.output(mean((kfCV$cvpred - kfCV$logk)^2), file = BS_summary_file, append = T)

# AS_model
AS_model <- lm(logk ~ st + nd + rd + st:nd + nd:rd, data = retention)
kfCV2 <- cv.lm(retention, AS_model, m = 27, plotit = F)
AS_summary_file <- './output/1_AS_model_summary.txt'

# AS model statistics
capture.output(summary(AS_model), file = AS_summary_file)
capture.output(print('ANOVA'), file = AS_summary_file, append = T)
capture.output(anova(AS_model), file = AS_summary_file, append = T)
capture.output(print('Shapiro-wilk test of residual errors'), append = T, file = AS_summary_file)
capture.output(shapiro.test(AS_model$residuals), file = AS_summary_file, append = T)
capture.output(print('MSE'), append = T, file = AS_summary_file)
capture.output(mean(AS_model$residuals^2), file = AS_summary_file, append = T)
capture.output(print('MSE crossvalidated'), append = T, file = AS_summary_file)
capture.output(mean((kfCV2$cvpred - kfCV2$logk)^2), file = AS_summary_file, append = T)

# AS BS submodel test
capture.output(anova(BS_model,AS_model), file = './output/1_BS_AS_model_ANOVA.txt')
# at this point, we see a strong dependency of retention factor to sequence of trinucletide.