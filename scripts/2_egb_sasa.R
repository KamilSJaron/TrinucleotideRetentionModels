# ======= CV physically correct model =======
# log(k) = a*E_gb + b*SASA + c

NNN_model <- lm(logk ~ egb + sasaALL, data = retention)
kfCV <- cv.lm(retention, NNN_model, m = 27, plotit = F)

capture.output(summary(NNN_model), file = './output/2_NNN_model_summary.txt')
capture.output(print('ANOVA'), file = './output/2_NNN_model_summary.txt', append = T)
capture.output(anova(NNN_model), file = './output/2_NNN_model_summary.txt', append = T)

pdf('./output/2_NNN_model_residuals.pdf')
  hist(NNN_model$residuals, main='residuals of the NNN model')
dev.off()

capture.output(print('MSE'), file = './output/2_NNN_model_summary.txt', append = T)
capture.output(mean(NNN_model$residuals^2), file = './output/2_NNN_model_summary.txt', append = T)
capture.output(print('MSE crossvalidated'), file = './output/2_NNN_model_summary.txt', append = T)
capture.output(mean((kfCV$cvpred - kfCV$logk)^2), file = './output/2_NNN_model_summary.txt', append = T)




