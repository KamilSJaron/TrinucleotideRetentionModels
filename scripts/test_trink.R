test_trink <- function(model1, retention, table_out = F, alpha = 0.05, cv_m = 27){
  stnks <- c("A_1","G_1","T_1")
  middlenks = c("A_2","G_2","T_2")
  lastnks <- c("A_3","G_3","T_3")
  total_tests <- length(stnks) * length(middlenks) + 
    length(lastnks) * length(middlenks) + 
    length(stnks) + length(middlenks) + length(lastnks)
  trink_buffer <- data.frame(component = 'Nothing significant', p_value = 1, MSE = Inf, cvMSE = Inf)
  
  for(nk2 in middlenks){ 
    for(nk in stnks){
      f <- as.character(model1$call)[2]
      if(grepl(paste(nk,nk2,sep = ':'),f)){ next }
      f <- paste(f, paste(nk,nk2, sep=":") , sep=" + ")
      model2 <- lm(f, data=retention)
      if(anova(model1,model2)$`Pr(>F)`[2] < alpha / total_tests){
        kfCV <- cv.lm(retention, model2, m = cv_m, plotit = F)
        trink_buffer <- rbind(trink_buffer, data.frame(component = paste(nk,nk2, sep=":"), 
                                                       p_value = anova(model1,model2)$`Pr(>F)`[2], 
                                                       MSE = mean(model2$residuals^2),
                                                       cvMSE = mean((kfCV$cvpred - kfCV$logk)^2)
                                                       )
                              )
      } 
    }
    for(nk3 in lastnks){
      f <- as.character(model1$call)[2]
      if(grepl(paste(nk2,nk3,sep = ':'),f)){ next }
      f <- paste(f, paste(nk2,nk3, sep=":") , sep=" + ")
      model2 <- lm(f, data=retention)
      if(anova(model1,model2)$`Pr(>F)`[2] < alpha / total_tests){
        kfCV <- cv.lm(retention, model2, m = cv_m, plotit = F)
        trink_buffer <- rbind(trink_buffer, data.frame(component = paste(nk2,nk3, sep=":"), 
                                                       p_value = anova(model1,model2)$`Pr(>F)`[2], 
                                                       MSE = mean(model2$residuals^2),
                                                       cvMSE = mean((kfCV$cvpred - kfCV$logk)^2)
                                                       )
                              )
      } 
    }
  }
  
  for(nk in c(stnks, middlenks, lastnks)){
    f <- as.character(model1$call)[2]
    if(grepl(paste(' ',nk,sep = ''),f)){ next }
    f <- paste(f, nk, sep=" + ")
    model2 <- lm(f, data=retention)
    if(anova(model1,model2)$`Pr(>F)`[2] < alpha / total_tests){
      kfCV <- cv.lm(retention, model2, m = cv_m, plotit = F)
      trink_buffer <- rbind(trink_buffer, 
                            data.frame(component = nk, 
                                       p_value = anova(model1,model2)$`Pr(>F)`[2],
                                       MSE = mean(model2$residuals^2),
                                       cvMSE = mean((kfCV$cvpred - kfCV$logk)^2)
                                       )
                            )
    }
  }
  
  if(table_out){
    return(trink_buffer[-1,]) 
  } else {
    print(trink_buffer[-1,])
    return(trink_buffer[which.min(trink_buffer$p_value),])
  }

}