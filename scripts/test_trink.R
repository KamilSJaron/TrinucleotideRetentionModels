test_trink <- function(model_formula, retention, table_out = F, alpha = 0.05, cv_m = 27){
  model1 <- lm(model_formula, data=retention)
  stnks <- c("A_1","G_1","T_1")
  middlenks = c("A_2","G_2","T_2")
  lastnks <- c("A_3","G_3","T_3")
  total_tests <- length(stnks) * length(middlenks) + 
    length(lastnks) * length(middlenks) + 
    length(stnks) + length(middlenks) + length(lastnks)
  trink_buffer <- data.frame(component = numeric(0), p_value = numeric(0), MSE = numeric(0), cvMSE = numeric(0))
  
  for(nk2 in middlenks){ 
    for(nk in stnks){
      f <- model_formula
      if(grepl(paste(nk,nk2,sep = ':'),f)){ next }
      f <- paste(f, paste(nk,nk2, sep=":") , sep=" + ")
      model2 <- lm(f, data=retention)
      if(anova(model1,model2)$`Pr(>F)`[2] < alpha / total_tests){
        kfCV <- suppressMessages(suppressWarnings(cv.lm(retention, model2, m = cv_m, printit = F)))
        trink_buffer <- rbind(trink_buffer, data.frame(component = paste(nk,nk2, sep=":"), 
                                                       p_value = anova(model1,model2)$`Pr(>F)`[2], 
                                                       MSE = mean(model2$residuals^2),
                                                       cvMSE = mean((kfCV$cvpred - kfCV$logk)^2)
                                                       )
                              )
      } 
    }
    for(nk3 in lastnks){
      f <- model_formula
      if(grepl(paste(nk2,nk3,sep = ':'),f)){ next }
      f <- paste(f, paste(nk2,nk3, sep=":") , sep=" + ")
      model2 <- lm(f, data=retention)
      if(anova(model1,model2)$`Pr(>F)`[2] < alpha / total_tests){
        kfCV <- suppressMessages(suppressWarnings(cv.lm(retention, model2, m = cv_m, printit = F)))
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
    f <- model_formula
    if(grepl(paste(' ',nk,sep = ''),f)){ next }
    f <- paste(f, nk, sep=" + ")
    model2 <- lm(f, data=retention)
    if(anova(model1,model2)$`Pr(>F)`[2] < alpha / total_tests){
      kfCV <- suppressWarnings(cv.lm(retention, model2, m = cv_m, printit = F))
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
    return(trink_buffer) 
  } else {
    print(trink_buffer)
    return(trink_buffer[which.min(trink_buffer$p_value),])
  }
}