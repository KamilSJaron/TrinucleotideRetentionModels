# Flyvbjerg & Petersen 1989, ESTIMATORS BASED ON THE CORRELATION FUNCTIONS

# load sasa data
sasaTS = data.frame(order = 1:40000)
for(st in c('A','T','G')){
  for(nd in c('A','T','G')){
    for(rd in c('A','T','G')){
      trink <- paste(st,nd,rd,sep = '')
      file <- paste('./data/sasa/',trink,'.sasa', sep = '')
      trink_values <- read.csv(file,header = F, sep = ' ', dec = '.')[,3]
      sasaTS[,as.character(trink)] <- trink_values
    }
  }
}

# smaller computing functions
setCutoff <- function(tau = 15, threshold = 0.05, Tmin = 0, Tmax = 100){
  ser <- exp(-seq(Tmin,Tmax) / tau)
  return(seq(Tmin,Tmax)[min(which(ser < threshold))])
}

getCt <- function(x,rs){
  n <- length(x)
  ct <- c()
  for(r in rs){
    ct <- c(ct,mean((x[1:(n-r)] - mean(x[1:(n-r)])) * (x[(1+r):n] - mean(x[(1+r):n]))))
  }
  return(ct)
}

getAutocorVar <- function(x, cutoff){
  n <- length(x)
  #upper and lower sides of fraction
  upper <- getCt(x, 0) + 2 * sum(1 - ((1:cutoff) / n) * getCt(x,1:cutoff))
  lower <- (n - 2 * cutoff - 1 + (cutoff * (cutoff + 1) / n))
  return((upper / lower))
}

#estimation of variance for all trinucleotides
vartab <- data.frame(trink = numeric(0), mean = numeric(0), tau = numeric(0), T = numeric(0), var = numeric(0))
for(st in c('A','T','G')){
  for(nd in c('A','T','G')){
    for(rd in c('A','T','G')){
      trink = paste(st,nd,rd,sep = '')
      x <- sasaTS[,trink]
      trink_acf <- acf(x, lag.max = NULL, type = c("correlation"), plot = F, na.action = na.fail, demean = TRUE)
      tau <- min(which(trink_acf$acf < 0.1))
      cutoff <- setCutoff(tau,0.05,1,200)
      vartab <- rbind(vartab, data.frame(trink = trink, mean = mean(x), tau = tau, T = cutoff, var = getAutocorVar(x, cutoff)))
    }
  }
}

vartab$sd <- sqrt(vartab$var)

write.csv(vartab, './output/s2_sasa_sd.csv')