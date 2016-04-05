# solve how to read is (one way out is to replace every pair of spaces by a space)

sasaTS = data.frame(order = 1:40000)
for(st in c('A','T','G')){
  for(nd in c('A','T','G')){
    for(rd in c('A','T','G')){
      trink <- paste(st,nd,rd,sep = '')
      file <- paste('/Volumes/dump/data/chemical_properties_nucl/sasa/',trink,'.sasa', sep = '')
      trink_values <- read.csv(file,header = F, sep = ' ', dec = '.')[,3]
      sasaTS[,as.character(trink)] <- trink_values
    }
  }
}

CIs <- data.frame(trink = '', lb = 0, mean = 0, rb = 0)
for(st in c('A','T','G')){
  for(nd in c('A','T','G')){
    for(rd in c('A','T','G')){
      trink = paste(st,nd,rd,sep = '')
      n = length(sasaTS[,trink]) 
      s = sd(sasaTS[,trink])        # sample standard deviation 
      SE = s/sqrt(n)             # standard error estimate 
      E = qt(.975, df=n-1)*SE     # margin of error 
      xbar = mean(sasaTS[,trink])   # sample mean 
      CIs <- rbind(CIs, data.frame(trink = trink, lb = xbar-E, mean = xbar, rb = xbar+E))
    }
  }
}
CIs <- CIs[-1,]

source('/Volumes/dump/scripts/R/start_up/table_print.R')
print.wikitable(CIs)

for(st in c('A','T','G')){
  for(nd in c('A','T','G')){
    for(rd in c('A','T','G')){
      trink <- paste(st,nd,rd,sep = '')
      acf(sasaTS[,trink], lag.max = NULL,
          type = c("correlation"),
          plot = TRUE, na.action = na.fail, demean = TRUE)
      readline();
    }
  }
}


AAAacf <- acf(sasaTS[,'AAA'], lag.max = NULL,
    type = c("correlation"),
    plot = TRUE, na.action = na.fail, demean = TRUE)


tau <- min(which(AAAacf$acf < 0.05))
x <- sasaTS[,'AAA']

setCutoff <- function(tau = 15, threshold = 0.05, Tmin = 0, Tmax = 100){
  ser <- exp(-seq(Tmin,Tmax) / tau)
  return(seq(Tmin,Tmax)[min(which(ser < threshold))])
}

cutoff <- setCutoff(tau)

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

vartab <- data.frame(trink = '', mean = 0, tau = 0, T = 0, var = 0)
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
vartab <- vartab[-1,]

vartab$sd <- sqrt(vartab$var)
vartab


source('/Volumes/dump/scripts/R/start_up/table_print.R')
print.wikitable(vartab)

write.csv(vartab,'~/Downloads/160311SASA_errors.csv')




