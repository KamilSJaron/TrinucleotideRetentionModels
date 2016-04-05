# To support following statement in the article we performed a permutation test of AS-model
# statement:
# Despite relatively low number of measurements (27) and high number of predictive variables, both models relate a position of the base to the overall retention of the trinucleotide with significance far beyond just a coincidence, confirming thus a sequence dependent retention, although it gives no further hints on the retention mechanism.

BS_model <- lm(formula = logk ~ st + nd + rd, data = retention)
AS_model <- lm(logk ~ st + nd + rd + st:nd + nd:rd, data = retention)

# computes our P value and r squared of BS and AS model
BS_pval <- anova(BS_model)$`Pr(>F)`[2]
AS_pval <- anova(AS_model)$`Pr(>F)`[2]

BS_r <- summary(BS_model)$adj.r.squared
AS_r <- summary(AS_model)$adj.r.squared

r_sqs <- c()
r_sqs_int <- c() 
p_vals <- c()
p_vals_int <- c()

# ensures reproductivility
set.seed(160215)

# saving log of retention factor out of non-randomized data frame
logk <- retention$logk

for(i in 1:10000){
  
  # ranodomly permutes the nucleotides on first second and trird position
  st <- sample(retention$st)
  nd <- sample(retention$nd)
  rd <- sample(retention$rd)
  
  # construction of randomized model
  r_mod <- lm(formula = logk ~ st + nd + rd)
  r_mod_i <- lm(formula = logk ~ st + nd + rd + st:nd + nd:rd)
   
  # capture of stats
  r_sqs <- c(r_sqs, summary(r_mod)$adj.r.squared)
  r_sqs_int <- c(r_sqs_int, summary(r_mod_i)$adj.r.squared) 
  p_vals <- c(p_vals, anova(r_mod)$`Pr(>F)`[2])
  p_vals_int <- c(p_vals_int, anova(r_mod_i)$`Pr(>F)`[2])
} 

pdf('./output/s1_rsq_hist_BSmodel.pdf')
hist(r_sqs, main = 'random models logk ~ st + nd + rd', xlab = 'adj. R srqared')
dev.off()

pdf('./output/s1_rsq_hist_ASmodel.pdf')
hist(r_sqs_int, main = 'logk ~ st + nd + rd + st:nd + nd:rd', xlab = 'adj. R srqared')
dev.off()

pdf('./output/s1_pval_hist_BSmodel.pdf')
  hist(p_vals, main = 'random models logk ~ st + nd + rd', xlab = 'p-val')
dev.off()

pdf('./output/s1_hist_pval_ASmodel.pdf')
hist(p_vals_int, main = 'logk ~ st + nd + rd + st:nd + nd:rd', xlab = 'p-val')
dev.off()
