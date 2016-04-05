setwd('160113_finalization/')
source('./160113_load_data.R')

# ------------LOG MODELS----------------
# position and nk
m_pos_nk <- lm(formula = logk ~ st + nd + rd, data = retention)
kfCV <- cv.lm(retention, m_pos_nk, m = 27, plotit = F)

summary(m_pos_nk)
mean(m_pos_nk$residuals^2)
mean((kfCV$cvpred - kfCV$logk)^2)

# position, nk and neighbours
m_pos_nk_int <- lm(logk ~ st + nd + rd + st:nd + nd:rd, data = retention)
kfCV2 <- cv.lm(retention, m_pos_nk_int, m = 27, plotit = F)

summary(m_pos_nk_int)
anova(m_pos_nk_int)
shapiro.test(m_pos_nk_int$residuals)
mean(m_pos_nk_int$residuals^2)
mean((kfCV2$cvpred - kfCV2$logk)^2)

anova(m_pos_nk,m_pos_nk_int)
# at this point, we see a strong dependency of retention factor to sequence of trinucletide.

# # position independent interactions, position of nk
# retention$AT <- 
#   (retention$st == 'A' & retention$nd == 'T') +
#   (retention$nd == 'A' & retention$rd == 'T') +
#   (retention$st == 'T' & retention$nd == 'A') +
#   (retention$nd == 'T' & retention$rd == 'A')
# retention$AG <- 
#   (retention$st == 'A' & retention$nd == 'G') +
#   (retention$nd == 'A' & retention$rd == 'G') +
#   (retention$st == 'G' & retention$nd == 'A') +
#   (retention$nd == 'G' & retention$rd == 'A')
# retention$TG <- 
#   (retention$st == 'G' & retention$nd == 'T') +
#   (retention$nd == 'G' & retention$rd == 'T') +
#   (retention$st == 'T' & retention$nd == 'G') +
#   (retention$nd == 'T' & retention$rd == 'G')
# retention$AA <- 
#   (retention$st == 'A' & retention$nd == 'A') +
#   (retention$nd == 'A' & retention$rd == 'A')
# retention$TT <- 
#   (retention$st == 'T' & retention$nd == 'T') +
#   (retention$nd == 'T' & retention$rd == 'T')
# retention$GG <- 
#   (retention$st == 'G' & retention$nd == 'G') +
#   (retention$nd == 'G' & retention$rd == 'G')
# 
# m_pos_nk_int <- lm(logk ~ st + nd + rd + TT + TG + GG, data = retention)
# kfCV2 <- cv.lm(retention, m_pos_nk_int, m = 27, plotit = F)
# 
# summary(m_pos_nk_int)
# anova(m_pos_nk_int)
# shapiro.test(m_pos_nk_int$residuals)
# mean(m_pos_nk_int$residuals^2)
# mean((kfCV2$cvpred - kfCV2$logk)^2)
