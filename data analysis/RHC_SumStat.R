### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~ PS Methods RHC Data Application ~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### Summary statistics and figures of the causal effects

### by Yi Liu

################ Weighting for Paper ################

# source("RHC_Causal_Effects.R")
load("RHC_Results.RData")
library(PSweight)
library(xtable)
library(tidyverse)

## PS and outcome models
ps.mult <- swang1 ~ cardiohx + chfhx + dementhx + psychhx + chrpulhx + renalhx + liverhx + 
  gibledhx + malighx + immunhx + transhx + amihx + age + sex + edu + surv2md1 + das2d3pc + 
  aps1 + scoma1 + meanbp1 + wblc1 + hrt1 + resp1 + temp1 + pafi1 + alb1 + hema1 + bili1 + 
  crea1 + sod1 + pot1 + paco21 + ph1 + wtkilo1 + dnr1 + resp + card + neuro + gastr +
  renal + meta + hema + seps + trauma + ortho + income1 + income2 + income3 + raceblack + 
  raceother + ins_care + ins_pcare + ins_carecaid + ins_caid + ins_no + cat1_copd + 
  cat1_mosfsep + cat1_mosfmal + cat1_chf + cat1_coma + cat1_cirr + cat1_lung + cat1_colon +
  cat2_mosfsep + cat2_mosfmal + cat2_coma + cat2_cirr + cat2_lung + cat2_colon + ca_yes +
  ca_meta + wt0

### --- For ATC using PSweight package
atc.mult <- swang2 ~ cardiohx + chfhx + dementhx + psychhx + chrpulhx + renalhx + liverhx + 
  gibledhx + malighx + immunhx + transhx + amihx + age + sex + edu + surv2md1 + das2d3pc + 
  aps1 + scoma1 + meanbp1 + wblc1 + hrt1 + resp1 + temp1 + pafi1 + alb1 + hema1 + bili1 + 
  crea1 + sod1 + pot1 + paco21 + ph1 + wtkilo1 + dnr1 + resp + card + neuro + gastr +
  renal + meta + hema + seps + trauma + ortho + income1 + income2 + income3 + raceblack + 
  raceother + ins_care + ins_pcare + ins_carecaid + ins_caid + ins_no + cat1_copd + 
  cat1_mosfsep + cat1_mosfmal + cat1_chf + cat1_coma + cat1_cirr + cat1_lung + cat1_colon +
  cat2_mosfsep + cat2_mosfmal + cat2_coma + cat2_cirr + cat2_lung + cat2_colon + ca_yes +
  ca_meta + wt0


## Propensity Score Plot
bal.mult <- SumStat(ps.formula = ps.mult, weight = c("IPW", "overlap", "matching", "entropy"), 
                    data = RHC, delta = 0)
png("ps_rhc.png", res=72*2, width = 1200, height = 800)
plot(bal.mult, type = "hist")
dev.off()

## Love Plot
rhc_trim1 <- PStrim(data = RHC, ps.formula = ps.mult, delta = 0.05)
rhc_trim2 <- PStrim(data = RHC, ps.formula = ps.mult, delta = 0.1)
rhc_trim3 <- PStrim(data = RHC, ps.formula = ps.mult, delta = 0.15)

covar <- names(SumStat(ps.formula = ps.mult, weight = "IPW", data = RHC)$IPW.sumstat[, "ASD weighted var 1-2"])

bal.IPW <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = RHC)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE", covar = covar)
colnames(bal.IPW) <- c("ASD", "Method", "covar")
bal.IPW.5 <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = rhc_trim1$data)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE (0.05)", covar = covar)
colnames(bal.IPW.5) <- c("ASD", "Method", "covar")
bal.IPW.10 <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = rhc_trim2$data)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE (0.1)", covar = covar)
colnames(bal.IPW.10) <- c("ASD", "Method", "covar")
bal.IPW.15 <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = rhc_trim3$data)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE (0.15)", covar = covar)
colnames(bal.IPW.15) <- c("ASD", "Method", "covar")
bal.ATO <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "overlap", data = RHC)$overlap.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATO", covar = covar)
colnames(bal.ATO) <- c("ASD", "Method", "covar")
bal.ATM <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "matching", data = RHC)$matching.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATM", covar = covar)
colnames(bal.ATM) <- c("ASD", "Method", "covar")
bal.ATEN <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "entropy", data = RHC)$entropy.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATEN", covar = covar)
colnames(bal.ATEN) <- c("ASD", "Method", "covar")
bal.ATT <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "treated", data = RHC)$treated.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATT", covar = covar)
colnames(bal.ATT) <- c("ASD", "Method", "covar")
bal.ATC <- as.data.frame(SumStat(ps.formula = atc.mult, weight = "treated", data = RHC)$treated.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATC", covar = covar)
colnames(bal.ATC) <- c("ASD", "Method", "covar")

df <- rbind(bal.IPW, bal.IPW.5, bal.IPW.10, bal.IPW.15, bal.ATO, bal.ATM, bal.ATEN, bal.ATT, bal.ATC)
df$Method <- factor(df$Method, levels = c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", "ATO", "ATM", "ATEN", "ATT", "ATC"))

png("asd_rhc.png", res = 200, width = 1050, height = 1350)
ggplot(df, aes(x = ASD, y = covar, col = Method)) + geom_point(alpha = 0.85) + 
  labs(x = "Standarized Mean Difference", y = "Covariates") + 
  scale_color_manual(values = c("royalblue", "deepskyblue", "darkseagreen1", 
                                "forestgreen", "sienna", "indianred2", "plum2", 
                                "goldenrod2", "dimgrey")) + 
  geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed", color = "black", size = 0.2) + 
  theme_bw()
dev.off()

## Summary Statistics
p.RHC
dim(RHC)

n1/n
n2/n
n3/n

cols <- c("Sample Size", "Weighted Est", "Std.Error", "p-value", 
          "DR Est", "Std.Error", "p-value")
rows <- c("IPW", "IPW (0.05)", "IPW (0.1)", "IPW (0.15)",
          "ATO", "ATM", "ATEN", "ATC", "ATT")

p.val <- function(p) ifelse(p>=0.001, round(p, 3), "<0.001")
f.table <- data.frame(SS = c(n,n1,n2,n3,n,n,n,n,n),
                      
                      PE = c(IPW.wt[1], IPW_wt_trim[1], IPW_wt_trim1[1], IPW_wt_trim2[1],
                             ATO.wt[1], ATM.wt[1], ATEN.wt[1], -ATC.wt[1], ATT.wt[1]),
                      SD = c(IPW.wt[2], IPW_wt_trim[2], IPW_wt_trim1[2], IPW_wt_trim2[2],
                              ATO.wt[2], ATM.wt[2], ATEN.wt[2], ATC.wt[2], ATT.wt[2]),
                      p = p.val(c(IPW.wt[6], IPW_wt_trim[6], IPW_wt_trim1[6], IPW_wt_trim2[6],
                            ATO.wt[6], ATM.wt[6], ATEN.wt[6], ATC.wt[6], ATT.wt[6])),
                      
                      DR.Est = c(IPW.aug[1], IPW_aug_trim[1], IPW_aug_trim1[1], IPW_aug_trim2[1],
                             ATO.aug[1], ATM.aug[1], ATEN.aug[1], -ATC.aug[1], ATT.aug[1]),
                      DR.SD = c(IPW.aug[2], IPW_aug_trim[2], IPW_aug_trim1[2], IPW_aug_trim2[2],
                                 ATO.aug[2], ATM.aug[2], ATEN.aug[2], ATC.aug[2], ATT.aug[2]),
                      DR.p = p.val(c(IPW.aug[6], IPW_aug_trim[6], IPW_aug_trim1[6], IPW_aug_trim2[6],
                                     ATO.aug[6], ATM.aug[6], ATEN.wt[6], ATC.aug[6], ATT.aug[6])))

rownames(f.table) <- rows
colnames(f.table) <- cols

xtable(f.table, digits = 3)
 


################ Var DR ATT/ATC Paper ################
load("RHC_VarDR.RData")


### Covariates Balancing
covar <- names(SumStat(ps.formula = ps.mult, weight = "IPW", data = RHC)$IPW.sumstat[, "ASD weighted var 1-2"])

bal.ATT <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "treated", data = RHC)$treated.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATT", covar = covar)
colnames(bal.ATT) <- c("ASD", "Estimand", "covar")
bal.ATC <- as.data.frame(SumStat(ps.formula = atc.mult, weight = "treated", data = RHC)$treated.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATC", covar = covar)
colnames(bal.ATC) <- c("ASD", "Estimand", "covar")

df <- rbind(bal.ATT, bal.ATC)

png("asd_rhc_VarDR.png", res = 200, width = 1050, height = 1350)
ggplot(df, aes(x = ASD, y = covar, col = Estimand)) + geom_point(alpha = 0.85, size=2) + 
  labs(x = "Standarized Mean Difference", y = "Covariates") + 
  scale_color_manual(values = c("royalblue", "goldenrod2")) + 
  geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed", color = "black", size = 0.2) + 
  theme_bw()
dev.off()


### Summary Statistics of Variance Estimations
cols <- c("Effect", "Estimation", "Variance Method", "Standard Error", "p-value")
p.val <- function(p) ifelse(p>=0.001, round(p, 3), "<0.001")
f.table <- data.frame(Effect = rep(c("ATC", "ATT"), 3),
                      Method = c(rep("Asymptotic (Sandwich)", 2), rep("Wild bootstrap (Rademacher)", 2), rep("Wild bootstrap (Exponential)", 2)),
                      Est = c(ATC.sand$Est, ATT.sand$Est, ATC.wb1$Est, ATT.wb1$Est, ATC.wb2$Est, ATT.wb2$Est),
                      SD = c(ATC.sand$SD, ATT.sand$SD, ATC.wb1$SD, ATT.wb1$SD, ATC.wb2$SD, ATT.wb2$SD),
                      p = p.val(c(ATC.sand$p, ATT.sand$p, ATC.wb1$p, ATT.wb1$p, ATC.wb2$p, ATT.wb2$p)))

colnames(f.table) <- cols

print(xtable(f.table, digits = 3), include.rownames=FALSE)
