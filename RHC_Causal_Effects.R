### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~ PS Methods RHC Data Application ~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### Data cleaning and causal estimates

### by Yi Liu

RHC.org <- read.csv("rhc.csv")

# install.packages("PSweight")
library(PSweight)
library(dplyr)

RHC <- RHC.org %>% select(-X, -t3d30, -dth30, -adld3p, -urin1, -ptid)

## Data cleaning and coding dummy variables
RHC$sex <- NA
RHC$sex <- as.factor(ifelse(RHC.org$sex == "Female", 1, 0))

RHC <- RHC %>% mutate(income1 = as.factor(ifelse(income == "$11-$25k", 1, 0)),
                      income2 = as.factor(ifelse(income == "$25-$50k", 1, 0)),
                      income3 = as.factor(ifelse(income == "> $50k", 1, 0))) %>%
  select(-income)

RHC <- RHC %>% mutate(raceblack = as.factor(ifelse(race == "black", 1, 0)),
                      raceother = as.factor(ifelse(race == "other", 1, 0))) %>% 
  select(-race)

RHC <- RHC %>% mutate(ins_care = as.factor(ifelse(ninsclas == "Medicare", 1, 0)),
                      ins_pcare = as.factor(ifelse(ninsclas == "Private & Medicare", 1, 0)),
                      ins_carecaid = as.factor(ifelse(ninsclas == "Medicare & Medicaid", 1, 0)),
                      ins_caid = as.factor(ifelse(ninsclas == "Medicaid", 1, 0)),
                      ins_no = as.factor(ifelse(ninsclas == "No insurance", 1, 0))) %>% 
  select(-ninsclas)

RHC <- RHC %>% mutate(cat1_copd = as.factor(ifelse(cat1 == "COPD", 1, 0)),
                      cat1_mosfsep = as.factor(ifelse(cat1 == "MOSF w/Sepsis", 1, 0)),
                      cat1_mosfmal = as.factor(ifelse(cat1 == "MOSF w/Malignancy", 1, 0)),
                      cat1_chf = as.factor(ifelse(cat1 == "CHF", 1, 0)),
                      cat1_coma = as.factor(ifelse(cat1 == "Coma", 1, 0)),
                      cat1_cirr = as.factor(ifelse(cat1 == "Cirrhosis", 1, 0)),
                      cat1_lung = as.factor(ifelse(cat1 == "Lung Cancer", 1, 0)),
                      cat1_colon = as.factor(ifelse(cat1 == "Colon Cancer", 1, 0))) %>% 
  select(-cat1)

RHC$cat2[is.na(RHC$cat2)==TRUE] <- 1
RHC <- RHC %>% mutate(cat2_mosfsep = as.factor(ifelse(cat2 == "MOSF w/Sepsis", 1, 0)),
                      cat2_mosfmal = as.factor(ifelse(cat2 == "MOSF w/Malignancy", 1, 0)),
                      cat2_coma = as.factor(ifelse(cat2 == "Coma", 1, 0)),
                      cat2_cirr = as.factor(ifelse(cat2 == "Cirrhosis", 1, 0)),
                      cat2_lung = as.factor(ifelse(cat2 == "Lung Cancer", 1, 0)),
                      cat2_colon = as.factor(ifelse(cat2 == "Colon Cancer", 1, 0))) %>%
  select(-cat2)

RHC <- RHC %>% mutate(ca_yes = as.factor(ifelse(ca == "Yes", 1, 0)),
                      ca_meta = as.factor(ifelse(ca == "Metastatic", 1, 0))) %>% 
  select(-ca)

RHC[, c("resp", "gastr", "card", "renal", "meta", "hema", 
        "seps", "trauma", "ortho", "dnr1", "neuro")] <- NA
RHC$resp <- as.factor(ifelse(RHC.org$resp == "Yes", 1, 0))
RHC$gastr <- as.factor(ifelse(RHC.org$gastr == "Yes", 1, 0))
RHC$card <- as.factor(ifelse(RHC.org$card == "Yes", 1, 0))
RHC$renal <- as.factor(ifelse(RHC.org$renal == "Yes", 1, 0))
RHC$meta <- as.factor(ifelse(RHC.org$meta == "Yes", 1, 0))
RHC$hema <- as.factor(ifelse(RHC.org$hema == "Yes", 1, 0))
RHC$seps <- as.factor(ifelse(RHC.org$seps == "Yes", 1, 0))
RHC$trauma <- as.factor(ifelse(RHC.org$trauma == "Yes", 1, 0))
RHC$ortho <- as.factor(ifelse(RHC.org$ortho == "Yes", 1, 0))
RHC$dnr1 <- as.factor(ifelse(RHC.org$dnr1 == "Yes", 1, 0))
RHC$neuro <- as.factor(ifelse(RHC.org$neuro == "Yes", 1, 0))

RHC <- RHC %>% mutate(wt0 = as.factor(ifelse(wtkilo1 == 0, 1, 0)))

RHC$aps1 <- as.numeric(RHC$aps1)
RHC$scoma1 <- as.numeric(RHC$scoma1)
RHC$hrt1 <- as.numeric(RHC$hrt1)
RHC$sod1 <- as.numeric(RHC$sod1)

## Proportion of the treatment
table(RHC$swang1)
p.RHC <- length(RHC$swang1[RHC$swang1=="RHC"])/length(RHC$swang1)

## Outcome
RHC <- RHC %>% mutate(icusty = log(dschdte - sadmdte))
RHC$icusty[is.na(RHC$dschdte)] <- 
  log(RHC$dthdte[is.na(RHC$dschdte)] - RHC$sadmdte[is.na(RHC$dschdte)])
sum(is.na(RHC$icusty)) # Check missing value


### Descriptive summary statistics
V.name <- c("cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx", "renalhx", "liverhx", "gibledhx", 
            "malighx", "immunhx", "transhx", "amihx", "age", "sex", "edu", "surv2md1", 
            "das2d3pc", "aps1", "scoma1", "meanbp1", "wblc1", "hrt1", "resp1", "temp1", "pafi1", "alb1", 
            "hema1", "bili1", "crea1", "sod1", "pot1", "paco21", "ph1", "wtkilo1", "dnr1", "resp", "card", 
            "neuro", "gastr", "renal", "meta", "hema", "seps", "trauma", "ortho", "income1", "income2", 
            "income3", "raceblack", "raceother", "ins_care", "ins_pcare", "ins_carecaid", "ins_caid", 
            "ins_no", "cat1_copd", "cat1_mosfsep", "cat1_mosfmal", "cat1_chf", "cat1_coma", "cat1_cirr", 
            "cat1_lung", "cat1_colon", "cat2_mosfsep", "cat2_mosfmal", "cat2_coma", "cat2_cirr", 
            "cat2_lung", "cat2_colon", "ca_yes", "ca_meta", "wt0", "swang1", "icusty")
ctrl <- RHC[,V.name] %>% filter(swang1=="No RHC") %>% select(-swang1)
trt <- RHC[,V.name] %>% filter(swang1=="RHC") %>% select(-swang1)
ctrl <- apply(ctrl, 2, as.numeric)
trt <- apply(trt, 2, as.numeric)

mean.c <- round(apply(ctrl, 2, mean), 2)
sd.c <- round(apply(ctrl, 2, sd), 2)

mean.t <- round(apply(trt, 2, mean), 2)
sd.t <- round(apply(trt, 2, sd), 2)

descr <- c("Acute MI, peripheral vascular disease, severe cardiovascular symptoms (NYHA-Class III), very severe cardiovascular symptoms (NYHA-Class IV)", 
           "Congestive heart failure",
           "Dementia, stroke or cerebral infarct, Parkinson¡¯s disease",
           "Psychiatric history, active psychosis or severe depression",
           "Pulmonary disease",
           "Chronic renal disease, chronic hemodialysis or peritoneal dialysis",
           "Cirrhosis, hepatic failure",
           "Upper GI bleeding",
           "Solid tumor, metastatic disease, chronic leukemia/myeloma, acute leukemia, lymphoma",
           "Immunosupperssion, organ transplant, HIV positivity, diabetes mellitus without end organ damage, diabetes mellitus with end organ damage, connective tissue disease",
           "Transfer (> 24 hours) from another hospital",
           "Definite myocardial infarction",
           "Age",
           "Sex with birth",
           "Years of eduction",
           "Support model estimate of the probability of surviving 2 months",
           "DASI (Duke Activity Status Index)",
           "APACHE score",
           "Glasgow Coma Score",
           "Mean blood pressure",
           "WBC",
           "Heart rate",
           "Respiratory rate",
           "Temperature",
           "PaO2/FIO2 ratio",
           "Albumin",
           "Hematocrit",
           "Bilirubin",
           "Creatinine",
           "Sodium",
           "Potassium",
           "PaCo2",
           "PH",
           "Weight",
           "DNR status on day 1",
           "Respiratory diagnosis",
           "Cardiovascular diagnosis",
           "Neurological diagnosis",
           "Gastrointestinal diagnosis",
           "Renal diagnosis",
           "Metabolic diagnosis",
           "Hematologic diagnosis",
           "Sepsis diagnosis",
           "Trauma diagnosis",
           "Orthopedic diagnosis",
           "Income level: $11-$25k",
           "Income level: $25-$50k",
           "Income level: > $50k",
           "Race is black",
           "Other race",
           "Medical insurance type: Medicare",
           "Medical insurance type: Private & Medicare",
           "Medical insurance type: Medicare & Medicaid",
           "Medical insurance type: Medicaid",
           "Medical insurance type: No insurance",
           "Primary disease category: COPD",
           "Primary disease category: MOSF with Sepsis",
           "Primary disease category: MOSF with Malignancy",
           "Primary disease category: CHF",
           "Primary disease category: Coma",
           "Primary disease category: Cirrhosis",
           "Primary disease category: Lung Cancer",
           "Primary disease category: Colon Cancer",
           "Secondary disease category: MOSF with Sepsis",
           "Secondary disease category: MOSF with Malignancy",
           "Secondary disease category: Coma",
           "Secondary disease category: Cirrhosis",
           "Secondary disease category: Lung Cancer",
           "Secondary disease category: Colon Cancer",
           "Having cancer",
           "Metastatic",
           "Indicator of weight is 0",
           "Log value of days of stay in ICU")


table <- data.frame(Description = descr, `No RHC` = mean.c, `RHC` = mean.t)
rownames(table) <- V.name[-73]

xtable(table)

## Design and Analysis

################ Weighting for Paper ################

### PS and outcome models
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
RHC$swang2 <- ifelse(RHC$swang1=="RHC", "No RHC", "RHC")
atc.mult <- swang2 ~ cardiohx + chfhx + dementhx + psychhx + chrpulhx + renalhx + liverhx + 
  gibledhx + malighx + immunhx + transhx + amihx + age + sex + edu + surv2md1 + das2d3pc + 
  aps1 + scoma1 + meanbp1 + wblc1 + hrt1 + resp1 + temp1 + pafi1 + alb1 + hema1 + bili1 + 
  crea1 + sod1 + pot1 + paco21 + ph1 + wtkilo1 + dnr1 + resp + card + neuro + gastr +
  renal + meta + hema + seps + trauma + ortho + income1 + income2 + income3 + raceblack + 
  raceother + ins_care + ins_pcare + ins_carecaid + ins_caid + ins_no + cat1_copd + 
  cat1_mosfsep + cat1_mosfmal + cat1_chf + cat1_coma + cat1_cirr + cat1_lung + cat1_colon +
  cat2_mosfsep + cat2_mosfmal + cat2_coma + cat2_cirr + cat2_lung + cat2_colon + ca_yes +
  ca_meta + wt0

ps.out.form <- icusty ~ cardiohx + chfhx + dementhx + psychhx + chrpulhx + renalhx + liverhx + 
  gibledhx + malighx + immunhx + transhx + amihx + age + sex + edu + surv2md1 + das2d3pc + 
  aps1 + scoma1 + meanbp1 + wblc1 + hrt1 + resp1 + temp1 + pafi1 + alb1 + hema1 + bili1 + 
  crea1 + sod1 + pot1 + paco21 + ph1 + wtkilo1 + dnr1 + resp + card + neuro + gastr +
  renal + meta + hema + seps + trauma + ortho + income1 + income2 + income3 + raceblack + 
  raceother + ins_care + ins_pcare + ins_carecaid + ins_caid + ins_no + cat1_copd + 
  cat1_mosfsep + cat1_mosfmal + cat1_chf + cat1_coma + cat1_cirr + cat1_lung + cat1_colon +
  cat2_mosfsep + cat2_mosfmal + cat2_coma + cat2_cirr + cat2_lung + cat2_colon + ca_yes +
  ca_meta + wt0


### Unadjusted
data.frame(
  Unadjusted = mean(RHC$icusty[RHC$swang1=="RHC"])-mean(RHC$icusty[RHC$swang1=="No RHC"]),
  `Std Error` = sqrt(var(RHC$icusty[RHC$swang1=="RHC"])/length(RHC$icusty[RHC$swang1=="RHC"]) + 
                       var(RHC$icusty[RHC$swang1=="No RHC"])/length(RHC$icusty[RHC$swang1=="No RHC"]))) 


### ATE (IPW)
IPW.wt <- PSweight(ps.formula = ps.mult, yname = "icusty", data = RHC, weight = "IPW")
IPW.wt <- summary(IPW.wt, type = "DIF")$estimates

IPW.aug <- PSweight(ps.formula = ps.mult, yname = "icusty", data = RHC, weight = "IPW",
                    augmentation = TRUE, out.formula = ps.out.form, family = "gaussian")
IPW.aug <- summary(IPW.aug, type = "DIF")$estimates


### IPW with 0.05 trimmed
rhc_trim <- PStrim(data = RHC, ps.formula = ps.mult, delta = 0.05)
rhc_trim

nrow(rhc_trim$data)
nrow(rhc_trim$data)/nrow(RHC)

IPW_wt_trim <- PSweight(ps.formula = ps.mult, yname = "icusty", data = rhc_trim$data, weight = "IPW")
IPW_wt_trim <- summary(IPW_wt_trim, type = "DIF")$estimates

IPW_aug_trim <- PSweight(ps.formula = ps.mult, yname = "icusty", data = rhc_trim$data, weight = "IPW",
                          augmentation = TRUE, out.formula = ps.out.form, family = "gaussian")
IPW_aug_trim <- summary(IPW_aug_trim, type = "DIF")$estimates


### IPW with 0.1 trimmed
rhc_trim1 <- PStrim(data = RHC, ps.formula = ps.mult, delta = 0.1)
rhc_trim1

nrow(rhc_trim1$data)
nrow(rhc_trim1$data)/nrow(RHC)

IPW_wt_trim1 <- PSweight(ps.formula = ps.mult, yname = "icusty", data = rhc_trim1$data, weight = "IPW")
IPW_wt_trim1 <- summary(IPW_wt_trim1, type = "DIF")$estimates

IPW_aug_trim1 <- PSweight(ps.formula = ps.mult, yname = "icusty", data = rhc_trim1$data, weight = "IPW",
                          augmentation = TRUE, out.formula = ps.out.form, family = "gaussian")
IPW_aug_trim1 <- summary(IPW_aug_trim1, type = "DIF")$estimates


### IPW with 0.15 trimmed
rhc_trim2 <- PStrim(data = RHC, ps.formula = ps.mult, delta = 0.15)
rhc_trim2

nrow(rhc_trim2$data)
nrow(rhc_trim2$data)/nrow(RHC)

IPW_wt_trim2 <- PSweight(ps.formula = ps.mult, yname = "icusty", data = rhc_trim2$data, weight = "IPW")
IPW_wt_trim2 <- summary(IPW_wt_trim2, type = "DIF")$estimates

IPW_aug_trim2 <- PSweight(ps.formula = ps.mult, yname = "icusty", data = rhc_trim2$data, weight = "IPW",
                          augmentation = TRUE, out.formula = ps.out.form, family = "gaussian")
IPW_aug_trim2 <- summary(IPW_aug_trim2, type = "DIF")$estimates


### ATO
ATO.wt <- PSweight(ps.formula = ps.mult, yname = "icusty", data = RHC, weight = "overlap")
ATO.wt <- summary(ATO.wt, type = "DIF")$estimates

ATO.aug <- PSweight(ps.formula = ps.mult, yname = "icusty", data = RHC, weight = "overlap",
                    augmentation = TRUE, out.formula = ps.out.form, family = "gaussian")
ATO.aug <- summary(ATO.aug, type = "DIF")$estimates


### ATM
ATM.wt <- PSweight(ps.formula = ps.mult, yname = "icusty", data = RHC, weight = "matching")
ATM.wt <- summary(ATM.wt, type = "DIF")$estimates

ATM.aug <- PSweight(ps.formula = ps.mult, yname = "icusty", data = RHC, weight = "matching",
                    augmentation = TRUE, out.formula = ps.out.form, family = "gaussian")
ATM.aug <- summary(ATM.aug, type = "DIF")$estimates


### ATEN
ATEN.wt <- PSweight(ps.formula = ps.mult, yname = "icusty", data = RHC, weight = "entropy")
ATEN.wt <- summary(ATEN.wt, type = "DIF")$estimates

ATEN.aug <- PSweight(ps.formula = ps.mult, yname = "icusty", data = RHC, weight = "entropy",
                     augmentation = TRUE, out.formula = ps.out.form, family = "gaussian")
ATEN.aug <- summary(ATEN.aug, type = "DIF")$estimates


### ATT
ATT.wt <- PSweight(ps.formula = ps.mult, yname = "icusty", data = RHC, weight = "treated")
ATT.wt <- summary(ATT.wt, type = "DIF")$estimates
  
ATT.aug <- PSweight(ps.formula = ps.mult, yname = "icusty", data = RHC, weight = "treated",
                    augmentation = TRUE, out.formula = ps.out.form, family = "gaussian")
ATT.aug <- summary(ATT.aug, type = "DIF")$estimates


### ATC
ATC.wt <- PSweight(ps.formula = atc.mult, yname = "icusty", data = RHC, weight = "treated")
ATC.wt <- summary(ATC.wt, type = "DIF")$estimates

ATC.aug <- PSweight(ps.formula = atc.mult, yname = "icusty", data = RHC, weight = "treated",
                    augmentation = TRUE, out.formula = ps.out.form, family = "gaussian")
ATC.aug <- summary(ATC.aug, type = "DIF")$estimates


## Save data
n <- nrow(RHC)
rhc_trim <- PStrim(data = RHC, ps.formula = ps.mult, delta = 0.05)
n1 <- nrow(rhc_trim$data)
rhc_trim <- PStrim(data = RHC, ps.formula = ps.mult, delta = 0.1)
n2 <- nrow(rhc_trim$data)
rhc_trim <- PStrim(data = RHC, ps.formula = ps.mult, delta = 0.15)
n3 <- nrow(rhc_trim$data)

save(file="RHC_Results.RData", 
     RHC, p.RHC, n, n1, n2, n3,
     IPW.wt, IPW.aug, IPW_wt_trim, IPW_aug_trim,
     IPW_wt_trim1, IPW_aug_trim1, IPW_wt_trim2, IPW_aug_trim2,
     ATO.wt, ATO.aug, ATM.wt, ATM.wt, ATEN.wt, ATEN.aug, 
     ATT.wt, ATT.aug, ATC.wt, ATC.aug)


################ Var DR ATT/ATC Paper ################ 
source("newSand_func.R")
source("wild_boot_func.R")

library(dplyr)

### PS and outcome models
z <- ifelse(RHC$swang1 == "No RHC", 0, 1)
y <- RHC$icusty
V.name <- c("cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx", "renalhx", "liverhx", "gibledhx", 
            "malighx", "immunhx", "transhx", "amihx", "age", "sex", "edu", "surv2md1", 
            "das2d3pc", "aps1", "scoma1", "meanbp1", "wblc1", "hrt1", "resp1", "temp1", "pafi1", "alb1", 
            "hema1", "bili1", "crea1", "sod1", "pot1", "paco21", "ph1", "wtkilo1", "dnr1", "resp", "card", 
            "neuro", "gastr", "renal", "meta", "hema", "seps", "trauma", "ortho", "income1", "income2", 
            "income3", "raceblack", "raceother", "ins_care", "ins_pcare", "ins_carecaid", "ins_caid", 
            "ins_no", "cat1_copd", "cat1_mosfsep", "cat1_mosfmal", "cat1_chf", "cat1_coma", "cat1_cirr", 
            "cat1_lung", "cat1_colon", "cat2_mosfsep", "cat2_mosfmal", "cat2_coma", "cat2_cirr", 
            "cat2_lung", "cat2_colon", "ca_yes", "ca_meta", "wt0")
X <- apply(RHC[, V.name], 2, as.numeric)
X.out <- X

### Asymptotic variance: ATE, ATC, ATT
result <- ATE(y = y, z = z, X = X, DR = TRUE, X.out = X.out)
ATE.sand <- data.frame(Est = result$tau, SD = result$se, p = 2*(1-pnorm(abs(result$tau/result$se))))

result <- ATC(y = y, z = z, X = X, DR = TRUE, X.out = X.out)
ATC.sand <- data.frame(Est = result$tau, SD = result$se, p = 2*(1-pnorm(abs(result$tau/result$se))))

result <- ATT(y = y, z = z, X = X, DR = TRUE, X.out = X.out) 
ATT.sand <- data.frame(Est = result$tau, SD = result$se, p = 2*(1-pnorm(abs(result$tau/result$se))))


### Wild bootstrap variance: ATE, ATC, ATT
result <- ATE.boot(y = y, z = z, X = X, X.out = X.out, R = 1000, RV = "Rad") 
ATE.wb1 <- data.frame(Est = result$PE, SD = result$Std.Boot, p = 2*(1-pnorm(abs(result$PE/result$Std.Boot))))

result <- ATE.boot(y = y, z = z, X = X, X.out = X.out, R = 1000, RV = "Exp") 
ATE.wb2 <- data.frame(Est = result$PE, SD = result$Std.Boot, p = 2*(1-pnorm(abs(result$PE/result$Std.Boot))))

result <- ATC.boot(y = y, z = z, X = X, X.out = X.out, R = 1000, RV = "Rad") 
ATC.wb1 <- data.frame(Est = result$PE, SD = result$Std.Boot, p = 2*(1-pnorm(abs(result$PE/result$Std.Boot))))

result <- ATC.boot(y = y, z = z, X = X, X.out = X.out, R = 1000, RV = "Exp")
ATC.wb2 <- data.frame(Est = result$PE, SD = result$Std.Boot, p = 2*(1-pnorm(abs(result$PE/result$Std.Boot))))

result <- ATT.boot(y = y, z = z, X = X, X.out = X.out, R = 1000, RV = "Rad") 
ATT.wb1 <- data.frame(Est = result$PE, SD = result$Std.Boot, p = 2*(1-pnorm(abs(result$PE/result$Std.Boot))))

result <- ATT.boot(y = y, z = z, X = X, X.out = X.out, R = 1000, RV = "Exp") 
ATT.wb2 <- data.frame(Est = result$PE, SD = result$Std.Boot, p = 2*(1-pnorm(abs(result$PE/result$Std.Boot))))

save(file = "RHC_VarDR.RData", ATE.sand, ATC.sand, ATT.sand, 
     ATE.wb1, ATE.wb2, ATC.wb1, ATC.wb2, ATT.wb1, ATT.wb2)

### TMLE
library(tmle)
z <- ifelse(RHC$swang1 == "No RHC", 0, 1)
y <- RHC$icusty
X <- apply(RHC[, V.name], 2, as.numeric)
ps.mult <- A ~ cardiohx + chfhx + dementhx + psychhx + chrpulhx + renalhx + liverhx + 
  gibledhx + malighx + immunhx + transhx + amihx + age + sex + edu + surv2md1 + das2d3pc + 
  aps1 + scoma1 + meanbp1 + wblc1 + hrt1 + resp1 + temp1 + pafi1 + alb1 + hema1 + bili1 + 
  crea1 + sod1 + pot1 + paco21 + ph1 + wtkilo1 + dnr1 + resp + card + neuro + gastr +
  renal + meta + hema + seps + trauma + ortho + income1 + income2 + income3 + raceblack + 
  raceother + ins_care + ins_pcare + ins_carecaid + ins_caid + ins_no + cat1_copd + 
  cat1_mosfsep + cat1_mosfmal + cat1_chf + cat1_coma + cat1_cirr + cat1_lung + cat1_colon +
  cat2_mosfsep + cat2_mosfmal + cat2_coma + cat2_cirr + cat2_lung + cat2_colon + ca_yes +
  ca_meta + wt0

ps.out.form <- Y ~ cardiohx + chfhx + dementhx + psychhx + chrpulhx + renalhx + liverhx + 
  gibledhx + malighx + immunhx + transhx + amihx + age + sex + edu + surv2md1 + das2d3pc + 
  aps1 + scoma1 + meanbp1 + wblc1 + hrt1 + resp1 + temp1 + pafi1 + alb1 + hema1 + bili1 + 
  crea1 + sod1 + pot1 + paco21 + ph1 + wtkilo1 + dnr1 + resp + card + neuro + gastr +
  renal + meta + hema + seps + trauma + ortho + income1 + income2 + income3 + raceblack + 
  raceother + ins_care + ins_pcare + ins_carecaid + ins_caid + ins_no + cat1_copd + 
  cat1_mosfsep + cat1_mosfmal + cat1_chf + cat1_coma + cat1_cirr + cat1_lung + cat1_colon +
  cat2_mosfsep + cat2_mosfmal + cat2_coma + cat2_cirr + cat2_lung + cat2_colon + ca_yes +
  ca_meta + wt0

tmle(Y=y, A=z, W=X, Qform=ps.out.form, family = "gaussian", gform=ps.mult)

### Standard bootstrap
source("newSand_func.R")
z <- ifelse(RHC$swang1 == "No RHC", 0, 1)
y <- RHC$icusty
X <- apply(RHC[, V.name], 2, as.numeric)
R <- 1000 # bootstrap replications
ATE.boot <- c()
ATT.boot <- c()
ATC.boot <- c()
n <- length(RHC$sadmdte)
for (i in 1:R) {
  RHC.boot <- RHC[sample(1:n, size=n, replace=T),]
  z <- ifelse(RHC.boot$swang1 == "No RHC", 0, 1)
  y <- RHC.boot$icusty
  X <- apply(RHC.boot[, V.name], 2, as.numeric)
  ATE.boot <- c(ATE.boot, ATE.PE(y=y, z=z, X=X, DR=TRUE, X.out=X)$tau)
  ATT.boot <- c(ATT.boot, ATT.PE(y=y, z=z, X=X, DR=TRUE, X.out=X)$tau)
  ATC.boot <- c(ATC.boot, ATC.PE(y=y, z=z, X=X, DR=TRUE, X.out=X)$tau)
}

# bootstrap variance and p-values
z <- ifelse(RHC$swang1 == "No RHC", 0, 1)
y <- RHC$icusty
X <- apply(RHC[, V.name], 2, as.numeric)
X.out <- X
ATT <- ATT.PE(y=y, z=z, X=X, DR=TRUE, X.out=X)$tau
ATC <- ATC.PE(y=y, z=z, X=X, DR=TRUE, X.out=X)$tau

sd(ATT.boot)
2*(1-pnorm(abs(ATT)/sd(ATT.boot)))
sd(ATC.boot)
2*(1-pnorm(abs(ATC)/sd(ATC.boot)))
