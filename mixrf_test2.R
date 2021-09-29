## ----setup, include=FALSE------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------------------------------------------
set.seed(1234)
source("MixRFb.R")

## ----message=FALSE-------------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(ggpubr)
library(randomForest)
library(vip)
library(pdp)
library(MixRF)
library(caret)
library(party)
library(lme4)
library(PresenceAbsence)

library(foreach)

## ------------------------------------------------------------------------------------------
all.df <- read.csv("./data/all.df.csv")


## ------------------------------------------------------------------------------------------
all.df$dot <- ymd(all.df$dot)
all.df$dor <- ymd(all.df$dor)
all.df$bdate <- ymd(all.df$bdate)
all.df$pdate <- ymd(all.df$pdate)


## ------------------------------------------------------------------------------------------
all.df <- all.df %>% mutate_if(is.character,as.factor)


## ------------------------------------------------------------------------------------------
all.df$rbin <- factor(all.df$rbin, levels = c("yes", "no"))


## ------------------------------------------------------------------------------------------
all.df$time <- as.numeric(substr(all.df$test, 2, 2))


## ------------------------------------------------------------------------------------------
all.df <- all.df[which(all.df$bdate < all.df$dor | is.na(all.df$dor)), ]


## ------------------------------------------------------------------------------------------
all.df <- all.df[which(all.df$rbin == "no" | all.df$rtime < 720),]


## ------------------------------------------------------------------------------------------
all.df <- all.df[!is.na(all.df$bmc_cdw) & !is.na(all.df$bmc_cd3) & 
                   !is.na(all.df$bmc_cd15) & !is.na(all.df$bmc_cd34) &
                   !is.na(all.df$pbc_cdw) & !is.na(all.df$pbc_cd3) & 
                   !is.na(all.df$pbc_cd15) & !is.na(all.df$pbc_cd34),]
all.df <<- all.df


## ------------------------------------------------------------------------------------------
f1 = rbin ~ txage + sex + rstatprtx + g + hla + 
  tbi + abd + ci + mtx + mmf + e + agvhd + cgvhd +
  bmc_cdw + bmc_cd3 + bmc_cd15 + bmc_cd34 + 
  pbc_cdw + pbc_cd3 + pbc_cd15 + pbc_cd34
f1 = rbin ~ txage + sex + rstatprtx + hla + 
  tbi + gmgp + agvhd + cgvhd +
  bmc_cdw + bmc_cd3 + bmc_cd15 + bmc_cd34 + 
  pbc_cdw + pbc_cd3 + pbc_cd15 + pbc_cd34

## Full model
Y <- all.df %>%
  select(rbin) %>%
  mutate(rbin = as.numeric(as.factor(rbin)) - 1)
X <- all.df %>%
  select(ID, txage, rbin, sex, rstatprtx, hla, 
         tbi, gmgp, agvhd, cgvhd, 
         bmc_cdw, bmc_cd3, bmc_cd15, bmc_cd34,
         pbc_cdw, pbc_cd3, pbc_cd15, pbc_cd34)

all.df$rbin2 <- as.numeric(all.df$rbin) - 1
f1 <- "txage + sex + rstatprtx + hla + 
  tbi + gmgp + agvhd + cgvhd +
  bmc_cdw + bmc_cd3 + bmc_cd15 + bmc_cd34 + 
  pbc_cdw + pbc_cd3 + pbc_cd15 + pbc_cd34"
# merf1 <- MixRFb(all.df$rbin2, f1, random = "(1 | ID)", 
#                data = all.df, verbose = TRUE)

## ------------------------------------------------------------------------------------------
## Cross validation
patients = unique(sort(all.df$ID))
npatients = length(patients)
nfolds = 5
pfolds = createFolds(all.df$relapse, k = nfolds)
pfolds


## ------------------------------------------------------------------------------------------
mod_results <- data.frame(folds = c(1:nfolds, "mean"),
                          auc = rep(NA, nfolds+1),
                          sens = rep(NA, nfolds+1),
                          spec = rep(NA, nfolds+1)
)


## ------------------------------------------------------------------------------------------
for (i in 1:nfolds) {
  print(i)
  patient_id <- pfolds[[i]]
  training <- all.df[-patient_id,]
  testing <-  all.df[ patient_id,]
  
  ## Select predictors
  train_x <- training %>% 
    select(ID, txage, rbin2, sex, rstatprtx, hla, 
           tbi, gmgp, agvhd, cgvhd, 
           bmc_cdw, bmc_cd3, bmc_cd15, bmc_cd34,
           pbc_cdw, pbc_cd3, pbc_cd15, pbc_cd34)
  
  test_x <- testing %>% 
    select(ID, txage, rbin2, sex, rstatprtx, hla, 
           tbi, gmgp, agvhd, cgvhd, 
           bmc_cdw, bmc_cd3, bmc_cd15, bmc_cd34,
           pbc_cdw, pbc_cd3, pbc_cd15, pbc_cd34)
  
  control<-list(se = TRUE, ntrees = 500, R = 100)
  
  # merf1 <- MixRFb(train_x$rbin2, f1, 
  #                 random = "(1 | ID)", 
  #                 data = train_x, 
  #                 ErrorTolerance=0.1, MaxIterations=50,
  #                 ErrorTolerance0=0.1, MaxIterations0=5, 
  #                 verbose = TRUE)
  # 
  # pred_test <- predict.MixRF(merf1, test_x)  ## Log odds
  # pred_test <- exp(pred_test) / (1 + exp(pred_test)) ## Prob
  pred_test <- runif(nrow(test_x))
  
  df <- data.frame(id = seq(1:nrow(testing)),
                   obs = as.numeric(testing$rbin)-1,
                   pred = pred_test)
  
  mod_results$auc[i] <- auc(df)$AUC[[1]]
  ## Get threshold
  opt_thresh <- optimal.thresholds(df, opt.methods = 3)[2]
  opt_cmx <- cmx(df, threshold = opt_thresh$pred)
  mod_results$sens[i] = sensitivity(opt_cmx)$sensitivity[[1]]
  mod_results$spec[i] = specificity(opt_cmx)$specificity[[1]]
  
}


## ------------------------------------------------------------------------------------------
mod_results[6, 2:4] = apply(mod_results[, 2:4], 2, mean, na.rm = TRUE)
knitr::kable(mod_results, digits = 4)

