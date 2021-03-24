library(tidyverse)
library(reshape2)
library(ggpubr)
library(lubridate)

myfun = function(x) {
  #as.numeric(x[1]) / as.numeric(x[2]) ## Proportion
  as.numeric(x[2]) - as.numeric(x[1]) ## Number mismatches
}

## Read data
chim_file = "./data/Chimerism prediction of relapse internal copy v8 g_mtx_mmf_ci_abd merged.csv"
dat = read.csv(chim_file)
names(dat)

## Add ID column
dat$ID = as.factor(seq(1:nrow(dat)))

## Convert hla to proportion
dat$HLA2 <- unlist(lapply(strsplit(dat$HLA, "/"), myfun))

## ------------------------------------------------------------------------- ##
## BMC values
## Extract dates
bmc.df = dat[,c("ID", names(dat[,grep(glob2rx("D*MC"), names(dat))]))]

## Extract 'constants'
bmc.df$dot = dat$DOT
bmc.df$dor = dat$DOR
bmc.df$dor[bmc.df$dor == ""] <- "7/1/20"
bmc.df$txage = dat$TXAGE
bmc.df$relapse = dat$R
bmc.df$sex = dat$G
bmc.df$rstatprtx = dat$RSTATPRTX
#bmc.df$g = dat$G.1
bmc.df$hla = dat$HLA2
bmc.df$tbi = dat$TBI
#bmc.df$abd = dat$aBD
#bmc.df$ci = dat$CI
#bmc.df$mtx = dat$MTX
#bmc.df$mmf = dat$MMF
#bmc.df$e = dat$E
bmc.df$abd_ci_mtx_mmi = dat$G.HLA.GVHD.prophy
bmc.df$agvhd = dat$aGVHD
bmc.df$cgvhd = dat$cGVHD

## Melt to long format
bmc.df = melt(bmc.df, id.vars = c("ID", "dot", "dor", "txage", "relapse",
                                  "sex", "rstatprtx", "hla", "tbi",
                                  "abd_ci_mtx_mmi", "agvhd", "cgvhd"), 
                variable.name = "test", value.name = "bdate")
bmc.df$bdate = mdy(bmc.df$bdate)
bmc.df$dot = mdy(bmc.df$dot)
bmc.df$dor = mdy(bmc.df$dor)
bmc.df$btime = bmc.df$bdate - bmc.df$dot
bmc.df$rtime = bmc.df$dor - bmc.df$dot

## Binary relapse variable for modeling
bmc.df$rbin = factor(bmc.df$relapse)

## Get chimerism data
bmc_cdw = dat[,c("ID", names(dat[,grep(glob2rx("*MCW"), names(dat))]))]
bmc_cd3 = dat[,c("ID", names(dat[,grep(glob2rx("*MC3"), names(dat))]))]
bmc_cd3 = bmc_cd3[, -ncol(bmc_cd3)]
bmc_cd15 = dat[,c("ID", names(dat[,grep(glob2rx("*MC15"), names(dat))]))]
bmc_cd34 = dat[,c("ID", names(dat[,grep(glob2rx("*MC34"), names(dat))]))]

## Melting
bmc_cdw = melt(bmc_cdw, id.vars = "ID", variable.name = "test", value.name = "prop")
bmc_cd3 = melt(bmc_cd3, id.vars = "ID", variable.name = "test", value.name = "prop")
bmc_cd15 = melt(bmc_cd15, id.vars = "ID", variable.name = "test", value.name = "prop")
bmc_cd34 = melt(bmc_cd34, id.vars = "ID", variable.name = "test", value.name = "prop")

bmc.df$cdw = bmc_cdw$prop
bmc.df$cd3 = bmc_cd3$prop
bmc.df$cd15 = bmc_cd15$prop
bmc.df$cd34 = bmc_cd34$prop

write.csv(bmc.df, "./data/bmc.df.csv", row.names = FALSE)

## Peripheral blood

dat = read.csv(chim_file)

names(dat)

## Add ID column
dat$ID = as.factor(seq(1:nrow(dat)))

## Convert hla to proportion
dat$HLA2 <- unlist(lapply(strsplit(dat$HLA, "/"), myfun))

## Extract dates
pbc.df = dat[,c("ID", names(dat[,grep(glob2rx("D*PBC"), names(dat))]))]

## Extract 'constants'
pbc.df$dot = dat$DOT
pbc.df$dor = dat$DOR
pbc.df$dor[pbc.df$dor == ""] <- "7/1/20"
pbc.df$txage = dat$TXAGE
pbc.df$relapse = dat$R
pbc.df$sex = dat$G
pbc.df$rstatprtx = dat$RSTATPRTX
#pbc.df$g = dat$G.1
pbc.df$hla = dat$HLA2
pbc.df$tbi = dat$TBI
#pbc.df$abd = dat$aBD
#pbc.df$ci = dat$CI
#pbc.df$mtx = dat$MTX
#pbc.df$mmf = dat$MMF
#pbc.df$e = dat$E
pbc.df$abd_ci_mtx_mmi = dat$G.HLA.GVHD.prophy
pbc.df$agvhd = dat$aGVHD
pbc.df$cgvhd = dat$cGVHD

## Melt to long format
pbc.df = melt(pbc.df, id.vars = c("ID", "dot", "dor", "txage", "relapse",
                                  "sex", "rstatprtx", "hla", "tbi",
                                  "abd_ci_mtx_mmi", "agvhd", "cgvhd"), 
              variable.name = "test", value.name = "pdate")
pbc.df$pdate = mdy(pbc.df$pdate)
pbc.df$dot = mdy(pbc.df$dot)
pbc.df$dor = mdy(pbc.df$dor)
pbc.df$ptime = pbc.df$pdate - pbc.df$dot
pbc.df$rtime = pbc.df$dor - pbc.df$dot

## Binary relapse variable for modeling
pbc.df$rbin = factor(pbc.df$relapse)

## Get chimerism data
pbc_cdw = dat[,c("ID", names(dat[,grep(glob2rx("*PBCW"), names(dat))]))]
pbc_cd3 = dat[,c("ID", names(dat[,grep(glob2rx("*PBC3"), names(dat))]))]
pbc_cd3 = pbc_cd3[, -ncol(pbc_cd3)]
pbc_cd15 = dat[,c("ID", names(dat[,grep(glob2rx("*PBC15"), names(dat))]))]
pbc_cd34 = dat[,c("ID", names(dat[,grep(glob2rx("*PBC34"), names(dat))]))]

## Melting
pbc_cdw = melt(pbc_cdw, id.vars = "ID", variable.name = "test", value.name = "prop")
pbc_cd3 = melt(pbc_cd3, id.vars = "ID", variable.name = "test", value.name = "prop")
pbc_cd15 = melt(pbc_cd15, id.vars = "ID", variable.name = "test", value.name = "prop")
pbc_cd34 = melt(pbc_cd34, id.vars = "ID", variable.name = "test", value.name = "prop")

pbc.df$cdw = pbc_cdw$prop
pbc.df$cd3 = pbc_cd3$prop
pbc.df$cd15 = pbc_cd15$prop
pbc.df$cd34 = pbc_cd34$prop

write.csv(pbc.df, "./data/pbc.df.csv", row.names = FALSE)

all.df = bmc.df
all.df$bmc_cdw = all.df$cdw
all.df$bmc_cd3 = all.df$cd3
all.df$bmc_cd15 = all.df$cd15
all.df$bmc_cd34 = all.df$cd34
all.df = all.df %>% 
  select(-cdw, -cd3, -cd15, -cd34)

all.df$pdate = pbc.df$pdate
all.df$ptime = pbc.df$ptime
all.df$pbc_cdw = pbc.df$cdw
all.df$pbc_cd3 = pbc.df$cd3
all.df$pbc_cd15 = pbc.df$cd15
all.df$pbc_cd34 = pbc.df$cd34

write.csv(all.df, "./data/all.df.csv", row.names = FALSE)

