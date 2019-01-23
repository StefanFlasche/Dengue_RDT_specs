
# Investigating Dengue RDT sensitivity and specificty requirements
# written by Stefan Flasche using R version 3.4.4

# load packages
require(tidyverse)

# data: cumulative incidence over 5yrs observed in the trial and its stratification into 
#       seropositive and seronegative using PRNT50 with multiple imputation (in per 100). 
df <- tibble(no=1:8,
             randomisation = c("vacc","vacc","control","control","vacc","vacc","control","control"),
             serostatus = c("pos","neg","pos","neg","pos","neg","pos","neg"),
             outcome = c(rep("hosp",4),rep("severe",4)),
             incidence.ph.mid = c(.375,1.571,1.883,1.093,0.075,0.404,0.48,0.174),
             incidence.ph.lo = c(.263,1.125,1.536,.526,0.034,0.218,0.335,0.036),
             incidence.ph.hi = c(.535,2.193,2.307,2.265,0.165,0.749,0.688,0.834))

# plot data
df %>% ggplot(aes(x= serostatus, y = incidence.ph.mid, ymin = incidence.ph.lo, ymax = incidence.ph.hi,color=randomisation)) +
  geom_linerange(position=position_dodge(width = 0.5)) +
  geom_point(position=position_dodge(width = 0.5)) +
  facet_grid(.~outcome, scales = "free") +
  coord_flip() + ylab("incidence per 100")

# sample from lognormal distribution in which the 50%, 2.5% and 97.5% quanitles fit the observed mean, CI.lo and CI.hi respectively
LnfitSample <- function(input= c(mean=1, lo=.5, hi=2), N=1000){
  mean = input[1] %>% as.numeric()
  lo = input[2] %>% as.numeric()
  hi = input[3] %>% as.numeric()
  rlnorm(N, meanlog = log(mean), sdlog = (log(mean/lo) +  log(hi/mean))/2/1.92) %>%
    return()
}

# calculate population impact of test and vaccinate strategy
## need to include resampling of incidence 
CasesAverted <- function(seroPrevalence = .7, sensitivity = 1, specificity = 0, 
                         cohortSize=100000, df.tmp = df, outcome = "hosp"){
  df.tmp = df.tmp %>% filter(outcome == "hosp")
  Inc.SeroPos.Vacc <- df.tmp %>% filter(randomisation=="vacc" & serostatus =="pos") %>% 
    select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
  Inc.SeroPos.Cont <- df.tmp %>% filter(randomisation=="control" & serostatus =="pos") %>% 
    select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
  Inc.SeroNeg.Vacc <- df.tmp %>% filter(randomisation=="vacc" & serostatus =="neg") %>% 
    select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
  Inc.SeroNeg.Cont <- df.tmp %>% filter(randomisation=="control" & serostatus =="neg") %>% 
    select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
  
  CasesAvertedSeroPos = cohortSize * seroPrevalence * sensitivity * (Inc.SeroPos.Cont - Inc.SeroPos.Vacc) / 100
  CasesAvertedSeroNeg = cohortSize * (1-seroPrevalence) * (1 - specificity) * (Inc.SeroNeg.Cont - Inc.SeroNeg.Vacc) /100
  
  df_res <- tibble(seroPrevalence = seroPrevalence, sensitivity = sensitivity,
                   specificity = specificity, outcome = outcome, 
                   CasesAvertedSeroPos.mid = median(CasesAvertedSeroPos), 
                   CasesAvertedSeroPos.lo = quantile(CasesAvertedSeroPos,.025), 
                   CasesAvertedSeroPos.hi = quantile(CasesAvertedSeroPos,.975),
                   CasesAvertedSeroNeg.mid = median(CasesAvertedSeroNeg), 
                   CasesAvertedSeroNeg.lo = quantile(CasesAvertedSeroNeg,.025), 
                   CasesAvertedSeroNeg.hi = quantile(CasesAvertedSeroNeg,.975), 
                   CasesAvertedTotal.mid =  median(CasesAvertedSeroPos + CasesAvertedSeroNeg), 
                   CasesAvertedTotal.lo = quantile(CasesAvertedSeroPos + CasesAvertedSeroNeg, .025),
                   CasesAvertedTotal.hi = quantile(CasesAvertedSeroPos + CasesAvertedSeroNeg, .975))
  return(df_res)
}

# plot impact of test and vaccinate strategy
df.plt = NULL
for(seroPrevalence in c(.6,.7,.8)){
  for(specificity in c(.9,.95,.98,.99,1)){
    for(sensitivity in seq(.8,1,by=.05)){
      df.plt <- df.plt %>% 
        rbind(CasesAverted(seroPrevalence = seroPrevalence, sensitivity = sensitivity, specificity = specificity))
    }
  }
}

df.plt <- df.plt %>%
  gather(key,value, -seroPrevalence, -sensitivity, -specificity, -outcome) %>%
  extract(key, c("outcome","conf"),"(.*)\\.(.*)") %>%
  spread(conf, value)
  
df.plt %>% ggplot(aes(x=sensitivity, y = mid, ymin=lo, ymax= hi, color = as.factor(seroPrevalence))) +
  geom_linerange(position = position_dodge(width = 0.04)) +
  geom_point(position = position_dodge(width = 0.04)) +
  facet_grid(outcome ~ specificity, scales = "free") 

  
# calculate sensitivity needed to avoid chosing to vaccinate test negative
SensitvityNeeded <- function(){
  
}

# calculate test sens and spec given specified risk threshold and serprevalence range

