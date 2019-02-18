
# Investigating Dengue RDT sensitivity and specificty requirements
# written by Stefan Flasche using R version 3.4.4

# load packages
require(tidyverse)


# Data --------------------------------------------------------------------

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
p.data <- df %>% ggplot(aes(x= serostatus, y = incidence.ph.mid, ymin = incidence.ph.lo, ymax = incidence.ph.hi,color=randomisation)) +
  geom_linerange(position=position_dodge(width = 0.5)) +
  geom_point(position=position_dodge(width = 0.5)) +
  facet_grid(.~outcome, scales = "free") +
  coord_flip() + ylab("incidence per 100")
ggsave(filename = "Pics\\Fig1_Data.tiff",p.data ,unit="cm", width = 14, height = 5, compression = "lzw", dpi = 300)


# sample from lognormal distribution in which the 50%, 2.5% and 97.5% quanitles fit the observed mean, CI.lo and CI.hi respectively
LnfitSample <- function(input= c(mean=1, lo=.5, hi=2), N=10000){
  mean = input[1] %>% as.numeric()
  lo = input[2] %>% as.numeric()
  hi = input[3] %>% as.numeric()
  rlnorm(N, meanlog = log(mean), sdlog = (log(mean/lo) +  log(hi/mean))/2/1.92) %>%
    return()
}


# Population impact of test and vaccinate ---------------------------------------------

# calculate population impact of test and vaccinate strategy
## need to include resampling of incidence 
CasesAverted <- function(seroPrevalence = .7, sensitivity = 1, specificity = 0, 
                         cohortSize=100000, df.tmp = df, outcm = "hosp"){
  df.tmp = df.tmp %>% filter(outcome == outcm)
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
                   specificity = specificity, outcome = outcm, 
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
  
p.tandv = df.plt %>% ggplot(aes(x=sensitivity, y = mid, ymin=lo, ymax= hi, color = as.factor(seroPrevalence))) +
  geom_linerange(position = position_dodge(width = 0.04)) +
  geom_point(position = position_dodge(width = 0.04)) +
  facet_grid(outcome ~ specificity, scales = "free") +
  scale_color_discrete(name="seroprevalence") + ylab("hospitalised dengue cases averted\nin a 100,000 cohort") +
  geom_hline(yintercept = 0, color="black", lty="dashed")
ggsave(filename = "Pics\\Fig2_Impact_TandV.tiff",p.tandv ,unit="cm", width = 25, height = 14, compression = "lzw", dpi = 300)


# Impact of sensitivity --------------------------------------------------------------------
  
# calculate sensitivity needed to avoid chosing to vaccinate test negative
TestNegRatio <- function(seroPrevalence = .7, sensitivity = .8, specificity = .95,
                             df.tmp = df, outcm = "hosp", midonly=F){
  df.tmp = df.tmp %>% filter(outcome == outcm)
  if(midonly){
    Inc.SeroPos.Vacc <- df.tmp %>% filter(randomisation=="vacc" & serostatus =="pos") %>% 
      select(incidence.ph.mid) %>% as.numeric()
    Inc.SeroPos.Cont <- df.tmp %>% filter(randomisation=="control" & serostatus =="pos") %>% 
      select(incidence.ph.mid) %>% as.numeric()
    Inc.SeroNeg.Vacc <- df.tmp %>% filter(randomisation=="vacc" & serostatus =="neg") %>% 
      select(incidence.ph.mid) %>% as.numeric()
    Inc.SeroNeg.Cont <- df.tmp %>% filter(randomisation=="control" & serostatus =="neg") %>% 
      select(incidence.ph.mid) %>% as.numeric()
  }else{
    Inc.SeroPos.Vacc <- df.tmp %>% filter(randomisation=="vacc" & serostatus =="pos") %>% 
      select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
    Inc.SeroPos.Cont <- df.tmp %>% filter(randomisation=="control" & serostatus =="pos") %>% 
      select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
    Inc.SeroNeg.Vacc <- df.tmp %>% filter(randomisation=="vacc" & serostatus =="neg") %>% 
      select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
    Inc.SeroNeg.Cont <- df.tmp %>% filter(randomisation=="control" & serostatus =="neg") %>% 
      select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
  }
   
  risk.if.testNeg.vacc <- ((1-seroPrevalence) * specificity *  Inc.SeroNeg.Vacc + 
    seroPrevalence * (1-sensitivity) * Inc.SeroPos.Vacc ) / 
    ((1-seroPrevalence) * specificity  + seroPrevalence * (1-sensitivity) ) 
  risk.if.testNeg.nova <- ((1-seroPrevalence) * specificity *  Inc.SeroNeg.Cont + 
    seroPrevalence * (1-sensitivity) * Inc.SeroPos.Cont ) / 
    ((1-seroPrevalence) * specificity  + seroPrevalence * (1-sensitivity) ) 
  
  df_res <- tibble(seroPrevalence = seroPrevalence, sensitivity = sensitivity,
                   specificity = specificity, outcome = outcm, 
                   riskiftestNegvacc.mid = median(risk.if.testNeg.vacc), 
                   riskiftestNegnova.mid = median(risk.if.testNeg.nova),
                   riskiftestNegvacc.lo = quantile(risk.if.testNeg.vacc,.025), 
                   riskiftestNegvacc.hi = quantile(risk.if.testNeg.vacc,.975),
                   riskiftestNegnova.lo = quantile(risk.if.testNeg.nova,.025), 
                   riskiftestNegnova.hi = quantile(risk.if.testNeg.nova,.975),
                   RR.mid =  median(risk.if.testNeg.vacc / risk.if.testNeg.nova), 
                   RR.lo = quantile(risk.if.testNeg.vacc / risk.if.testNeg.nova, .025),
                   RR.hi = quantile(risk.if.testNeg.vacc / risk.if.testNeg.nova, .975))
  return(df_res)
}

# analyse where test nagative vacc vs novacc risks are for a number of scenarios
df.plt = NULL
for(seroPrevalence in c(.5,.7,.9)){
  for(specificity in c(.95,1)){
    for(sensitivity in seq(.6,1,by=.1)){
      df.plt <- df.plt %>% 
        rbind(TestNegRatio(seroPrevalence = seroPrevalence, sensitivity = sensitivity, specificity = specificity)) %>%
        rbind(TestNegRatio(seroPrevalence = seroPrevalence, sensitivity = sensitivity, specificity = specificity, outcm = "severe"))
      
    }
  }
}

p.sens <- df.plt %>% 
  ggplot(aes(x=sensitivity, y = RR.mid, ymin=RR.lo, ymax= RR.hi, color=outcome)) +
  geom_linerange(position = position_dodge(width = 0.04)) +
  geom_point(position = position_dodge(width = 0.04)) +
  facet_grid(seroPrevalence ~ specificity, scales = "free") +
  scale_y_log10() + geom_hline(yintercept = 1, color="black", lty="dashed") +
  scale_color_discrete(name="outcome") +
  ylab("RR for hospitalised dengue\nin test-negative for\nvaccination vs no vaccination")
ggsave(filename = "Pics\\FigX_SensitivityImpact.tiff",p.sens ,unit="cm", width = 20, height = 12, compression = "lzw", dpi = 300)

# analyse where test nagative vacc vs novacc risks - map
df.plt2 = NULL
for(seroPrevalence in c(.5,.6,.7,.8,.9)){
  for(specificity in seq(.9,1,by=.01)){
    for(sensitivity in seq(.6,1,by=.01)){
      df.plt2 <- df.plt2 %>% 
        rbind(TestNegRatio(seroPrevalence = seroPrevalence, sensitivity = sensitivity, specificity = specificity, midonly=T)) %>%
        rbind(TestNegRatio(seroPrevalence = seroPrevalence, sensitivity = sensitivity, specificity = specificity, outcm = "severe", midonly=T))
      
    }
  }
}

p.sens.b <- df.plt2 %>% 
  ggplot(aes(x=sensitivity, y=specificity, fill = RR.mid)) +
  facet_grid(seroPrevalence ~ outcome) + 
  geom_tile() +
  scale_fill_gradient2(name="RR in vaccinees\nvs control", low = "green", mid = "yellow",
                       high = "red", midpoint = 0, breaks=round(c(.5,1/1.2,1/1.1,1,1.1,1.2,2),2),trans = "log", limits=c(1/1.21,1.2)) +
  theme_bw() 
ggsave("Pics\\Fig3a_SensitivityImpactMap.tiff", p.sens.b, width = 20, height = 13, units = "cm", compression="lzw", dpi =300)


# importance of NPV ---------------------------------------------
TestNegRatio.NPV <- function(NPV = .7, df.tmp = df, outcm = "hosp"){
  df.tmp = df.tmp %>% filter(outcome == outcm)

  Inc.SeroPos.Vacc <- df.tmp %>% filter(randomisation=="vacc" & serostatus =="pos") %>% 
    select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
  Inc.SeroPos.Cont <- df.tmp %>% filter(randomisation=="control" & serostatus =="pos") %>% 
    select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
  Inc.SeroNeg.Vacc <- df.tmp %>% filter(randomisation=="vacc" & serostatus =="neg") %>% 
    select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
  Inc.SeroNeg.Cont <- df.tmp %>% filter(randomisation=="control" & serostatus =="neg") %>% 
    select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()

  risk.if.testNeg.vacc <- ( NPV * Inc.SeroNeg.Vacc + 
                            (1-NPV) * Inc.SeroPos.Vacc ) 
  risk.if.testNeg.nova <- ( NPV * Inc.SeroNeg.Cont + 
                            (1-NPV) * Inc.SeroPos.Cont ) 
  
  df_res <- tibble(NPV = NPV, outcome = outcm, 
                   riskiftestNegvacc.mid = median(risk.if.testNeg.vacc), 
                   riskiftestNegnova.mid = median(risk.if.testNeg.nova),
                   riskiftestNegvacc.lo = quantile(risk.if.testNeg.vacc,.025), 
                   riskiftestNegvacc.hi = quantile(risk.if.testNeg.vacc,.975),
                   riskiftestNegnova.lo = quantile(risk.if.testNeg.nova,.025), 
                   riskiftestNegnova.hi = quantile(risk.if.testNeg.nova,.975),
                   RR.mid =  median(risk.if.testNeg.vacc / risk.if.testNeg.nova), 
                   RR.lo = quantile(risk.if.testNeg.vacc / risk.if.testNeg.nova, .025),
                   RR.hi = quantile(risk.if.testNeg.vacc / risk.if.testNeg.nova, .975))
  return(df_res)
}

df.plt.NPV = NULL
for( NPV in seq(.5,1,by=0.01)){
  df.plt.NPV = df.plt.NPV %>% 
    rbind(TestNegRatio.NPV(NPV)) %>%
    rbind(TestNegRatio.NPV(NPV, outcm = "severe"))
}
df.plt.NPV %>%
  ggplot(aes(x=NPV, y=RR.mid, ymin=RR.lo, ymax=RR.hi, group=outcome, color=outcome, fill=outcome)) +
  geom_ribbon(alpha=0.4, color=NA) +
  geom_line() +
  geom_hline(yintercept = 1, color = "black", lty = "dashed") + 
  facet_grid(outcome~., scales="free") +
  scale_y_log10() +
  xlab("Negative predictive value") + ylab("RR for dengue disease \nin test-negative for\nvaccination vs no vaccination") +
  theme_bw()+
  guides(fill = FALSE, color=F, group=F) 
ggsave("Pics\\Fig3b_NPVimpact.tiff", width = 15, height = 13, units = "cm", compression="lzw", dpi =300)


# calculate min sens needed ---------------------------------------------

get_PV <- function(specificity=.9, sensitivity=.7, seroprevalence=0.7 ){
  tp= sensitivity * seroprevalence
  fn= (1-sensitivity) * seroprevalence
  tn= specificity * (1-seroprevalence)
  fp= (1-specificity) * (1-seroprevalence)
  
  PPV = tp / (tp+fp) #prop of true positive among all pos
  NPV = tn / (tn+fn) #prop of negative positive among all negative
  
  df = data.frame(specificity=specificity, sensitivity=sensitivity, seroprevalence=seroprevalence, NPV=NPV, PPV=PPV) %>%
    as.tibble()
  
  return(df)
  
}

res = NULL
for(specificity in c(.9,.95,.99)){
  for(seroprevalence in seq(0.5,1,by=0.05)){
    for(sensitivity in seq(0,1,by=0.01)){
      res = res %>%
        rbind(get_PV(specificity, sensitivity, seroprevalence))
    }
  }
}
min.PPV = .9
min.NPV = .75

res %>% group_by(specificity,seroprevalence) %>% mutate(keep_PPV = ( PPV == min(PPV + 5*(PPV<min.PPV))),
                                                        keep_NPV = ( NPV == min(NPV + 5*(NPV<min.NPV)))) %>%
  filter(keep_PPV | keep_NPV)

