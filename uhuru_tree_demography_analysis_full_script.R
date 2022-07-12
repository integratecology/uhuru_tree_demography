# PACKAGES ####
library(broom)  # converting statistical objects into tidy tibbles
library(broom.mixed)
library(car) # calculating VIFs
library(data.table)
library(lme4) # mixed effects models
library(lmerTest)
library(lubridate)
library(MASS)
library(modelr)  # calculating RMSE
library(MuMIn) # model selection
library(plotrix)
library(pROC) # ROC curve validation
library(lubridate)
library(dplyr)

conflicted::conflict_prefer('select', winner = 'dplyr')
conflicted::conflict_prefer('filter', winner = 'dplyr')
conflicted::conflict_prefer('summarize', winner = 'dplyr')
conflicted::conflict_prefer('lmer', winner = 'lmerTest')
conflicted::conflict_prefer('select', winner = 'dplyr')

# Global Options #
options(na.action = "na.fail")

# Data ####
setwd("C:\\files\\publications_academic\\in_prep\\hays_et_al_tree_defense\\final_analysis\\")

growth <- read.csv('xgrowth_dung.csv', na.strings = 'NA')
growth$elephants <- ifelse(growth$treatment=='OPEN', 1,0)
growth$impala <- ifelse(growth$treatment == 'OPEN' | growth$treatment == 'MEGA',1,0)
growth$dikdik <- ifelse(growth$treatment == 'TOTAL', 0, 1)
growth <- growth %>% 
  mutate(blockNEW = paste(site, block, sep = '')) %>% 
  select(-block) %>% 
  rename(block = blockNEW) %>% 
  filter(year != 2018) #get rid of single row w/ 2018, not sure if typo or what. 
growth$date <- lubridate::mdy(growth$date)
rain <- read.csv('xlong_rain_yearly_avg.csv', na.strings = 'NA') 
rain$date <- lubridate::ymd(rain$date)
rain$rainfall_py <- scale(rain$rainfall_py, center = F, scale = T)
# large range causes scaling issues in some regressions, particularly logistic. Remember to rescale when interpreting/visualizing!!!

growth <- left_join(growth, rain, by = c('site','date')) %>% select(-1) %>% select(-X.y) # get rid of duplicate site and weird X column introduced by excel
growth <- growth %>% filter(year<2015) # Time of data collection changed after 2014, such that there were not consistent 1 year time intervals between data collection. Can only work w/ 2009-2014 = 5 transition year


# Data Cleaning ####
# Data Prep for Growth and Survival ###

tree_survey_wide <- 
  growth %>% 
  ungroup() %>% 
  dplyr::select(year, site, treatment, rainfall_py, block, species, number, dead, elephants, impala, dikdik, ht_m, elephant_dung_avg, impala_dung_avg, dikdik_dung_avg) %>% 
  rename(ht = ht_m,  elephant_dung = elephant_dung_avg, impala_dung = impala_dung_avg, dikdik_dung = dikdik_dung_avg) %>% 
  mutate(elephants = as.factor(elephants),
         impala = as.factor(impala),
         dikdik = as.factor(dikdik)) %>% 
  drop_na(dead, number)

tree_survey_wide$surv<-ifelse(tree_survey_wide$dead=='Y',0,1)

#newer dataset has a lot of non-repeat trees, data collected only once for fertility. Need to remove those from growth and survival analyses.
tree_survey_wide$uniquenumber<- !tree_survey_wide$number %in% tree_survey_wide$number[duplicated(tree_survey_wide$number)]
tree_survey_wide<- tree_survey_wide %>% filter(uniquenumber==FALSE) %>% select(-uniquenumber, -dead) %>% filter(species != 'Acacia_drepanolobium')

#need to be able to drop NAs in height w/out removing dead trees. Create placeholder value for dead tree height, change back later
tree_survey_wide$ht[tree_survey_wide$surv==0]<- 999 #211 dead trees total

tree_survey_wide <- tree_survey_wide %>% 
  pivot_wider(names_from = year, values_from = c(ht, surv, rainfall_py, elephant_dung, impala_dung, dikdik_dung), values_fn = list(ht = mean, rainfall_py = mean, elephant_dung = mean, impala_dung = mean, dikdik_dung = mean)) #get a new column for every year for each value of height, rainfall, and dungs. In cases where there are multiple vlaues in a year, take the mean 

tree_survey_long <- 
  bind_rows(tree_survey_wide %>% 
      dplyr::select(site, treatment, species, block, number, elephants, impala, dikdik, surv_2010, ht_2009, ht_2010, rainfall_py_2010,
                    elephant_dung_2010, impala_dung_2010, dikdik_dung_2010) %>% 
      rename(ht_t = ht_2009, ht_t1 = ht_2010, surv = surv_2010, rainfall_py = rainfall_py_2010, elephant_dung_py = elephant_dung_2010,
             impala_dung_py = impala_dung_2010, dikdik_dung_py = dikdik_dung_2010) %>% 
      drop_na() %>% 
      mutate(transition = "t09-t10", 
             year = 2010), 
    
    tree_survey_wide %>% 
      dplyr::select(site, treatment, species, block, number, elephants, impala, dikdik, surv_2011, ht_2010, ht_2011, rainfall_py_2011,
                    elephant_dung_2011, impala_dung_2011, dikdik_dung_2011) %>% 
      rename(ht_t = ht_2010, ht_t1 = ht_2011,  surv = surv_2011, rainfall_py = rainfall_py_2011, elephant_dung_py = elephant_dung_2011,
             impala_dung_py = impala_dung_2011, dikdik_dung_py = dikdik_dung_2011) %>% 
      drop_na() %>% 
      mutate(transition = "t10-t11", 
             year = 2011), 
    
    tree_survey_wide %>% 
      dplyr::select(site, treatment, species, block, number, elephants, impala, dikdik, surv_2012, ht_2011, ht_2012, rainfall_py_2012,
                    elephant_dung_2012, impala_dung_2012, dikdik_dung_2012) %>% 
      rename(ht_t = ht_2011, ht_t1 = ht_2012,  surv = surv_2012, rainfall_py = rainfall_py_2012, elephant_dung_py = elephant_dung_2012,
             impala_dung_py = impala_dung_2012, dikdik_dung_py = dikdik_dung_2012) %>% 
      drop_na() %>% 
      mutate(transition = "t11-t12", 
             year = 2012),
    
    tree_survey_wide %>% 
      dplyr::select(site, treatment, species, block, number, elephants, impala, dikdik, surv_2013, ht_2012, ht_2013, rainfall_py_2013,
                    elephant_dung_2013, impala_dung_2013, dikdik_dung_2013) %>% 
      rename(ht_t = ht_2012, ht_t1 = ht_2013,  surv = surv_2013, rainfall_py = rainfall_py_2013, elephant_dung_py = elephant_dung_2013,
             impala_dung_py = impala_dung_2013, dikdik_dung_py = dikdik_dung_2013) %>% 
      drop_na() %>% 
      mutate(transition = "t12-t13", 
             year = 2013), 
    
    tree_survey_wide %>% 
      dplyr::select(site, treatment, species, block, number, elephants, impala, dikdik, surv_2014, ht_2013, ht_2014, rainfall_py_2014,
                    elephant_dung_2014, impala_dung_2014, dikdik_dung_2014) %>% 
      rename(ht_t = ht_2013, ht_t1 = ht_2014,  surv = surv_2014, rainfall_py = rainfall_py_2014, elephant_dung_py = elephant_dung_2014,
             impala_dung_py = impala_dung_2014, dikdik_dung_py = dikdik_dung_2014) %>% 
      drop_na() %>% 
      mutate(transition = "t13-t14", 
             year = 2014))  %>%  
  mutate_if(is.character, as.factor) %>% 
  mutate(block = as.factor(block), 
         elephants = as.factor(elephants), 
         impala = as.factor(impala), 
         dikdik = as.factor(dikdik)
  ) %>% 
  dplyr::select(-c(year)) %>% 
  rename(year = transition) %>%
  arrange(species)

tree_survey_long$ht_t1[tree_survey_long$surv==0]<- 0 #set ht_t1 for dead trees back to 0!

# There are 41 instances of trees growing by more than 2m in a year. That seems more likely to be an artefact of confusion around data collection than representing natural life history of trees. Remove from analyes
tree_survey_long<- tree_survey_long %>% mutate(ht_diff= ht_t1 - ht_t) %>% filter(ht_diff<=2) %>% select(-ht_diff)


# 34 trees weren't recorded for 1 or two years and then marked dead - result is that they have NA for ht_m_py making them impossible to include in regression and lowering sample sizes for mortality events. Assuming here that these trees couldn't be found and therefore should have been marked dead in the year immmediately after the last recorded height. So, if the gap doesn't start at 2018 (when no survey was conducted and therefore it's impossible to say if the tree was still alive at that point), I'm going back into data file (growth_dung) and deleting those intervening years and changing the tree to be dead in the first subsequent year. Trees changed are:
# 1005, 2081, 2423, 3265, 3282, 3283, 3284, 3285b, 3289, 749 (Only tree 1005 relevant after limiting date range to 2009-2014). CRDI 65 has 1 initial year of data, 4 years of NA, then shows up dead. Not changing and leaving in so it gets filtered out when NA in prev year are dropped b/c that's a bullshit datapoint

# tree_survey_long %>% filter(surv==0 & is.na(ht_t)==TRUE)
# tree_survey_long %>% group_by(number) %>% filter(any(surv==0 & is.na(ht_t)==TRUE))

# split out into growth and survival datasets by species

# Survival species specific subsets
sACBR <- tree_survey_long %>% filter(species == "ACBR") 
sACET <- tree_survey_long %>% filter(species == "ACET")
sACME <- tree_survey_long %>% filter(species == "ACME")
sBARO <- tree_survey_long %>% filter(species == "BARO")
sCRDI <- tree_survey_long %>% filter(species == "CRDI")

# Growth Species-specific subsets
ACBRg <- tree_survey_long %>% filter(species == "ACBR" & surv == 1) %>% select(-surv)
ACETg <- tree_survey_long %>% filter(species == "ACET" & surv == 1) %>% select(-surv)
ACMEg <- tree_survey_long %>% filter(species == "ACME" & surv == 1) %>% select(-surv)
BAROg <- tree_survey_long %>% filter(species == "BARO" & surv == 1) %>% select(-surv)
CRDIg <- tree_survey_long %>% filter(species == "CRDI" & surv == 1) %>% select(-surv)


# ACBR ####
# Check number of deaths
sACBR %>% group_by(surv, treatment) %>% tally()

# Check sample size 
nrow(sACBR) > 10 * (5 / (nrow(sACBR[sACBR$surv==0,]) / nrow(sACBR[sACBR$surv==1,]) ) ) # Low sample size
nrow(sACBR %>% filter(surv == 0)) / 5 #max 4 predictors [max predictos  = smaller number of success/failure events / 5] 
# !!!! but going to forgo constraints on predictors and simply accept any models that converge

# Fit full model
full_model_ACBR <-  glmer(surv ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), family = binomial(link = 'logit'), data = sACBR)
summary(full_model_ACBR)

# Calculate VIFs 
vif(full_model_ACBR) # All < 4

# Model selection with AICc
options(na.action = "na.fail")
m_list_ACBRs <- dredge(full_model_ACBR, rank="AICc")
m_list_ACBRs[1:10,] # 6 top models, no changes in sign, some singular fitting but global model and top model are good. No changes in signs of parameters

ACBRs.best<- get.models(m_list_ACBRs, 1)[[1]]  #AML CHANGED

# Model Validation with AUC
ROC_model_ACBR <- roc(sACBR$surv, predict(full_model_ACBR, backTransform=TRUE), plot=T)
ROC_model_ACBR # predictive power = not good, looks somewhat overfit

r.squaredGLMM(ACBRs.best) #R2m = 0.09 R2C = 0.11
summary(ACBRs.best)


# ACET ####
# Check number of deaths
sACET %>% group_by(surv, treatment) %>% tally()

# Check Sample size
nrow(sACET) > 10 * (5 / (nrow(sACET[sACET$surv==0,]) / nrow(sACET[sACET$surv==1,]) ) ) #low sample size
nrow(sACET %>% filter(surv == 0)) / 5 #max 4 predictors

# Fit full model
full_model_ACET <- glmer(surv ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), family = binomial(link = 'logit'), data = sACET)

# Check VIFs
vif(full_model_ACET) # All < 4

# Model selection with AICc
options(na.action = "na.fail")
m_list_ACETs <- dredge(full_model_ACET, rank="AICc")
m_list_ACETs[1:10,] # 3 top models, several singularly fit models but global/top aren't 

ACETs.best<-  get.models(m_list_ACETs, 1)[[1]] # AML CHANGED

# Model Validation 
ROC_model_ACET <- roc(sACET$surv, predict(full_model_ACET, backTransform=TRUE), plot=T)
ROC_model_ACET # predictive power = good, not overfit

r.squaredGLMM(ACETs.best) #R2m = 0.79, R2c = 0.82
summary(ACETs.best)


# ACME ####
# Check number of deaths
sACME %>% group_by(surv,  treatment) %>% tally()

# Check sample size
nrow(sACME) > 10 * (5 / (nrow(sACME[sACME$surv==0,]) / nrow(sACME[sACME$surv==1,]) ) ) #ok
nrow(sACME %>% filter(surv == 0)) / 5 #max 14 predictors

# Fit full model
full_model_ACME <-  glmer(surv ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), family = binomial(link = "logit"), data = sACME)
# Calculate VIFs 
vif(full_model_ACME) # All < 4

# Model selection with AICc
options(na.action = "na.fail")
m_list_ACMEs <- dredge(full_model_ACME, rank="AICc")
m_list_ACMEs[1:10,] # 5 top models, no change in sign. 

ACMEs.best<- get.models(m_list_ACMEs, 1)[[1]] # AML CHANGED

# Model Validation with AUC
ROC_model_ACME <- roc(sACME$surv, predict(full_model_ACME, backTransform=TRUE), plot=T)
ROC_model_ACME # predictive power =  good, not overfit

r.squaredGLMM(ACMEs.best) # R2m = 0.15, R2c = 0.22
summary(ACMEs.best)


# BARO ####
# Check number of deaths
sBARO %>% group_by(surv, site, treatment) %>% tally()
# only 4 deaths across all treatments -- cannot model death rate (2 deaths not shown here b/c have NA in ht_t)

sBARO %>% group_by(surv, year) %>% tally()
length(unique(sBARO$number))
# instead, going to assume a constant mortality rate - 4 deaths in 11 years out of a total of 116 trees is 1-((4/116)/11) = 0.9968652 (i.e. a 99.7% chance of survival annually)

nrow(sBARO) > 10 * (5 / (nrow(sBARO[sBARO$surv==0,]) / nrow(sBARO[sBARO$surv==1,]) ) )


# CRDI ####
# Check number of deaths
sCRDI %>% group_by(surv, treatment) %>% tally()

# Check sample size
nrow(sCRDI) > 10 * (5 / (nrow(sCRDI[sCRDI$surv==0,]) / nrow(sCRDI[sCRDI$surv==1,]) ) ) #low sample size
nrow(sCRDI %>% filter(surv == 0)) / 5 #4.2, so max 4 predictors

# Fit full model
full_model_CRDI <-  glmer(surv ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), family = binomial(link = "logit"), data = sCRDI)

# Calculate VIFs 
vif(full_model_CRDI) # All < 4.0

# Model selection with AICc
options(na.action = "na.fail")
m_list_CRDIs <- dredge(full_model_CRDI, rank="AICc")
m_list_CRDIs[1:10,] # 4 top models no change in sign. some models fail to converge, global/top are fine

CRDIs.best<- get.models(m_list_CRDIs,1)[[1]] # AML CHANGED

# Model Validation with AUC
ROC_model_CRDI <- roc(sCRDI$surv, predict(full_model_CRDI, backTransform=TRUE), plot=T)
ROC_model_CRDI # predictive power = pretty good, not overfit

r.squaredGLMM(CRDIs.best) #  R2m = 0.27, R2c = 0.35
summary(CRDIs.best)


# GROWTH REGRESSIONS ####
# ACBR Growth ####
nrow(ACBRg)

# Fit global model
ACBR.fit.global <-  lmer(ht_t1 ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), REML = F, data = ACBRg)

# Calculate VIFs (don't include interactions--they're collinear by definition)
vif(ACBR.fit.global) # All < 4

# Selection
m_list_ACBRg <- dredge(ACBR.fit.global, rank="AICc")
head(m_list_ACBRg, 10) # 2 top models

ACBRg.best <- update(get.models(m_list_ACBRg, 1)[[1]], REML=T)  #AML CHANGED

# Summary of best model
r.squaredGLMM(ACBRg.best) #R2m = 0.71, R2c = 0.71
summary(ACBRg.best)


# ACET Growth ####
nrow(ACETg)

# Fit global model
ACET.fit.global <-  lmer(ht_t1 ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), REML = F, data = ACETg)

# Calculate VIFs (don't include interactions--they're collinear by definition)
vif(ACET.fit.global) # all < 4

# Selection
m_list_ACETg <- dredge(ACET.fit.global, rank="AICc")
head(m_list_ACETg, 10) # 3 top models

ACETg.best <- update(get.models(m_list_ACETg, 1)[[1]], REML=T) #AML CHANGED

r.squaredGLMM(ACETg.best) #R2m = 0.80, R2c = 0.81
summary(ACETg.best)


# ACME Growth ####
nrow(ACMEg)

# Fit global model
ACME.fit.global <-  lmer(ht_t1 ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), REML = F, data = ACMEg)

# Calculate VIFs (don't include interactions--they're collinear by definition)
vif(ACME.fit.global) # all < 4

# Selection
m_list_ACMEg <- dredge(ACME.fit.global, rank="AICc")
head(m_list_ACMEg, 10) # 4 top models, dikdik changes in sign not in top models, several singularly fit models but not best one

ACMEg.best <-  update(get.models(m_list_ACMEg, 1)[[1]], REML=T)  #AML CHANGED
r.squaredGLMM(ACMEg.best) # R2m = 0.80, R2c = 0.80
summary(ACMEg.best)


# BARO Growth ####
nrow(BAROg)

# Fit global model
BARO.fit.global <-  glm(ht_t1 ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + block, data = BAROg) #global model singularly fit w/ random effect of block, runs with fixed effect

# Calculate VIFs (don't include interactions--they're collinear by definition)
vif(BARO.fit.global) # all < 4

# Selection
m_list_BAROg <- dredge(BARO.fit.global, rank="AICc")
head(m_list_BAROg, 10) # 6 top models, impala changes sign but not in top models, no singular fitting, block not included in any of the top models

BAROg.best <- get.models(m_list_BAROg, 1)[[1]] #AML CHANGED -- NB no singular effect

r.squaredGLMM(BAROg.best) #R2m -= 0.86, R2c = 0.86
summary(BAROg.best)


# CRDI Growth ####
nrow(CRDIg)

# Fit global model
CRDI.fit.global <-  lmer(ht_t1 ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), REML = F, data = CRDIg)

# Calculate VIFs (don't include interactions--they're collinear by definition)
vif(CRDI.fit.global) # all < 4

# Selection
m_list_CRDIg <- dredge(CRDI.fit.global, rank="AICc")
head(m_list_CRDIg, 10) # 4 top models,  no change in sign, some singularly fit but not top model

CRDIg.best <- update(get.models(m_list_CRDIg, 1)[[1]], REML=T)
r.squaredGLMM(CRDIg.best) #R2m = 0.71, R2c = 0.71
summary(CRDIg.best)


# GROWTH VARIANCE
#r Data Cleaning III########
## Bryan's old code to Calculate residuals for growth variance regressions
# (ht_regs <- 
#   tibble(species = c("ACBR","ACET","ACME","BARO","CRDI"), 
#          models = list(ACBRg.best, 
#                        ACETg.best, 
#                        ACMEg.best, 
#                        BAROg.best,
#                        CRDIg.best
#                        )) %>% 
#   mutate(tidy = map(models, broom.mixed::tidy), 
#          aug = map(models, broom.mixed::augment)) 
# )
# 
# ht_regs %>%
#   transmute(species, tidy) %>%
#   unnest(cols = c(tidy)) %>% 
#   print(n = Inf)
# 
# # Calculate residuals
# ht_resids <- 
#   ht_regs %>%
#   transmute(species, aug) %>%
#   unnest(cols = c(aug)) %>% 
#   select(species, ht_t, .resid) %>% 
#   rename(residuals = .resid) 
#   bind_cols(tree_survey_long[,c("treatment","rainfall_py","block", "year","elephants","impala","dikdik", "elephant_dung_py", "impala_dung_py", "dikdik_dung_py")])
# 
# ht_resids %>% print(n = 20

#default behavior of lme4 package's override for the residuals function w/ lme is 'response' (i.e. observed - fitted), same as base residuals function treatment of glm objects. Thus no problem using same function for glm and lme models (though lme4 override does have different default behavior for glmer models)
ACBRvar<- ACBRg %>% mutate(residuals = resid(ACBRg.best)^2)
ACETvar<- ACETg %>% mutate(residuals = resid(ACETg.best)^2)
ACMEvar<- ACMEg %>% mutate(residuals = resid(ACMEg.best)^2)
BAROvar<- BAROg %>% mutate(residuals = resid(BAROg.best)^2)
CRDIvar<- CRDIg %>% mutate(residuals = resid(CRDIg.best)^2)

# Check residuals' distribution
graphvar<- rbind(ACBRvar, ACETvar, ACMEvar, BAROvar, CRDIvar)
graphvar %>%
  ggplot(aes(residuals)) +
  geom_histogram(color = "black", fill = "grey")+
  facet_grid(species ~ treatment)


# ACBR Variance in Growth ####
# Fit global model
ACBR.resid.global <-  lmer(residuals ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), REML = F, data = ACBRvar)

# Calculate VIFs (don't include interactions--they're collinear by definition)
vif(ACBR.resid.global) # all < 4

# Selection
m_list_ACBRv <- dredge(ACBR.resid.global, rank="AICc")
head(m_list_ACBRv, 10) #1 top supported models, no change in signs, spotty variable inclusion = poor selection, top model is intercept only

ACBRres.best<- update(get.models(m_list_ACBRv,1)[[1]], REML=T) # AML CHANGED

r.squaredGLMM(ACBRres.best) #R2m = 0.025 R2c = 0.035
summary(ACBRres.best)


# ACET Variance in Growth####
# Fit global model
ACET.resid.global <-   glm(residuals ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + block, data = ACETvar) #singular fit w/ random effects

# Calculate VIFs (don't include interactions--they're collinear by definition)
vif(ACET.resid.global) # all < 4

# Selection
m_list_ACETv <- dredge(ACET.resid.global, rank="AICc")
head(m_list_ACETv, 10) # 5 top models, again intercept only model is best, elephant changes sign though not in top models

ACETres.best<- get.models(m_list_ACETv,1)[[1]]# AML CHANGED

r.squaredGLMM(ACETres.best) #R2m = 0.015  R2c = 0.015
summary(ACETres.best)


# ACME Variance in Growth ####
# Fit global model
ACME.resid.global <-   glm(residuals ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + block, data = ACMEvar) #singular fit w/ random effects

#Calculate VIF
vif(ACME.resid.global) # all < 4

# Selection
m_list_ACMEv <- dredge(ACME.resid.global, rank="AICc")
head(m_list_ACMEv, 10) #3 top models, no changes in sign,

ACMEres.best<- get.models(m_list_ACMEv,1)[[1]]# AML CHANGED

r.squaredGLMM(ACMEres.best)# R2m = 0.1  R2c = 0.1
summary(ACMEres.best)


# BARO Variance in Growth ####
# Fit global model
BARO.resid.global <- lmer(residuals ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py +  dikdik_dung_py + (1|block), REML = F, data = BAROvar) #singular fit w/ rnadom effects

# Calculate VIFs
vif(BARO.resid.global) # all < 4

# Selection
m_list_BAROv <- dredge(BARO.resid.global, rank="AICc")
head(m_list_BAROv, 10) #5 top models, impala changes sign but not in top model, some models singularly fit but top and best are fine

BAROres.best<- update(get.models(m_list_BAROv,1)[[1]], REML=T) # AML CHANGED

r.squaredGLMM(BAROres.best)# R2m = 0.03   R2c = 0.05
summary(BAROres.best)


# CRDI Variance in Growth ####
# Fit global model
CRDI.resid.global <-   lmer(residuals ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), REML = F, data = CRDIvar) #singular fit w/ random effects

# Calculate VIFs (don't include interactions--they're collinear by definition)
vif(CRDI.resid.global) # all < 4

# Selection
m_list_CRDIv <- dredge(CRDI.resid.global, rank="AICc")
head(m_list_CRDIv, 10) # 4 top models,no change in signs

CRDIres.best<- update(get.models(m_list_CRDIv,1)[[1]], REML=T) # AML CHANGED

r.squaredGLMM(CRDIres.best) #R2m = 0.01, R2c = 0.01
summary(CRDIres.best)


# FERTILITY ####
# Data Cleaning IV ####

#pull initial data from overall database
fert <- growth %>% 
  filter(year == 2013 | year == 2014) %>% filter (dead =='N' & species != 'Acacia_drepanolobium') %>% 
  ungroup() %>% 
  dplyr::select(-c(survey, plot:date, old_number:dead,  l_m:no_stems, new_old:interval)) %>% 
  rename(ht = ht_m,  elephant_dung = elephant_dung_avg, impala_dung = impala_dung_avg, dikdik_dung = dikdik_dung_avg) %>% 
  mutate(elephants = as.factor(elephants),
         impala = as.factor(impala),
         dikdik = as.factor(dikdik), 
  ) %>% 
  drop_na(number)

# widen so can get 2013-2014 transition data
fert_wide<-  fert %>% 
  pivot_wider(names_from = year, values_from = c(ht, no_flowers, no_buds, no_fruits, rainfall_py, elephant_dung, impala_dung, dikdik_dung), values_fn = list(ht = mean, rainfall_py = mean, elephant_dung = mean, impala_dung = mean, dikdik_dung = mean)) #get a new column for every year for each value of height, rainfall, and dungs. In cases where there are multiple vlaues in a year, take the mean (though there shouldn't be any after accounting for repeat trees in first chunk)

# Deciding whether to use 2013-2014 trasniion or just 2014 data
nrow(fert_wide%>% filter(ht_2013>0 & ht_2014>0) %>% filter(no_fruits_2014 >0 | no_buds_2014 >0 | no_flowers_2014 > 0)) #1329 trees with data for both 2013 and 2014, of which 736 are reproductive and 716 have fruit
nrow(fert_wide%>% filter(ht_2014>0 & no_fruits_2014)) #2928 trees w/ data in 2014, of which 1108 have fruits
#enough data to use transition year, much more rigorous than just using 2014

# recapitulate to isolate transition year
fert_long <- fert_wide %>% 
  dplyr::select(site:ht_2014, no_flowers_2014, no_buds_2014, no_fruits_2014, rainfall_py_2014, elephant_dung_2014,
                impala_dung_2014, dikdik_dung_2014) %>% 
  rename(ht_t = ht_2013, ht_t1 = ht_2014, flowers = no_flowers_2014,  buds = no_buds_2014, fruits = no_fruits_2014,
         rainfall_py = rainfall_py_2014, elephant_dung_py = elephant_dung_2014,impala_dung_py = impala_dung_2014,
         dikdik_dung_py = dikdik_dung_2014) %>% 
  mutate(repro = ifelse(fruits >0 | buds >0 | flowers > 0, 1, 0)) %>% 
  drop_na(fruits) #leaves in one entry w/ na for flowers but data for fruits

fert_long2<- fert_long %>% drop_na(ht_t1)#just 2014 data
fert_long<-fert_long %>% drop_na(ht_t) #2013-2014 transition data

#View(fert_long%>% filter(site == 'N') %>%  group_by(species, treatment, repro) %>%  summarize(n()))  
#have 5 or fewer trees w fruits for ACME in S open, ACBR C/S/N open, and BARO C/S open + all of north. Otherwise pretty good ample size
#View(fert_long2%>% filter(site =='N') %>%  group_by(species, treatment, repro) %>%  summarize(n()))
#using 2014 data mostly increases sample sizes where sample sizes are already adequate, doesn't help low sample sizes for BARO, still no reproductive trees in ACBR N open. Think it's more ecologically defensible to use height in 2013 as predictor of fruiting in 2014, will model that way unless model reuslts are bad.

nrow(fert_long %>% filter(ht_t1 <=1 & repro ==1)) #only one reroductive tree that's 1m tall or less (only has flowers too), justified in using as sapling class

# 2 trees w/ height change >=2 in 2013-2014 transition, filter them out but keep single tree w/ ht_2014==NA b/c it has fruits and 2014 height isn't included in regression
fert_long<- fert_long %>% mutate(ht_diff= ht_t1 - ht_t) %>% filter(ht_diff<=2|is.na(ht_diff)==TRUE) %>% select(-ht_diff, -ht_t1)

nrow(fert_long %>% filter(flowers>0 & fruits<1)) # 13 trees w/ flowers but no fruit. 
nrow(fert_long %>% filter(buds >0 & fruits<1)) # 12 trees w/ buds but no fruit.  25 out of total of 1330, not a problem. For consistency's sake filter them out and only keep trees w/ fruits >1
fert_long<- fert_long %>% mutate(repro = ifelse(fruits>0,1,0))

# Determined that using transformation of fruit number data is best for seed production and produces mostly normally distributed data, bit of a left skew. There are 20 instances of reproductive trees with # fruits < height, resulting in negative log values. So added 1 before taking log. Preserves same distribution - just increases left skew a bit
fert_long<- fert_long %>% mutate(seed = log((fruits / ht_t)+1))
plot(seed ~ ht_t, data = fert_long)  #horizontal lines in data due to estimation of fruit #? All data are in increments of 5/10, not problematic but worth noting
hist(fert_long$seed[fert_long$repro==1])

# create separate dataset by species
ACBRfert<- fert_long %>% filter(species == 'ACBR')
ACETfert<- fert_long %>% filter(species == 'ACET')
ACMEfert<- fert_long %>% filter(species == 'ACME')
BAROfert<- fert_long %>% filter(species == 'BARO') #BARO has no trees w/ data in N-Mega
CRDIfert<- fert_long %>% filter(species == 'CRDI')

ACBRseed<- ACBRfert %>% filter(repro == 1)
ACETseed<- ACETfert %>% filter(repro == 1)
ACMEseed<- ACMEfert %>% filter(repro == 1)
BAROseed<- BAROfert %>% filter(repro == 1)
CRDIseed<- CRDIfert %>% filter(repro == 1)

# look at height distributions of trees across treatment/site combos
dummy<-data.frame( BAROfert %>% group_by(site, treatment) %>% summarize(avg = mean(ht_t), CIup = mean(ht_t) + 1.96*(sd(ht_t)/sqrt(n())), CIlow = mean(ht_t) - 1.96*(sd(ht_t)/sqrt(n())) ) )
ggplot(BAROfert, aes(x=ht_t)) + 
  geom_histogram() + 
  facet_grid( rows = vars(site), cols = vars(treatment))+
  geom_vline(data = dummy, aes(xintercept = avg), 
             color = 'red') +
  geom_vline(data = dummy, aes(xintercept = CIup), 
             color = 'black')+
  geom_vline(data = dummy, aes(xintercept = CIlow), 
             color = 'black')
# only small reproductive trees (<1m) are in South, hard to say, but looks like central distributions are further to the right/taller. Almost no reproductive trees in the north



# ACBR Probability of Reproduction ####
# Check number of reproductive trees
ACBRfert %>% group_by(repro, site, treatment) %>% tally()

# Check sample size
nrow(ACBRfert) > 10 * (5 / (nrow(ACBRfert[ACBRfert$repro==1,]) / nrow(ACBRfert[ACBRfert$repro==0,]) ) ) #ok
nrow(ACBRfert[ACBRfert$repro==1,])/5 #26 max predictors

# Fit full model
full_model_ACBRpr <-  glmer(repro ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), family = binomial(link = "logit"), data = ACBRfert)
summary(full_model_ACBRpr) 

# Calculate VIFs 
vif(full_model_ACBRpr) # all < 4

# Model selection with AICc
m_list_ACBRpr <- dredge(full_model_ACBRpr, rank="AICc")
m_list_ACBRpr[1:10,] # 4 top models, no change in signs for covariates, change in intercept not a problem. 

ACBRpr.best<-get.models(m_list_ACBRpr,1)[[1]] # AML CHANGED
summary(ACBRpr.best)
r.squaredGLMM(ACBRpr.best) #R2m = 0.32, R2c = 0.40

# Model Validation with AUC
ROC_model_ACBRpr <- roc(ACBRfert$repro, predict(full_model_ACBRpr, backTransform=TRUE), plot=T)
ROC_model_ACBRpr # predictive power = pretty good, doesn't look overfit


# ACET Probability of Reproduction ####
#Check number of reproductive trees
ACETfert %>% group_by(repro, site, treatment) %>% tally()

# Check sample size
nrow(ACETfert) > 10 * (5 / (nrow(ACETfert[ACETfert$repro==1,]) / nrow(ACETfert[ACETfert$repro==0,]) ) ) #ok
nrow(ACETfert[ACETfert$repro==1,])/5 #29 max predictors

# Fit full model
full_model_ACETpr <-  glmer(repro ~ ht_t + rainfall_py  + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), family = binomial(link = "logit"), data = ACETfert)
summary(full_model_ACETpr) 

# Calculate VIFs 
vif(full_model_ACETpr) # all < 4

# Model selection with AICc
m_list_ACETpr <- dredge(full_model_ACETpr, rank="AICc")
m_list_ACETpr[1:10,] # 5 top models, no changes in sign. 

ACETpr.best <- get.models(m_list_ACETpr,1)[[1]] # AML CHANGED
summary(ACETpr.best)
r.squaredGLMM(ACETpr.best) #R2m = 0.42, R2c = 0.44

# Model Validation with AUC
ROC_model_ACETpr <- roc(ACETfert$repro, predict(full_model_ACETpr, backTransform=TRUE), plot=T)
ROC_model_ACETpr # predictive power = good, doesn't look overfit


#r ACME Probability of Reproduction ####
#Check number of reproductive trees
ACMEfert %>% group_by(repro, site, treatment) %>% tally()

# Check sample size
nrow(ACMEfert) > 10 * (5 / (nrow(ACMEfert[ACMEfert$repro==1,]) / nrow(ACMEfert[ACMEfert$repro==0,]) ) ) #ok
nrow(ACMEfert[ACMEfert$repro==1,])/5 #30 max predictors

# Fit full model
full_model_ACMEpr <-  glmer(repro ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), family = binomial(link = "logit"), data = ACMEfert)
summary(full_model_ACMEpr) 

# Calculate VIFs 
vif(full_model_ACMEpr) # all < 4

# Model selection with AICc
m_list_ACMEpr <- dredge(full_model_ACMEpr, rank="AICc")
m_list_ACMEpr[1:10,] # 3 top models, no changes in sign for covariates, change in sign for intercept but not in top models. 

ACMEpr.best <- get.models(m_list_ACMEpr,1)[[1]] # AML CHANGED

summary(ACMEpr.best)
r.squaredGLMM(ACMEpr.best) #R2m = 0.35, R2c = 0.43

# Model Validation with AUC
ROC_model_ACMEpr <- roc(ACMEfert$repro, predict(full_model_ACMEpr, backTransform=TRUE), plot=T)
ROC_model_ACMEpr # predictive power = good, doesn't look overfit


# BARO Probability of Reproduction ####
# Check number of reproductive trees
BAROfert %>% group_by(repro, site, treatment) %>% tally()

# Check sample size
nrow(BAROfert) > 10 * (5 / (nrow(BAROfert[BAROfert$repro==0,]) / nrow(BAROfert[BAROfert$repro==1,]) ) ) #low sample size, actually have more reproductive than non reproductive trees
nrow(BAROfert[BAROfert$repro==0,])/5 #4 max predictors

# Fit full model
full_model_BAROpr <-  glm(repro ~ ht_t + block +impala_dung_py + dikdik_dung_py, family = binomial(link = "logit"), data = BAROfert) #singular effect using block as random effect, if include as fixed effect highly correlated w/ rainfall and elephant, leaving elephant and rainfall out per Allisons suggestion
summary(full_model_BAROpr) 

# Calculate VIFs 
vif(full_model_BAROpr) # all < 4

# Model selection with AICc
m_list_BAROpr <- dredge(full_model_BAROpr, rank="AICc", m.lim = c(0,4))
m_list_BAROpr[1:10,] # 3 top models, no changes in sign

BAROpr.best<- get.models(m_list_BAROpr,1)[[1]] # AML CHANGED
summary(BAROpr.best)
r.squaredGLMM(BAROpr.best) #R2m = 0.42, R2c = 0.42 

# Model Validation with AUC
ROC_model_BAROpr <- roc(BAROfert$repro, predict(full_model_BAROpr, backTransform=TRUE), plot=T)
ROC_model_BAROpr # predictive power =  good, doesn't look overfit


# CRDI Probability of Reproduction ####
# Check number of reproductive trees
CRDIfert %>% group_by(repro, site, treatment) %>% tally()

# Check sample size
nrow(CRDIfert) > 10 * (5 / (nrow(CRDIfert[CRDIfert$repro==0,]) / nrow(CRDIfert[CRDIfert$repro==1,]) ) ) #ok, here also have more reproductive than non reproductive trees
nrow(CRDIfert[CRDIfert$repro==0,])/5 #11 max predictors

# Fit full model
full_model_CRDIpr <-  glmer(repro ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), family = binomial(link = "logit"), data = CRDIfert)

# Calculate VIFs 
vif(full_model_CRDIpr) # all < 4

# Model selection with AICc
m_list_CRDIpr <- dredge(full_model_CRDIpr, rank="AICc")
m_list_CRDIpr[1:10,] # 3 top models, changes in sign for elephant_dung and rainfall, but not in top models

CRDIpr.best <- get.models(m_list_CRDIpr,1)[[1]] # AML CHANGED

summary(CRDIpr.best)
r.squaredGLMM(CRDIpr.best) # R2m = 0.07, R2c = 0.24

# Model Validation with AUC
ROC_model_CRDIpr <- roc(CRDIfert$repro, predict(full_model_CRDIpr, backTransform=TRUE), plot=T)
ROC_model_CRDIpr # predictive power =  decent, doesn't look overfit


# ACBR Seed Production ####
# Check sample sizes
ACBRseed %>% group_by(site, treatment) %>% tally() 

# Fit full model
full_model_ACBRsp <- lmer(seed ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), REML = F, data = ACBRseed)

# Calculate VIFs 
vif(full_model_ACBRsp) # all < 4

# Model selection with AICc
m_list_ACBRsp <- dredge(full_model_ACBRsp, rank="AICc")
m_list_ACBRsp[1:10,] # 7 top models, no change sin sign for covariates, change in intercept sign

ACBRsp.best<- update(get.models(m_list_ACBRsp,1)[[1]], REML=T)# AML CHANGED
summary(ACBRsp.best)
r.squaredGLMM(ACBRsp.best) #R2m = 0.05, R2c = 0.18


# ACET Seed Production ####
# Check sample sizes
ACETseed %>% group_by(site, treatment) %>% tally() 

# Fit full model
full_model_ACETsp <- lmer(seed ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), REML = F, data = ACETseed)

# Calculate VIFs 
vif(full_model_ACETsp) # all < 4

# Model selection with AICc
m_list_ACETsp <- dredge(full_model_ACETsp, rank="AICc")
m_list_ACETsp[1:10,] # 2 top models, some models produce singular fitting, but global and top models are fine.

ACETsp.best <- update(get.models(m_list_ACETsp,1)[[1]], REML=T)# AML CHANGED

summary(ACETsp.best)
r.squaredGLMM(ACETsp.best) #R2m = 0.23, R2c = 0.26


# ACME Seed Production ####
#Check sample sizes
ACMEseed %>% group_by(site, treatment) %>% tally() 

# Fit full model
full_model_ACMEsp <- glm(seed ~ ht_t + elephant_dung_py + impala_dung_py + dikdik_dung_py + block, data = ACMEseed) #singular fit w/ random effects, if block included as fixed effect it's collinear w/ rainfall, removed rainfall consistent w/ Allison's suggestion

# Calculate VIFs 
vif(full_model_ACMEsp) # all < 4

# Model selection with AICc
m_list_ACMEsp <- dredge(full_model_ACMEsp, rank="AICc")
m_list_ACMEsp[1:10,] # 3 top models,

ACMEsp.best <- get.models(m_list_ACMEsp,1)[[1]]# AML CHANGED

summary(ACMEsp.best)
r.squaredGLMM(ACMEsp.best) #R2 = 0.1


# BARO Seed Production ####
# Check sample sizes
BAROseed %>% group_by(site, treatment) %>% tally() 

# Fit full model
full_model_BAROsp <- lmer(seed ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), REML = F, data = BAROseed)

# Calculate VIFs 
vif(full_model_BAROsp) # all < 4

# Model selection with AICc
m_list_BAROsp <- dredge(full_model_BAROsp, rank="AICc")
m_list_BAROsp[1:10,] # 5 top models, no changes in sign

BAROsp.best <- get.models(m_list_BAROsp,1)[[1]]# AML CHANGED

summary(BAROsp.best)
r.squaredGLMM(BAROsp.best) #R2m = 0.15, R2c = 0.26


# CRDI Seed Production ####
CRDIseed %>% group_by(site, treatment) %>% tally() 

# Fit full model
full_model_CRDIsp <- lmer(seed ~ ht_t + rainfall_py + elephant_dung_py + impala_dung_py + dikdik_dung_py + (1|block), REML = F, data = CRDIseed)

# Calculate VIFs 
vif(full_model_CRDIsp) # all < 4

# Model selection with AICc
m_list_CRDIsp <- dredge(full_model_CRDIsp, rank="AICc")
m_list_CRDIsp[1:11,] # 10 top models, no changes in sign. Some models produce singular fitting but global model and top model are fine.

CRDIsp.best <- get.models(m_list_CRDIsp,1)[[1]]# AML CHANGED

summary(CRDIsp.best)
r.squaredGLMM(CRDIsp.best) # R2m = 0.09, R2c = 0.1

# MATRICES ####
# setting up input data to use for predictions
vec_treatment<- c('TOTAL', 'MEGA', 'MESO', 'OPEN')
vec_year<- c(2010, 2013, 2014, 2011, 2012)
vec_year_repro<- c(2010, 2013, 2014)
vec_species<- c('ACBR', 'ACET', 'ACME', 'BARO', 'CRDI')
vec_block<- c("C1","C2","C3","N1","N2","N3","S1","S2","S3")
# ACBRg/s and growth have filtering applied. Use sACBR b/c don't want to extrapolate to heights way outside range of vital rate regression inputs. Numbers in comments were obtained from max(ht_t) in growth dataframe prior to filtering year<2015 (most max values occ ur in year==2020)
# max for CRDI in growth is actually 9.58, but that's almost double the height of next highest tree. Seems like an outlier so included 2nd tallest at 5.0
height_ACBR<-max(max(sACBR$ht_t), max(sACBR$ht_t1)) #8.1
height_ACET<-max(max(sACET$ht_t), max(sACET$ht_t1)) #7.5
height_ACME<-max(max(sACME$ht_t), max(sACME$ht_t1)) #10
height_CRDI<-max(max(sCRDI$ht_t), max(sCRDI$ht_t1)) #5.0 / 9.58
height_BARO<-max(max(sBARO$ht_t), max(sBARO$ht_t1)) #8

heights<-list(height_ACBR, height_ACET, height_ACME, height_CRDI, height_BARO)
names(heights)<- c('ACBR', 'ACET', 'ACME', 'CRDI', 'BARO')

binedges_ACBR <- c(1, 2, 3, 4, height_ACBR)
binedges_ACET <- c(1, 2, 3, 4, height_ACET)
binedges_ACME <- c(1, 2, 3, 4, height_ACME)
binedges_BARO <- c(1, 2, 3, 4, height_BARO)
binedges_CRDI <- c(1, 2, 3, height_CRDI) # height_CRDI = 3.92, don't include bin edge 4

# binedges_ACBR <- seq(from = 1, to = plyr::round_any(height_ACBR, 0.1, ceiling), by = 0.1) # AML CHANGED
binmids_ACBR <- (binedges_ACBR[1:(length(binedges_ACBR)-1)] + binedges_ACBR[2:(length(binedges_ACBR))])/2# AML CHANGED
# binedges_ACET <- seq(from = 1, to = plyr::round_any(height_ACET, 0.1, ceiling), by = 0.1)# AML CHANGED
binmids_ACET <- (binedges_ACET[1:(length(binedges_ACET)-1)] + binedges_ACET[2:(length(binedges_ACET))])/2# AML CHANGED
# binedges_ACME <- seq(from = 1, to = plyr::round_any(height_ACME, 0.1, ceiling), by = 0.1)# AML CHANGED
binmids_ACME <- (binedges_ACME[1:(length(binedges_ACME)-1)] + binedges_ACME[2:(length(binedges_ACME))])/2# AML CHANGED
# binedges_BARO <- seq(from = 1, to = plyr::round_any(height_BARO, 0.1, ceiling), by = 0.1)# AML CHANGED
binmids_BARO <- (binedges_BARO[1:(length(binedges_BARO)-1)] + binedges_BARO[2:(length(binedges_BARO))])/2# AML CHANGED
# binedges_CRDI <- seq(from = 1, to = plyr::round_any(height_CRDI, 0.1, ceiling), by = 0.1)# AML CHANGED
binmids_CRDI <- (binedges_CRDI[1:(length(binedges_CRDI)-1)] + binedges_CRDI[2:(length(binedges_CRDI))])/2# AML CHANGED

all.binedges <- list(binedges_ACBR, binedges_ACET, binedges_ACME, binedges_BARO, binedges_CRDI)
all.binmids <- list(binmids_ACBR, binmids_ACET, binmids_ACME, binmids_BARO, binmids_CRDI)
names(all.binedges) <- names(all.binmids)  <- vec_species

binedges_sap <- c(0,1) # AML CHANGED
binmids_sap <- 0.5 # NB: this is a methodological change-- I am just going to assume that all saplings survive at a similar rate. Brandon's way of doing this was unnecesarily complex, imo

# getting #'s of recruits and adults 



# Number of recruits at time t+1 #
# First step, read in census dates, lengthen, as.date, replace mixed treatment vocab to be consistent w/ rest of analyses
cendat<- data.table::fread('tree_census_dates_bestguess.csv', header = T)
cendat_long<- reshape2::melt(cendat, id.vars = 1, measure.vars = 2:length(names(cendat)), variable.name = 'year', variable.factor = F, value.name = 'CensusDate')
cendat_long$year<- as.numeric(as.character(cendat_long$year))
cendat_long$CensusDate<-lubridate::dmy(cendat_long$CensusDate)
cendat_long$plot<- gsub(pattern = 'CTL', replacement = 'OPEN', cendat_long$plot)
cendat_long$plot<- gsub(pattern = 'LMH', replacement = 'TOTAL', cendat_long$plot)

# Second step, read in census data, filter to species of interest, merge w/ dates
census<- read.csv('xtree_census_2009-2020_detailed_lab_paper_dec20.csv') %>% 
  mutate(blockNEW = paste(site, block, sep = '')) %>% 
  select(-block) %>% 
  rename(block = blockNEW, sap = X.0.5M, c1 = X0.5_1M, combosap = X.1M, c2 = X1_2M, c3 = X2_3M, c4 = X3_4M, c5 = X.4M) 

census$species<- gsub('Acacia_brevispica', 'ACBR', census$species)
census$species<- gsub('Acacia_etbaica', 'ACET', census$species)
census$species<- gsub('Acacia_mellifera', 'ACME', census$species)
census$species<- gsub('Balanites_rotundifolia', 'BARO', census$species)
census$species<- gsub('Croton_dichogamus', 'CRDI', census$species)
                                     
census2<-census %>% 
  filter(species %in% c('ACBR','ACET','ACME','BARO','CRDI')) %>% 
  select(-c(census, sap, c1, c2, c3, c4, c5, total))

census2<- merge(census2, cendat_long, by = c('plot', 'year'), all.y = T)

# Third step, limit to only trees <1m tall, sum across subplots, calculate within block change in saplings and intercensus intervals, convert negative numbers to 0's (only interested in recruitment, not net change in saplings), calculate time adjusted recruitment (since inter-census intervals vary by 1-3 months)
census3<- census2 %>% 
  filter(year<2015) %>% 
  group_by(species, site, block, treatment, year, CensusDate) %>% 
  summarize(recruit = sum(combosap))

# calculate difference in number of saplings and date between years, filtering out rows for years 2009 and 2012 after calculation
census3<- census3 %>% 
  group_by(species, treatment, block) %>% 
 filter(year != 2011) 
census3 <- tidyr::pivot_wider(census3, id_cols= c(species, block, treatment), names_from= year, values_from= c(CensusDate, recruit))
df1 <- census3[, c("species", "block", "treatment", "CensusDate_2009", "CensusDate_2010", "recruit_2009", "recruit_2010")]
df2 <- census3[, c("species", "block", "treatment", "CensusDate_2012", "CensusDate_2013", "recruit_2012", "recruit_2013")]
df3 <- census3[, c("species", "block", "treatment", "CensusDate_2013", "CensusDate_2014", "recruit_2013", "recruit_2014")]
names(df1) <-names(df2) <- names(df3) <-c("species", "block","treatment", "initial_census_date", "next_census_date", "initial_recruit_no", "next_recruit_no")
census3 <- rbind(cbind(rep(2009, length.out= dim(census3)[1]), df1), 
                 cbind(rep(2012, length.out= dim(census3)[1]), df2), 
                 cbind(rep(2013, length.out= dim(census3)[1]),  df3))
names(census3) <- c("initial_year", "species", "block","treatment", "initial_census_date", "next_census_date", "initial_recruit_no", "next_recruit_no")
census3$diff.date <- as.numeric(census3$next_census_date - census3$initial_census_date)

# Fourth step - no data for BARO from North open, manually add rows for all 3 years w/ 0's for recruitment
  df1 <- as.data.frame(c(2009, 2009, 2009, 2012, 2012, 2012, 2013,2013,2013))
  names(df1) <- "initial_year"
  df1$species <- "BARO"
  df1$block <- c("N1", "N2", "N3", "S1", "S2", "S3", "C1", "C2","C3")
  df1$treatment <- "OPEN"
  df1$initial_census_date <-   df1$next_census_date <- NA
  df1$initial_recruit_no <- df1$next_recruit_no <- 0
  df1$diff.date <- 365
recruit_numbers <- rbind(census3, df1)

# write function to gather number of adults in preceding year (i.e., growth$year -1) by height class for each year-site-treatment-block combo
adulting <- function(species, binedges){
  no_adults <- array(NA, dim= c(length(vec_year_repro), length(vec_treatment), length(vec_block), (length(binedges)-1)))
  
  for(y in 1:length(vec_year_repro)){
    for(t in 1:length(vec_treatment)){
        for(b in 1:length(vec_block)){
        my_y <- vec_year_repro[y]
        my_t <- vec_treatment[t]
        my_b <- vec_block[b]
          binmids <- (binedges[1:(length(binedges)-1)] + binedges[2:(length(binedges))])/2# AML CHANGED
          temp_mat <- growth %>% filter(species == species & dead == 'N' & year == my_y-1 & treatment == my_t & block == my_b & ht_m >= 1 & ht_m <= max(binedges)) %>% select(ht_m) %>% tidyr::drop_na()  # there appear to be data transcription errors, so that sometimes the max height in "growth" is larger than the max height in the sACBR database... hence the max() argumnent
          adults<- hist(temp_mat$ht_m, # AML CHANGED
                        breaks= binedges, plot=FALSE)$counts# AML CHANGED
          no_adults[y, t, b,] <- adults
          dimnames(no_adults) <- list(vec_year_repro, vec_treatment, vec_block, binmids)
        }
    }
  }#end for loops
  return(no_adults)
}#end adulting function. 
# NB this function is unstable-- uses the dataframe "growth" WITHIN the function, but "growth" is not specified in the arguments. 
# I don't think the function should work b/c of this issue, but it does... 
#***#/// BRH  - growth is a global object that the function is pulling on. It's the initial data compilation at the beginning of the 1)vital_rate_regressions_v5 file. Since 'growth' itself isn't modified the function is stable.
n_adults_ACBR<- adulting(species = 'ACBR', binedges = all.binedges[['ACBR']])
n_adults_ACET<- adulting(species = 'ACET', binedges = all.binedges[['ACET']])
n_adults_ACME<- adulting(species = 'ACME', binedges = all.binedges[['ACME']])
n_adults_BARO<- adulting(species = 'BARO', binedges = all.binedges[['BARO']])
n_adults_CRDI<- adulting(species = 'CRDI', binedges = all.binedges[['CRDI']])
n_adults_all <- list(n_adults_ACBR, n_adults_ACET, n_adults_ACME, n_adults_BARO, n_adults_CRDI)

# predictor variables----- 
# Second calculate average rainfall and animal dung for each year-site and each year-site-treatment combo to differentiate different years in matrix construction. Remember that rainfall_py is scaled. Values are cumulative rainfall for 365 days Feb 1st in each year
rain_avg <- read.csv("xlong_rain_yearly_avg.csv", na.strings = 'NA') %>% 
  mutate(date = ymd(date)) %>% 
  mutate(year = year(date), month = month(date), day = day(date)) %>% 
  filter(month == 2 & day == 1) %>% 
  select(year, site, rainfall_py) %>% 
  filter(year<2015 & year>2009) %>% 
  mutate(rainfall_py = scale(rainfall_py, center = F, scale = T))

# calculate average dung in preceding year from Feb 1st. dung averaged across blocks in regression inputs, do so here too
dung_prep<-read.csv('dung/dung_surveys_2009-2020_clean_dec2020.csv') %>%
  select(survey:line,dikdik_old:elephant_new) %>% filter(treatment != 'OUT') %>% filter(year<2020)#NA's in 2020
dung_prep$treatment[dung_prep$treatment=='CTL']<-'OPEN'
dung_prep$treatment[dung_prep$treatment=='LMH']<-'TOTAL'
dung_prep_sum<- dung_prep %>% 
  mutate(year2 = ifelse(day>1 & month > 2, year + 1, year)) %>% 
  group_by(year2, site, block, treatment) %>%
  summarize(dikdik_dung_py = mean(dikdik_new), impala_dung_py = mean(impala_new), elephant_dung_py = mean(elephant_new)) %>%
  filter(year2 < 2015, block != 'CNA') %>% 
  rename(year = year2) %>% 
  mutate(year = as.character(year), block2 = paste0(site, block)) %>%
  ungroup() %>% 
  select(-c(block, site)) %>% 
  rename(block = block2)

all.random.s.coeff <- list( #note I have to use lists, bcuase the dimension of each of these is different 
 cbind( mvrnorm(n=1000, mu = fixef(ACBRs.best), Sigma= as.matrix(vcov(ACBRs.best))), as.data.frame(t(ranef(ACBRs.best)$block)), row.names = NULL), 
 cbind( mvrnorm(n=1000, mu = fixef(ACETs.best), Sigma= as.matrix(vcov(ACETs.best))), as.data.frame(t(ranef(ACETs.best)$block)), row.names = NULL), 
 cbind(  mvrnorm(n=1000, mu = fixef(ACMEs.best), Sigma= as.matrix(vcov(ACMEs.best))), as.data.frame(t(ranef(ACMEs.best)$block)), row.names = NULL), 
 NA, # for BARO, just used the mean 0.99.... there is no best fit model, either a glmer or glm)
 cbind(  mvrnorm(n=1000, mu = fixef(CRDIs.best), Sigma= as.matrix(vcov(CRDIs.best))), as.data.frame(t(ranef(CRDIs.best)$block)), row.names = NULL))

all.random.g.coeff <- list(
  cbind(mvrnorm(n=1000, mu = fixef(ACBRg.best), Sigma= as.matrix(vcov(ACBRg.best))), as.data.frame(t(ranef(ACBRg.best)$block)), row.names = NULL), 
  cbind(mvrnorm(n=1000, mu = fixef(ACETg.best), Sigma= as.matrix(vcov(ACETg.best))), as.data.frame(t(ranef(ACETg.best)$block)), row.names = NULL), 
  cbind(mvrnorm(n=1000, mu = fixef(ACMEg.best), Sigma= as.matrix(vcov(ACMEg.best))), as.data.frame(t(ranef(ACMEg.best)$block)), row.names = NULL), 
  mvrnorm(n=1000, mu = coefficients(BAROg.best), Sigma= as.matrix(vcov(BAROg.best))), # this is a lm
cbind(mvrnorm(n=1000, mu = fixef(CRDIg.best), Sigma= as.matrix(vcov(CRDIg.best))), as.data.frame(t(ranef(CRDIg.best)$block)), row.names = NULL))


all.random.res.coeff <- list(
cbind(mvrnorm(n=1000, mu = fixef(ACBRres.best), Sigma= as.matrix(vcov(ACBRres.best))), as.data.frame(t(ranef(ACBRres.best)$block)), row.names = NULL),  
mvrnorm(n=1000, mu = coefficients(ACETres.best), Sigma= as.matrix(vcov(ACETres.best))),# this is a lm
mvrnorm(n=1000, mu = coefficients(ACMEres.best), Sigma= as.matrix(vcov(ACMEres.best))),# this is a lm  
cbind(mvrnorm(n=1000, mu = fixef(BAROres.best), Sigma= as.matrix(vcov(BAROres.best))), as.data.frame(t(ranef(BAROres.best)$block)), row.names = NULL),
cbind(mvrnorm(n=1000, mu = fixef(CRDIres.best), Sigma= as.matrix(vcov(CRDIres.best))), as.data.frame(t(ranef(CRDIres.best)$block)), row.names = NULL))

all.random.pr.coeff <- list(
cbind(mvrnorm(n=1000, mu = fixef(ACBRpr.best), Sigma= as.matrix(vcov(ACBRpr.best))), as.data.frame(t(ranef(ACBRpr.best)$block)), row.names = NULL),   
cbind(mvrnorm(n=1000, mu = fixef(ACETpr.best), Sigma= as.matrix(vcov(ACETpr.best))), as.data.frame(t(ranef(ACETpr.best)$block)), row.names = NULL),  
cbind(mvrnorm(n=1000, mu = fixef(ACMEpr.best), Sigma= as.matrix(vcov(ACMEpr.best))), as.data.frame(t(ranef(ACMEpr.best)$block)), row.names = NULL),  
mvrnorm(n=1000, mu = coefficients(BAROpr.best), Sigma= as.matrix(vcov(BAROpr.best))),# this is a lm
cbind(mvrnorm(n=1000, mu = fixef(CRDIpr.best), Sigma= as.matrix(vcov(CRDIpr.best))), as.data.frame(t(ranef(CRDIpr.best)$block)), row.names = NULL))

all.random.sp.coeff <- list(
cbind(mvrnorm(n=1000, mu = fixef(ACBRsp.best), Sigma= as.matrix(vcov(ACBRsp.best))), as.data.frame(t(ranef(ACBRsp.best)$block)), row.names = NULL),
cbind(mvrnorm(n=1000, mu = fixef(ACETsp.best), Sigma= as.matrix(vcov(ACETsp.best))),as.data.frame(t(ranef(ACETsp.best)$block)), row.names = NULL),
mvrnorm(n=1000, mu = coefficients(ACMEsp.best), Sigma= as.matrix(vcov(ACMEsp.best))),# this is a glm
cbind(mvrnorm(n=1000, mu = fixef(BAROsp.best), Sigma= as.matrix(vcov(BAROsp.best))), as.data.frame(t(ranef(BAROsp.best)$block)), row.names = NULL),
cbind(mvrnorm(n=1000, mu = fixef(CRDIsp.best), Sigma= as.matrix(vcov(CRDIsp.best))),as.data.frame(t(ranef(CRDIsp.best)$block)), row.names = NULL))

min.seed.production <- c(min(ACBRsp.best@frame$seed), min(ACETsp.best@frame$seed), min(ACMEsp.best$data$seed), min(BAROsp.best@frame$seed), min(CRDIsp.best@frame$seed))
names(all.random.sp.coeff) <- names(all.random.pr.coeff) <-  names(all.random.res.coeff) <- names(all.random.g.coeff) <- names(all.random.s.coeff) <- vec_species


custom_predict_function <- function(forpredict, coeff, my_block){
  if (class(coeff) == "numeric") {coeff <- as.data.frame(t(coeff))}
  if("blockC1" %in% names(coeff)){ names(coeff)[which(names(coeff) == "blockC1")] <- "C1" }# if glm/ lm has block as an effect, the name of that coefficient is "blockC1" (and similar), so we need to replace it with C1 (or similar)
  if("blockC2" %in% names(coeff)){ names(coeff)[which(names(coeff) == "blockC2")] <- "C2" }
  if("blockC3" %in% names(coeff)){ names(coeff)[which(names(coeff) == "blockC3")] <- "C3" }
  if("blockN1" %in% names(coeff)){ names(coeff)[which(names(coeff) == "blockN1")] <- "N1" }
  if("blockN2" %in% names(coeff)){ names(coeff)[which(names(coeff) == "blockN2")] <- "N2" }
  if("blockN3" %in% names(coeff)){ names(coeff)[which(names(coeff) == "blockN3")] <- "N3" }
  if("blockS1" %in% names(coeff)){ names(coeff)[which(names(coeff) == "blockS1")] <- "S1" }
  if("blockS2" %in% names(coeff)){ names(coeff)[which(names(coeff) == "blockS2")] <- "S2" }
  if("blockS3" %in% names(coeff)){ names(coeff)[which(names(coeff) == "blockS3")] <- "S3" }
  
  possible_summation_terms <- cbind(forpredict$ht_t* ifelse("ht_t" %in% names(coeff), coeff$ht_t, 0), 
                                    ifelse("(Intercept)" %in% names(coeff), rep(coeff$"(Intercept)", length.out= length(forpredict$ht_t)), rep(0, length.out= length(forpredict$ht_t))), 
                                    ifelse("C1" %in% names(coeff) & my_block== "C1", rep(coeff$C1, length.out= length(forpredict$ht_t)), rep(0, length.out= length(forpredict$ht_t))), 
                                    ifelse("C2" %in% names(coeff) & my_block== "C2", rep(coeff$C2, length.out= length(forpredict$ht_t)), rep(0, length.out= length(forpredict$ht_t))), 
                                    ifelse("C3" %in% names(coeff) & my_block== "C3", rep(coeff$C3, length.out= length(forpredict$ht_t)), rep(0, length.out= length(forpredict$ht_t))), 
                                    ifelse("S1" %in% names(coeff) & my_block== "S1", rep(coeff$S1, length.out= length(forpredict$ht_t)), rep(0, length.out= length(forpredict$ht_t))), 
                                    ifelse("S2" %in% names(coeff) & my_block== "S2", rep(coeff$S2, length.out= length(forpredict$ht_t)), rep(0, length.out= length(forpredict$ht_t))), 
                                    ifelse("S3" %in% names(coeff) & my_block== "S3", rep(coeff$S3, length.out= length(forpredict$ht_t)), rep(0, length.out= length(forpredict$ht_t))), 
                                    ifelse("N1" %in% names(coeff) & my_block== "N1", rep(coeff$N1, length.out= length(forpredict$ht_t)), rep(0, length.out= length(forpredict$ht_t))), 
                                    ifelse("N2" %in% names(coeff) & my_block== "N2", rep(coeff$N2, length.out= length(forpredict$ht_t)), rep(0, length.out= length(forpredict$ht_t))), 
                                    ifelse("N3" %in% names(coeff) & my_block== "N3", rep(coeff$N3, length.out= length(forpredict$ht_t)), rep(0, length.out= length(forpredict$ht_t))), 
                                    forpredict$elephant_dung_py* ifelse("elephant_dung_py" %in% names(coeff), coeff$elephant_dung_py, 0), 
                                    forpredict$impala_dung_py* ifelse("impala_dung_py" %in% names(coeff), coeff$impala_dung_py, 0), 
                                    forpredict$dikdik_dung_py* ifelse("dikdik_dung_py" %in% names(coeff), coeff$dikdik_dung_py, 0), 
                                    forpredict$rainfall_py* ifelse("rainfall_py" %in% names(coeff), coeff$rainfall_py, 0))
  return(apply(possible_summation_terms, 1, sum))
}

 
#Now run the function

noreps <- 25

all_stlambdas <- array(NA, dim= c(length(vec_species), length(vec_treatment), length(vec_block), noreps))

tmax <- 10000

for (sp in 1:length(vec_species)){
  sp_name <- vec_species[sp]
  binedges_i <- all.binedges[[sp_name]]
  binmids_i <-  all.binmids[[sp_name]]
  all_mxes <- array(NA, dim= c(length(vec_treatment), length(vec_block), length(vec_year), noreps, (length(binmids_i)+1), (length(binmids_i)+1)))

for(t in 1:length(vec_treatment)){
  t_name <- vec_treatment[t]
  
for(b in 1:length(vec_block)){
  b_name <- vec_block[b]
  s <- substring(b_name, 1, 1)
  size_specific_fertility<- predicted_survival_all <- prob_surv_and_transition_to_other_classes <- array(NA, dim=c(length(vec_year), noreps, length(binmids_i)))
  prob_surv_and_stay_in_sap_class <- array(NA, dim=c(length(vec_year), noreps))
  survgmx <- array(NA, dim=c(length(vec_year), noreps, length(binmids_i), length(binmids_i)))
  
  for(y in 1:length(vec_year)){ # putting this year as the innermost loop because we'll have to average over some years later on in the loop
    y_name <- vec_year[y]
    
  for (i in 1: noreps){
    forpredict_i <- as.data.frame(c(binmids_sap, binmids_i))
    names(forpredict_i) <- "ht_t"
      
    forpredict_i$dikdik_dung_py = dung_prep_sum %>% filter(year == y_name & treatment == t_name & block == b_name ) %>% ungroup() %>% select(dikdik_dung_py) %>% nth(1)
    forpredict_i$impala_dung_py = dung_prep_sum %>% filter(year == y_name & treatment == t_name & block == b_name) %>% ungroup() %>% select(impala_dung_py) %>% nth(1) 
    forpredict_i$elephant_dung_py = dung_prep_sum %>% filter(year == y_name & treatment == t_name & block == b_name) %>% ungroup() %>% select(elephant_dung_py) %>% nth(1)
    forpredict_i$rainfall_py = as.numeric(rain_avg %>% filter(year == y_name & site == s) %>% ungroup() %>%  select(rainfall_py) %>% nth(1))
    
    if (sp_name== "BARO"){ predicted_survival <-  rep(0.9968652, length.out= length(forpredict_i$ht_t))} else {
      predicted_survival_eq <- custom_predict_function(forpredict = forpredict_i, coeff= all.random.s.coeff[[sp_name]][i,], my_block = b_name)
      predicted_survival <- 1/(1+exp(-predicted_survival_eq))}
  predicted_growth <- custom_predict_function(forpredict = forpredict_i, coeff= all.random.g.coeff[[sp_name]][i,], my_block = b_name)
  predicted_res <- custom_predict_function(forpredict = forpredict_i, coeff= all.random.res.coeff[[sp_name]][i,], my_block = b_name)
  predicted_psp_eq <- custom_predict_function(forpredict = forpredict_i, coeff= all.random.pr.coeff[[sp_name]][i,], my_block = b_name)
  predicted_psp <- 1/(1+exp(-predicted_psp_eq))
  predicted_sp_eq <- custom_predict_function(forpredict = forpredict_i, coeff= all.random.sp.coeff[[sp_name]][i,], my_block = b_name)
  predicted_sp <- (exp(predicted_sp_eq)-1)*forpredict_i$ht_t
  predicted_sp[which(predicted_sp <= 0)] <- ((exp(min.seed.production[sp])-1)*forpredict_i$ht_t)[which(predicted_sp <= 0)] 
  # methodological change: replacing any predicted seed numbers that are <= 0 with the minimum observed number for that species. this might not be required after hte reproductive vital rates are fit, but I'm leaving it in here for safety
  if (length( which(predicted_growth<0))>0) {stop("growth or res is negative")}
  predicted_res[which(predicted_res <0)] <- 0
   
  # making the survgrwoth matrix for binmids_sap
  predicted_survival_sap <- predicted_survival[1:length(binmids_sap)] 
  predicted_survival <- predicted_survival[(length(binmids_sap)+ 1 ): length(predicted_survival)] 
  predicted_growth_sap <- predicted_growth[1:length(binmids_sap)] 
  predicted_growth <- predicted_growth[(length(binmids_sap)+ 1 ): length(predicted_growth)] 
  predicted_res_sap <- predicted_res[1:length(binmids_sap)] 
  predicted_res <- predicted_res[(length(binmids_sap)+ 1 ): length(predicted_res)] 
  predicted_sp_sap <- predicted_sp[1:length(binmids_sap)] 
  predicted_sp <- predicted_sp[(length(binmids_sap)+ 1 ): length(predicted_sp)] 
  predicted_psp_sap <- predicted_psp[1:length(binmids_sap)] 
  predicted_psp <- predicted_psp[(length(binmids_sap)+ 1 ): length(predicted_psp)] 
  
  gmx <- matrix(data=NA, nrow=length(binmids_sap), ncol=length(binmids_sap))

  # get vital rates for saplings
  growcdf <- pnorm(c(binedges_sap, binedges_i[-1]), predicted_growth_sap, sqrt(predicted_res_sap)) 
  grows <- growcdf[2:length(growcdf)]-growcdf[1:(length(growcdf)-1)]
  grows <- grows/sum(grows)
  prob_surv_and_stay_in_sap_class[y, i] <- grows[1]*predicted_survival_sap
  prob_surv_and_transition_to_other_classes[y, i,] <- grows[-1]*predicted_survival_sap
 
  # making the matrix for the size-structed component of the matrix (i.e., heights >1)
  gmx <- matrix(data=NA, nrow=length(binmids_i), ncol=length(binmids_i))
  
  # get transition matrix using cdf fn
  for (ss in 1:length(binmids_i)) {
    growcdf <- pnorm(binedges_i,predicted_growth[ss],sqrt(predicted_res[ss]))
    grows <- growcdf[2:length(binedges_i)]-growcdf[1:(length(binedges_i)-1)]
    if(sum(grows)>0){grows <- grows/sum(grows)
    gmx[,ss] <- grows} else {gmx[,ss] <- NA} 
     # this if statement breaks the code (puts NA's into the matrix) if the sum of the PDF is zero (which happens if all the probability is outside of the size bounds)
  } # end ss loop
  
  # make the surv*growth mx 
  survgmx[y, i,,] <- gmx*t(matrix( rep(predicted_survival,length(binmids_i)),length(binmids_i))) # surv* growth matrix
  
  predicted_survival_all[y, i,] <- predicted_survival
if (y_name %in% c(2010, 2013,2014)){  
  initial_adults <- n_adults_all[[sp]][as.character(y_name), t_name, b_name,] # y_name is correct, because this n_adults array references the # of adults in the recruit_numbers$initial year
  sapling_row <- recruit_numbers[which(recruit_numbers$species== sp_name & (recruit_numbers$initial_year+1) == y_name & recruit_numbers$treatment == t_name & recruit_numbers$block == b_name),]
  
  no_new_saplings <- max((sapling_row$next_recruit_no - prob_surv_and_stay_in_sap_class[y, i]* sapling_row$initial_recruit_no) * #methodological change-- this change accounts for the fact that some individuals in the sapling class are OLD from last year; plus it calculates sapling number for each random draw of coeff estimates
              (365/sapling_row$diff.date),# corrects for diff date issues
              0) # max argument replaces any negative sapling # with a zero... 
  
  relative_fertility <- ifelse(sum(predicted_survival*predicted_psp *predicted_sp) ==0, # if no adult size classes are reproducing
                               rep(1/length(predicted_survival),length.out= length(predicted_survival) ), # make all relative fertilities equal
                               (predicted_survival*predicted_psp *predicted_sp)/ sum(predicted_survival*predicted_psp *predicted_sp) )# otherwise, calc relative fertility in the normal way (which will give you an NaN if no adult size classes are reproducting)
  # methodolgoical change: i have included survival of adults n the relative fertility estimates
  mean_fert <- no_new_saplings/ sum(relative_fertility*initial_adults)
  size_specific_fertility[y, i,] <- relative_fertility*mean_fert } # ends year if statement
    } # end i loop
  } # end  year loop
  
  size_specific_fertility[which(vec_year== 2011),  ,] <-  size_specific_fertility[which(vec_year== 2012),  ,] <- 
    (size_specific_fertility[which(vec_year== 2010),  ,] + size_specific_fertility[which(vec_year== 2013),  ,] + size_specific_fertility[which(vec_year== 2014),  ,] )/3 
  # replacing years 2011, 2012 fertiilty (2010-2011 interval and 2011-2012 interval) with the average over the other years 
  
  for(y in 1:length(vec_year)){ 
    for (i in 1: noreps){
      y_name <- vec_year[y]
  reprow <-  c(prob_surv_and_stay_in_sap_class[y,i], size_specific_fertility[y, i,]*predicted_survival_all[y,i, ])
  mx <- rbind(reprow, cbind(prob_surv_and_transition_to_other_classes[y, i,], survgmx[y,i,,]))
  all_mxes[t, b, y, i,,] <- mx
  
}} # ends y and i loops
  # now we have all our matrices stored so we can now project forward
 
  for (i in 1:noreps){
    
    n <- popbio::stable.stage(apply(all_mxes[t, b, , i,,] , MARGIN= c(2,3), FUN=mean)) # gets SSD of across-year mean matrix
    r <- rep(NA, length= tmax)
    for (t_run in 1:tmax){
      n <- all_mxes[t, b, sample(1:length(vec_year), 1), i,,]%*%n
      N <- sum(n)
      r[t_run] <-log(N)
      n <-n/N
    }
    st_lambda <- mean(r[2000:tmax])
  all_stlambdas[sp, t, b, i] <- st_lambda
    } # closes i loop
  print(paste(sp_name, t_name, b_name, "all reps done", sep= " "))
}}
}


# Sensitivities Code ####
# NB this code uses the same mvrnorm selected values, which is v important


all_stlambdas.sens <- array(NA, dim = c(length(vec_species), length(vec_treatment), length(vec_block), noreps, 13,2))
all_change_in_rates <- array(NA, dim= c(13, length(vec_species), length(vec_treatment), length(vec_block), noreps, length(vec_year),2))
all_rates <- array(NA, dim= c(13, length(vec_species), length(vec_treatment), length(vec_block), noreps, length(vec_year)))
                           
for (v in 1:13){ # order here is: 
      # v= 1 rainfall
      # v= 2 dik dik dung
      # v= 3 impala dung
      # v= 4 elephant dung
      # v=5 survival of adults
      # v= 6 growth of adults
      # v= 7 res of adults
      # v = 8 psp of adults
      # v= 9 sp of adults
      # v= 10 no.new.saplings
      # v=11 survival of saps
      # v= 12 growth of saps
      # v= 13 res of saps 
    for (p in 1:2){ # first is 0.05 up, second is 0.05 down
      perturb_fraction= c(1.05, 0.95)[p]
      for (sp in 1:length(vec_species)){
        sp_name <- vec_species[sp]
        binedges_i <- all.binedges[[sp_name]]
        binmids_i <-  all.binmids[[sp_name]]
        all_mxes <- array(NA, dim= c(length(vec_treatment), length(vec_block), length(vec_year), noreps, (length(binmids_i)+1), (length(binmids_i)+1)))
        
        for(t in 1:length(vec_treatment)){
          t_name <- vec_treatment[t]
          
          for(b in 1:length(vec_block)){
            b_name <- vec_block[b]
            s <- substring(b_name, 1, 1)
            size_specific_fertility<- predicted_survival_all <- prob_surv_and_transition_to_other_classes <- array(NA, dim=c(length(vec_year), noreps, length(binmids_i)))
            prob_surv_and_stay_in_sap_class <- array(NA, dim=c(length(vec_year), noreps))
            survgmx <- array(NA, dim=c(length(vec_year), noreps, length(binmids_i), length(binmids_i)))
            
            number_new_saplings_out<- array(NA, dim = c(length(vec_year), noreps))
            
            for(y in 1:length(vec_year)){ # putting this year as the innermost loop because we'll have to average over some years later on in the loop
              y_name <- vec_year[y]
              
              for (i in 1: (noreps)){
                forpredict_i <- as.data.frame(c(binmids_sap, binmids_i))
                names(forpredict_i) <- "ht_t"
                
                forpredict_i$dikdik_dung_py = dung_prep_sum %>% filter(year == y_name & treatment == t_name & block == b_name ) %>% ungroup() %>% select(dikdik_dung_py) %>% nth(1)
                forpredict_i$impala_dung_py = dung_prep_sum %>% filter(year == y_name & treatment == t_name & block == b_name) %>% ungroup() %>% select(impala_dung_py) %>% nth(1) 
                forpredict_i$elephant_dung_py = dung_prep_sum %>% filter(year == y_name & treatment == t_name & block == b_name) %>% ungroup() %>% select(elephant_dung_py) %>% nth(1)
                forpredict_i$rainfall_py = as.numeric(rain_avg %>% filter(year == y_name & site == s) %>% ungroup() %>%  select(rainfall_py) %>% nth(1))
                # elasticity/ sens code (environmental covariates)
                { 
                  if (v== 1){
                    all_rates[v,sp, t, b,i,y] <- unique(forpredict_i$rainfall_py)
                    forpredict_i$rainfall_py <- forpredict_i$rainfall_py*perturb_fraction
                    all_change_in_rates[v, sp, t, b,i, y, p] <- unique(forpredict_i$rainfall_py)-all_rates[v,sp, t, b,i,y]
                    }
                  if (v== 2){all_rates[v,sp, t, b,i,y] <- unique(forpredict_i$dikdik_dung_py)
                    forpredict_i$dikdik_dung_py <- forpredict_i$dikdik_dung_py*perturb_fraction
                    all_change_in_rates[v, sp, t, b,i, y,p] <- unique(forpredict_i$dikdik_dung_py)-all_rates[v,sp, t, b,i,y] 
                      if(all_change_in_rates[v, sp, t, b,i, y, p] == 0){
                      all_change_in_rates[v, sp, t, b,i, y, p]<-
                        dung_prep_sum %>% filter(treatment == 'OPEN') %>% 
                        summarize(mean = mean(dikdik_dung_py)*0.05) %>% 
                        pull()
                        if(p==1){
                          forpredict_i$dikdik_dung_py<- 
                            all_change_in_rates[v, sp, t, b,i, y, p]
                        }}#make small artificial perturbation for cases where there's no change in rates in data. Only perform perturb up, perturb down remains 0
                    }
                  if (v== 3){all_rates[v,sp, t, b,i,y] <- unique(forpredict_i$impala_dung_py)
                    forpredict_i$impala_dung_py <- forpredict_i$impala_dung_py*perturb_fraction
                    all_change_in_rates[v, sp, t, b,i, y,p] <- unique(forpredict_i$impala_dung_py)-all_rates[v,sp, t, b,i,y] 
                    if(all_change_in_rates[v, sp, t, b,i, y, p] == 0){
                      all_change_in_rates[v, sp, t, b,i, y, p]<-
                        dung_prep_sum %>% filter(treatment == 'OPEN') %>% 
                        summarize(mean = mean(impala_dung_py)*0.05) %>% 
                        pull()
                      if(p==1){
                        forpredict_i$impala_dung_py<- 
                          all_change_in_rates[v, sp, t, b,i, y, p]
                      }}
                    }
                  if (v== 4){all_rates[v,sp, t, b,i,y] <- unique(forpredict_i$elephant_dung_py)
                    forpredict_i$elephant_dung_py <- forpredict_i$elephant_dung_py*perturb_fraction
                    all_change_in_rates[v, sp, t, b,i, y,p] <- unique(forpredict_i$elephant_dung_py)-all_rates[v,sp, t, b,i,y] 
                    if(all_change_in_rates[v, sp, t, b,i, y, p] == 0){
                      all_change_in_rates[v, sp, t, b,i, y, p]<-
                        dung_prep_sum %>% filter(treatment == 'OPEN') %>% 
                        summarize(mean = mean(elephant_dung_py)*0.05) %>% 
                        pull()
                      if(p==1){
                        forpredict_i$elephant_dung_py<- 
                          all_change_in_rates[v, sp, t, b,i, y, p]
                      }}
                  }
                }#v1:V4 grouping 
                if (sp_name== "BARO"){ predicted_survival <-  rep(0.9968652, length.out= length(forpredict_i$ht_t))} else {
                  predicted_survival_eq <- custom_predict_function(forpredict = forpredict_i, coeff= all.random.s.coeff[[sp_name]][i,], my_block = b_name)
                  predicted_survival <- 1/(1+exp(-predicted_survival_eq))}
                predicted_growth <- custom_predict_function(forpredict = forpredict_i, coeff= all.random.g.coeff[[sp_name]][i,], my_block = b_name)
                predicted_res <- custom_predict_function(forpredict = forpredict_i, coeff= all.random.res.coeff[[sp_name]][i,], my_block = b_name)
                predicted_psp_eq <- custom_predict_function(forpredict = forpredict_i, coeff= all.random.pr.coeff[[sp_name]][i,], my_block = b_name)
                predicted_psp <- 1/(1+exp(-predicted_psp_eq))
                predicted_sp_eq <- custom_predict_function(forpredict = forpredict_i, coeff= all.random.sp.coeff[[sp_name]][i,], my_block = b_name)
                predicted_sp <- (exp(predicted_sp_eq)-1)*forpredict_i$ht_t
                predicted_sp[which(predicted_sp <= 0)] <- ((exp(min.seed.production[sp])-1)*forpredict_i$ht_t)[which(predicted_sp <= 0)] 
                # methodological change: replacing any predicted seed numbers that are <= 0 with the minimum observed number for that species. this might not be required after the reproductive vital rates are fit, but I'm leaving it in here for safety
                if (length( which(predicted_growth<0))>0) {stop("growth or res is negative")}
                predicted_res[which(predicted_res <0)] <- 0
                
                # making the survgrwoth matrix for binmids_sap
                predicted_survival_sap <- predicted_survival[1:length(binmids_sap)] 
                predicted_survival <- predicted_survival[(length(binmids_sap)+ 1 ): length(predicted_survival)] 
                predicted_growth_sap <- predicted_growth[1:length(binmids_sap)] 
                predicted_growth <- predicted_growth[(length(binmids_sap)+ 1 ): length(predicted_growth)] 
                predicted_res_sap <- predicted_res[1:length(binmids_sap)] 
                predicted_res <- predicted_res[(length(binmids_sap)+ 1 ): length(predicted_res)] 
                predicted_sp <- predicted_sp[(length(binmids_sap)+ 1 ): length(predicted_sp)] 
                predicted_psp <- predicted_psp[(length(binmids_sap)+ 1 ): length(predicted_psp)] 
                # implicit in these matrices is that individuasl <1 m don't reproduce
                
# elasticities/ sensitivites perturbation for vital rates
                {             
                if (v== 5){ 
                  store_predicted_survival <- predicted_survival 
                  all_rates[v,sp, t, b,i,y]  <- mean(predicted_survival) # average survival across all adult size classes 
                  predicted_survival <- predicted_survival*perturb_fraction
                  predicted_survival[ which(predicted_survival > 1)] <- 1 # just in case the perturbation pushes survival to >1
                  all_change_in_rates[v, sp, t, b,i, y,p] <- mean(predicted_survival-store_predicted_survival) # average perturbation across all adult size classes 
                }
                if (v== 6){
                  store_predicted_growth <- predicted_growth
                  all_rates[v,sp, t, b,i,y]  <- mean(predicted_growth)# average growth across all adult size classes 
                  predicted_growth <- predicted_growth*perturb_fraction
                  all_change_in_rates[v, sp, t, b,i, y,p] <- mean(predicted_growth- store_predicted_growth) # average perturbation across all adult size classes 
                }
                if (v== 7){
                  store_predicted_res <- predicted_res
                  all_rates[v,sp, t, b,i,y]  <- mean(predicted_res)# average res across all adult size classes 
                  predicted_res <- predicted_res*perturb_fraction
                  all_change_in_rates[v, sp, t, b,i, y,p] <- mean(predicted_res-store_predicted_res) # average perturbation across all adult size classes 
                }
                if (v== 8){
                  store_predicted_psp <- predicted_psp
                  all_rates[v,sp, t, b,i,y]  <- mean(predicted_psp)# average psp across all adult size classes 
                  predicted_psp <- predicted_psp*perturb_fraction
                  predicted_psp[ which(predicted_psp > 1)] <- 1 # just in case the perturbation pushes survival to >1
                  all_change_in_rates[v, sp, t, b, i, y, p] <- mean(predicted_psp-store_predicted_psp) # average perturbation across all adult size classes 
                }
                if (v== 9){
                  store_predicted_sp <- predicted_sp
                  all_rates[v, sp, t, b, i, y]  <- mean(predicted_sp)# average sp across all adult size classes 
                  predicted_sp <- predicted_sp*perturb_fraction
                  all_change_in_rates[v, sp, t, b, i, y, p] <- mean(predicted_sp-store_predicted_sp) # average perturbation across all adult size classes 
                }
                if (v== 11){
                  store_predicted_survival_sap <- predicted_survival_sap
                  all_rates[v,sp, t, b,i,y]  <- predicted_survival_sap
                  predicted_survival_sap <- max(1, predicted_survival_sap*perturb_fraction) # limits survival of saplings to one
                  all_change_in_rates[v, sp, t, b,i, y,p] <- predicted_survival_sap- store_predicted_survival_sap
                }
                if (v== 12){
                  store_predicted_growth_sap <- predicted_growth_sap
                  all_rates[v,sp, t, b,i,y]  <- predicted_growth_sap
                  predicted_growth_sap <- predicted_growth_sap*perturb_fraction
                  all_change_in_rates[v, sp, t, b,i, y,p] <- predicted_growth_sap-store_predicted_growth_sap
                }
                if (v== 13){
                  store_predicted_res_sap <- predicted_res_sap
                  all_rates[v,sp, t, b,i,y]  <- predicted_res_sap
                  predicted_res_sap <- predicted_res_sap*perturb_fraction
                  all_change_in_rates[v, sp, t, b,i, y,p] <- predicted_res_sap-store_predicted_res_sap
                  if(all_change_in_rates[v, sp, t, b,i, y, p] == 0){
                    all_change_in_rates[v, sp, t, b,i, y, p]<-
                      min(predicted_res[which(predicted_res>0)])*0.05
                    if(p==1){
                      predicted_res_sap<-
                        all_change_in_rates[v, sp, t, b,i, y, p]
                    }}#provide small artificial perturbation in cases where there's no change in data. Only perform perturb for up, perturb down remains 0
                }
                
              }
                gmx <- matrix(data=NA, nrow=length(binmids_sap), ncol=length(binmids_sap))
                
                # get vital rates for saplings: methodological change
                growcdf <- pnorm(c(binedges_sap, binedges_i[-1]),predicted_growth_sap,sqrt(predicted_res_sap)) # check that res^2 is right here
                grows <- growcdf[2:length(growcdf)]-growcdf[1:(length(growcdf)-1)]
                grows <- grows/sum(grows)
                prob_surv_and_stay_in_sap_class[y, i] <- grows[1]*predicted_survival_sap
                prob_surv_and_transition_to_other_classes[y, i,] <- grows[-1]*predicted_survival_sap
                
                # making the matrix for the size-structed component of the matrix (i.e., heights >1)
                gmx <- matrix(data=NA, nrow=length(binmids_i), ncol=length(binmids_i))
                
                # get transition matrix using cdf fn
                for (ss in 1:length(binmids_i)) {
                  growcdf <- pnorm(binedges_i,predicted_growth[ss],sqrt(predicted_res[ss]))# check that res^2 is right here
                  grows <- growcdf[2:length(binedges_i)]-growcdf[1:(length(binedges_i)-1)]
                  if(sum(grows)>0){grows <- grows/sum(grows)
                  gmx[,ss] <- grows} else {gmx[,ss] <- NA} 
                  # this if statement breaks the code (puts NA's into the matrix) if the sum of the PDF is zero (which happens if all the probability is outside of the size bounds)
                } # end ss loop
                
                # make the surv*growth mx 
                survgmx[y, i,,] <- gmx*t(matrix( rep(predicted_survival,length(binmids_i)),length(binmids_i))) # surv* growth matrix
                
                predicted_survival_all[y, i,] <- predicted_survival
                
                if (y_name %in% c(2010, 2013,2014)){ 
                  initial_adults <- n_adults_all[[sp]][as.character(y_name), t_name, b_name,] # y_name is correct, because this n_adults array references the # of adults in the recruit_numbers$initial year
                  sapling_row <- recruit_numbers[which(recruit_numbers$species== sp_name & (recruit_numbers$initial_year+1) == y_name & recruit_numbers$treatment == t_name & recruit_numbers$block == b_name ),]
                  
                  no_new_saplings <- max((sapling_row$next_recruit_no - prob_surv_and_stay_in_sap_class[y, i]* sapling_row$initial_recruit_no) * #methodological change-- this change accounts for the fact that some individuals in the sapling class are OLD from last year; plus it calculates sapling number for each random draw of coeff estimates
                                           (365/sapling_row$diff.date),# corrects for diff date issues
                                         0) # max argument replaces any negative sapling # with a zero... 
                  
                  number_new_saplings_out[y,i]<-no_new_saplings #collect # new saplings so years w/ census data can be averaged to substitute for non-census years outside of loop
                if (v== 10){
                  store_no_new_saplings <- no_new_saplings
                    all_rates[v,sp, t, b,i,y] <- no_new_saplings
                    no_new_saplings <- no_new_saplings*perturb_fraction
                    all_change_in_rates[v, sp, t, b,i, y,p] <- no_new_saplings - store_no_new_saplings
                    if(all_change_in_rates[v, sp, t, b,i, y,p] == 0){
                      sapling_row0 <- recruit_numbers[which(recruit_numbers$species== sp_name & (recruit_numbers$initial_year+1) == y_name & recruit_numbers$treatment == t_name),] #if no recruits to perturb in a given species-treatment-year-block, take average of all blocks in the year after removing negative differences (tried taking average of all years in block but still had results w/ no recruitment)
                      no_new_saplings0 <- sapling_row0$next_recruit_no - prob_surv_and_stay_in_sap_class[y, i]* sapling_row0$initial_recruit_no *
                        (365/sapling_row0$diff.date)
                      no_new_saplings0 <- mean(no_new_saplings0[which(no_new_saplings0>0)]) * 0.05                     
                      if(is.nan(no_new_saplings0)){ #for cases where there's still no recruitment after averaging across blocks, average across years too
                        sapling_row0 <- recruit_numbers[which(recruit_numbers$species== sp_name & recruit_numbers$treatment == t_name),] 
                        no_new_saplings0 <- sapling_row0$next_recruit_no - prob_surv_and_stay_in_sap_class[y, i]* sapling_row0$initial_recruit_no *
                          (365/sapling_row0$diff.date)
                        no_new_saplings0 <- mean(no_new_saplings0[which(no_new_saplings0>0)]) * 0.05
                      }
                      all_change_in_rates[v, sp, t, b,i, y,p] <- 
                          no_new_saplings0
                      if(p==1){
                        no_new_saplings<- no_new_saplings0
                        number_new_saplings_out[y,i]<-no_new_saplings
                    }}
                  }
                  
                  relative_fertility <- ifelse(sum(predicted_survival*predicted_psp *predicted_sp) == 0, # if no adult size classes are reproducing
                                               rep(1/length(predicted_survival), length.out = length(predicted_survival) ), # make all relative fertilities equal
                                               (predicted_survival*predicted_psp *predicted_sp)/ sum(predicted_survival*predicted_psp *predicted_sp) )# otherwise, calc relative fertility in the normal way (which will give you an NaN if no adult size classes are reproducing)
                  # methodological change: i have included survival of adults n the relative fertility estimates
                  mean_fert <- no_new_saplings/ sum(relative_fertility*initial_adults)
                  size_specific_fertility[y, i,] <- relative_fertility*mean_fert } # ends year if statement
                
              } # end first i loop
            } #end first y loop

            size_specific_fertility[which(vec_year== 2011),  ,] <-  size_specific_fertility[which(vec_year== 2012),  ,] <- 
              (size_specific_fertility[which(vec_year== 2010),  ,] + size_specific_fertility[which(vec_year== 2013),  ,] + size_specific_fertility[which(vec_year== 2014),  ,] )/3 
            # replacing years 2011, 2012 fertiilty (2010-2011 interval and 2011-2012 interval) with the average over the other others 
            
            number_new_saplings_out[which(vec_year== 2011),] <-    
              number_new_saplings_out[which(vec_year== 2012),] <- 
              (number_new_saplings_out[which(vec_year== 2010),] + number_new_saplings_out[which(vec_year== 2013),] + number_new_saplings_out[which(vec_year== 2014),])/3  # do the same for # new saplings
            
            if (v== 10){ #now replace empty values in vital rates matrix from non-census year w/ averaged values of census years just created
              all_rates[v,sp, t, b, , which(vec_year == 2011)] <-
              all_rates[v,sp, t, b, , which(vec_year == 2012)] <- 
                number_new_saplings_out[which(vec_year == 2011),]
              #non-census year change in rate should be average of census-year change in rates to accomodate earlier alterations made to combos w/ 0 recruitment, rather than calculating perturbance here 
          all_change_in_rates[v, sp, t, b, , which(vec_year == 2011), p] <-
          all_change_in_rates[v, sp, t, b, , which(vec_year == 2012), p] <- 
           (all_change_in_rates[v, sp, t, b, , which(vec_year == 2010), p] +
            all_change_in_rates[v, sp, t, b, , which(vec_year == 2013), p] +
            all_change_in_rates[v, sp, t, b, , which(vec_year == 2014), p])/3
            }
           
            for(y in 1:length(vec_year)){ 
              for (i in 1: noreps){
                y_name <- vec_year[y]
                # finally-- do sensitivities, keeping careful track of the issue with parameter estimates
                reprow <-  c(prob_surv_and_stay_in_sap_class[y,i], size_specific_fertility[y, i,]*predicted_survival_all[y,i, ])
                mx <- rbind(reprow, cbind(prob_surv_and_transition_to_other_classes[y, i,], survgmx[y,i,,]))
                all_mxes[t, b, y, i,,] <- mx
                
              }} # ends y and i loops
            # now we have all our matrices stored so we can now project forward
              for (i in 1:noreps){
              
              n <- popbio::stable.stage(apply(all_mxes[t, b, , i,,] , MARGIN= c(2,3), FUN=mean)) # gets SSD of across-year mean matrix
              r <- rep(NA, length= tmax)
              for (t_run in 1:tmax){
                n <-all_mxes[t, b, sample(1:length(vec_year), 1), i,,]%*%n
                N <- sum(n)
                r[t_run] <-log(N)
                n <-n/N
              }
              st_lambda <- mean(r[2000:tmax])
              all_stlambdas.sens[sp, t, b, i, v, p] <- st_lambda
            } # closes final i loop
            print(paste(v, sp_name, t_name, #b_name, 
                        "all reps done", sep= " "))
          } #end b loop
        } #end t loop
      } #end sp loop
    } #end p loop
  } # end v loop
 

# Results ####
# to make figures, here is how to get mean & CI on stochastic lambdas for a given treatment x species
 
all_stlambdas_avrg_across_blocks <- apply(all_stlambdas, c(1,2,4), mean) # I recommend averaging stochastic lambda across blocks first
all_stlambdas_medianCI <- apply(all_stlambdas_avrg_across_blocks, c(1,2), quantile, c(0.05/2, 0.5, 1-0.05/2)) # then calculating the mean and CI across the parameter estimates
# the dimensions of this should correspond to species  and treatment

stlambda_results <- reshape2::melt(all_stlambdas_medianCI)
str(stlambda_results)

names(stlambda_results) <- c("parameter", "species", "treatment", "value")
stlambda_results$parameter <- dplyr::case_when(stlambda_results$parameter == "2.5%" ~ "LCL",
                                           stlambda_results$parameter == "50%" ~ "EST",
                                           stlambda_results$parameter == "97.5%" ~ "UCL")
stlambda_results$species <- dplyr::case_when(stlambda_results$species == 1 ~ "ACBR",
                                         stlambda_results$species == 2 ~ "ACET",
                                         stlambda_results$species == 3 ~ "ACME",
                                         stlambda_results$species == 4 ~ "BARO",
                                         stlambda_results$species == 5 ~ "CRDI")
stlambda_results$treatment <- dplyr::case_when(stlambda_results$treatment == 1 ~ "TOTAL",
                                           stlambda_results$treatment == 2 ~ "MEGA",
                                           stlambda_results$treatment == 3 ~ "MESO",
                                           stlambda_results$treatment == 4 ~ "OPEN")
 
# here is how to get mean & CI on sensitivities
all_rates2 <-  apply(all_rates, c(1, 2, 3, 4, 5), mean) # average across years. 
all_rates2 <- aperm(all_rates2, c(2, 3, 4, 5, 1)) # permutes array to be in same dimensions as the lambda vector
#repeat for changes in rates
all_change_in_rates2<- apply(all_change_in_rates, c(1, 2, 3, 4, 5, 7), mean)
all_change_in_rates2<- aperm(all_change_in_rates2, c(2, 3, 4, 5, 1, 6))
 
diff_up <- (sweep(x= all_stlambdas.sens[,,,,,1], # dim is species, treatment, block, rep, vital rate - all within perturbation up
       MARGIN= c(1,2,3,4) ,
       STATS= all_stlambdas, # dim is species, treatment, block, rep
       FUN= "-"))/all_change_in_rates2[,,,,,1] # the diff_up array has dimensions species, treatment, block, reps, vital rates
# diff_up[which(is.nan(diff_up))] <- 0 
# instances of no change in rates result in 0/0 = NaN, replace with 0's
# diff_up[which(is.infinite(diff_up))] <- 0 
# have instances where there is a difference in lambda values but no difference in rates, resulting in infinite values in diff_up/diff_down. Chalking this up to rounding error in lambda calulations, since 
# sum(length(diff_down[which(is.nan(diff_down))]) + length(diff_down[which(is.infinite(diff_down))])) == length(all_change_in_rates2[which(all_change_in_rates2[,,,,,2] == 0)] )

 
diff_down <- (sweep(x= all_stlambdas.sens[,,,,,2], # repeat for perturbation down
                   MARGIN= c(1,2,3,4) ,
                   STATS= all_stlambdas, 
                   FUN= "-"))/all_change_in_rates2[,,,,,2]  
# diff_down[which(is.nan(diff_down))] <- 0
# diff_down[which(is.infinite(diff_down))] <- 0 

sens_avrg_across_blocks <- apply( (diff_up + diff_down)/2, c(1,2,4,5), mean) # getting across-block average sensitivity; dims are species, treatment, reps, vital rates
sens_meanCI<- apply(sens_avrg_across_blocks, c(1,2,4), quantile, c(0.05/2, 0.5, 1-0.05/2))# then calculating the median and CI on sensitivities across the parameter estimates
# the dimensions of this should correspond to species, treatment, then vital rates 1-13

sens_results <- reshape2::melt(sens_meanCI)
str(sens_results)

names(sens_results) <- c("parameter", "species", "treatment", "vital_rate", "value")
sens_results$parameter <- dplyr::case_when(sens_results$parameter == "2.5%" ~ "LCL",
                                           sens_results$parameter == "50%" ~ "EST",
                                           sens_results$parameter == "97.5%" ~ "UCL")
sens_results$species <- dplyr::case_when(sens_results$species == 1 ~ "ACBR",
                                         sens_results$species == 2 ~ "ACET",
                                         sens_results$species == 3 ~ "ACME",
                                         sens_results$species == 4 ~ "BARO",
                                         sens_results$species == 5 ~ "CRDI")
sens_results$treatment <- dplyr::case_when(sens_results$treatment == 1 ~ "TOTAL",
                                         sens_results$treatment == 2 ~ "MEGA",
                                         sens_results$treatment == 3 ~ "MESO",
                                         sens_results$treatment == 4 ~ "OPEN")
sens_results$vital_rate <- dplyr::case_when(sens_results$vital_rate == 1 ~ "rainfall",
                                         sens_results$vital_rate == 2 ~ "dik-dik dung",
                                         sens_results$vital_rate == 3 ~ "impala dung",
                                         sens_results$vital_rate == 4 ~ "elephant dung",
                                         sens_results$vital_rate == 5 ~ "adult survival",
                                         sens_results$vital_rate == 6 ~ "adult growth",
                                         sens_results$vital_rate == 7 ~ "adult res",
                                         sens_results$vital_rate == 8 ~ "adult psp",
                                         sens_results$vital_rate == 9 ~ "adult sp",
                                         sens_results$vital_rate == 10 ~ "no of new saplings",
                                         sens_results$vital_rate == 11 ~ "sapling survival",
                                         sens_results$vital_rate == 12 ~ "sapling growth",
                                         sens_results$vital_rate == 13 ~ "sapling res")
 
# here is how to get mean & CI on elasticities
diff_up_elas <- diff_up *sweep(all_rates2, c(1,2,3,4), all_stlambdas, "/")
diff_down_elas <- diff_down *sweep(all_rates2, c(1,2,3,4), all_stlambdas, "/")

elas_avrg_across_blocks <- apply( (diff_up_elas + diff_down_elas)/2, c(1,2,4,5), mean) # getting across-block average elas; dims are species, treatment, reps, vital rates
elas_meanCI<- apply(elas_avrg_across_blocks, c(1,2,4), quantile, c(0.05/2, 0.5, 1-0.05/2))# then calculating the median and CI on elas across the parameter estimates
# the dimensions of this should correspond to species, treatment, then vital rates 1-13

elas_results <- reshape2::melt(elas_meanCI)
str(elas_results)

names(elas_results) <- c("parameter", "species", "treatment", "vital_rate", "value")
elas_results$parameter <- dplyr::case_when(elas_results$parameter == "2.5%" ~ "LCL",
                                           elas_results$parameter == "50%" ~ "EST",
                                           elas_results$parameter == "97.5%" ~ "UCL")
elas_results$species <- dplyr::case_when(elas_results$species == 1 ~ "ACBR",
                                         elas_results$species == 2 ~ "ACET",
                                         elas_results$species == 3 ~ "ACME",
                                         elas_results$species == 4 ~ "BARO",
                                         elas_results$species == 5 ~ "CRDI")
elas_results$treatment <- dplyr::case_when(elas_results$treatment == 1 ~ "TOTAL",
                                           elas_results$treatment == 2 ~ "MEGA",
                                           elas_results$treatment == 3 ~ "MESO",
                                           elas_results$treatment == 4 ~ "OPEN")
elas_results$vital_rate <- dplyr::case_when(elas_results$vital_rate == 1 ~ "rainfall",
                                            elas_results$vital_rate == 2 ~ "dik-dik dung",
                                            elas_results$vital_rate == 3 ~ "impala dung",
                                            elas_results$vital_rate == 4 ~ "elephant dung",
                                            elas_results$vital_rate == 5 ~ "adult survival",
                                            elas_results$vital_rate == 6 ~ "adult growth",
                                            elas_results$vital_rate == 7 ~ "adult res",
                                            elas_results$vital_rate == 8 ~ "adult psp",
                                            elas_results$vital_rate == 9 ~ "adult sp",
                                            elas_results$vital_rate == 10 ~ "no of new saplings",
                                            elas_results$vital_rate == 11 ~ "sapling survival",
                                            elas_results$vital_rate == 12 ~ "sapling growth",
                                            elas_results$vital_rate == 13 ~ "sapling res")

write.table(stlambda_results, 'results/stochastic_lambdas.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',')

write.table(sens_results, 'results/sensitivities.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',')

write.table(elas_results, 'results/elasticities.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',')
