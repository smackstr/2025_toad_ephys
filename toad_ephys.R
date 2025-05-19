# Script To-Do List --------------------------------------------------------------
#1. Decide on analysis approach, specifically model-selection approach
#2. Create exportable table for Supp Material -- use R Markdown for this? stargazer seems to be a thing
# see https://debyeeneuro.com/wp-content/uploads/2016/03/table_workshop.pdf
#stargazer(lmm.full.1, lmm.1intx.1, lmm.1intx.2, lmm.1intx.3, type = "html", out="test.doc", intercept.bottom = F, intercept.top = T, digits = 2)

# Load required packages and set working directory --------------------------------------------------------------

rm(list=ls())
required_pckg = c("dplyr", "readr", "NatParksPalettes", "ggpubr", "RColorBrewer",
                  "lme4", "MuMIn", "DHARMa", "stats", "latticeExtra", "cAIC4")
lapply(required_pckg, library, character.only=TRUE)

wd = "~/Desktop/R Working Directory/Databases"
setwd(wd)

# Quick fix for stargazer <= 5.2.3 is.na() issue with long model names in R >= 4.2 from https://gist.github.com/alexeyknorre/b0780836f4cec04d41a863a683f91b53 -------------
# Unload stargazer if loaded
detach("package:stargazer",unload=T)
# Delete it
remove.packages("stargazer")
# Download the source
download.file("https://cran.r-project.org/src/contrib/stargazer_5.2.3.tar.gz", destfile = "stargazer_5.2.3.tar.gz")
# Unpack
untar("stargazer_5.2.3.tar.gz")
# Read the sourcefile with .inside.bracket fun
stargazer_src <- readLines("stargazer/R/stargazer-internal.R")
# Move the length check 5 lines up so it precedes is.na(.)
stargazer_src[1990] <- stargazer_src[1995]
stargazer_src[1995] <- ""
# Save back
writeLines(stargazer_src, con="stargazer/R/stargazer-internal.R")
# Compile and install the patched package
install.packages("stargazer", repos = NULL, type="source")
load("stargazer")


# Load Datasets for Analysis -------------------------------------------
vib.thresh <- read.csv("Database_Ephys - Vibration.csv", header = TRUE, skip = 0, na.strings = "NA")
devo.data <- read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA")
morph.data.juv <- read.csv("Database_Morphometrics - Froglet_Toadlet Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")
hear.thresh <- read.csv("Database_Ephys - Hearing.csv", header = TRUE, skip = 0, na.strings = "NA")

# Clean/Join Datasets for Analysis -------------------------------------------

devo.data <- devo.data %>%
  mutate(unique.id.juv = factor(paste(gs.code, clutch, juv.tank.id, sep = "_")))

vib.thresh <- vib.thresh %>%
  mutate_at(vars(freq.hz, thresh.db, clip.thresh.db.plus3), as.numeric) %>% #change columns to numeric
  mutate_at(vars(gs, treatment), factor) %>% #change columns to factor
  #create new columns for different levels of identification
  mutate(combined.id = factor(paste(gs.code, tank.id, animal.id.num, sep = "_")), 
         unique.id.juv = factor(paste(gs.code, clutch, tank.id, sep = "_")),
         life.stage.num.sampling = factor(paste(life.stage, num.sampling, sep = "_")),
         id.life.stage.num.sampling = factor(paste(combined.id, life.stage, num.sampling, sep = "_")),
         # convert vibration threshold to m/s2 [[acceleration (m/s2) = (10^((rawdB-120)/20))*9.81]]
         thresh.ms2 = (10^((thresh.db-120)/20))*9.81,
         #thresh.cms2 = 100*((10^((thresh.db-120)/20))*9.81),
         #thresh.mms2 = 1000*((10^((thresh.db-120)/20))*9.81),
         rel.dB = (thresh.db-120),
         clip.thresh.db.plus3.ms2 = (10^((clip.thresh.db.plus3-120)/20))*9.81,
         # assign >dB "thresholds" where I couldn't get OVER threshold due to clipping
         clip.thresh.db.plus3 = if_else(is.na(clip.thresh.db.plus3) == TRUE & is.na(thresh.db) == FALSE, 
                                               thresh.db,
                                               clip.thresh.db.plus3)
         ) %>%
  #add mean developmental time for each juvenile tank, leaving NAs for adult tanks
  left_join(devo.data %>% 
              dplyr::group_by(unique.id.juv) %>%
              summarize(mean.days.forelimb = mean(days.forelimb, na.rm = TRUE)),
            join_by(unique.id.juv)
  ) %>%
  #remove individuals not included in final dataset
  dplyr::filter(date != "2025-04-24", #weird data day for RM_C_J029_6 
    date != "2025-04-14",
    !(date == "2025-05-06" & combined.id != "RM_J041_8"),
    !(date == "2025-05-12" & combined.id != "RM_J041_8"),
    !(date == "2025-05-05" & combined.id == "RM_J041_9"),
    combined.id != "RM_J027_4", #overflow individual
    is.na(svl.mm) == FALSE #adults that I started but didn't finish
  ) 


hear.thresh <- hear.thresh %>%
  mutate_at(vars(freq.hz, thresh.db, clip.thresh.db.plus3), as.numeric) %>% #change columns to numeric
  mutate_at(vars(gs, treatment), factor) %>% #change columns to factor
  #create new columns for different levels of identification
  mutate(combined.id = factor(paste(gs.code, tank.id, animal.id.num, sep = "_")),
         unique.id.juv = factor(paste(gs.code, clutch, tank.id, sep = "_")),
         life.stage.num.sampling = factor(paste(life.stage, num.sampling, sep = "_")),
         id.life.stage.num.sampling = factor(paste(combined.id, life.stage, num.sampling, sep = "_")),
         clip.thresh.db.plus3 = if_else(is.na(clip.thresh.db.plus3) == TRUE & is.na(thresh.db) == FALSE, 
                                        thresh.db,
                                        clip.thresh.db.plus3)
  ) %>%
  #add mean developmental time for each juvenile tank, leaving NAs for adult tanks
  left_join(devo.data %>% 
              dplyr::group_by(unique.id.juv) %>%
              summarize(mean.days.forelimb = mean(days.forelimb, na.rm = TRUE)),
            join_by(unique.id.juv)
            ) %>%
  #remove individuals not included in final dataset
  dplyr::filter(date != "2025-04-24", #weird data day for RM_C_J029_6
    freq.hz != 400,
    !(date == "2025-05-06" & combined.id != "RM_J041_8"),
    !(date == "2025-05-12" & combined.id != "RM_J041_8"),
    !(date == "2025-05-05" & combined.id == "RM_J041_9"),
    combined.id != "RM_J027_4") #overflow individual


morph.data.juv <- morph.data.juv %>%
  #filter out data that does not meet criteria for analysis
  filter(gs.code == "RM" & 
           (clutch == "B" | clutch == "C" | clutch == "D" | clutch == "E") &
           data.type == "ephys" &
           post.mm.weeks != "8-10") %>%
  #create new columns for different levels of identification
  mutate(combined.id = factor(paste(gs.code, juv.tank.id, animal.id.num, sep = "_")),
         unique.id = factor(paste(gs.code, clutch, juv.tank.id, sep = "_")),
         combined.id = factor(paste(gs.code, juv.tank.id, animal.id.num, sep = "_")),
         id.life.stage.num.sampling = factor(paste(gs.code, juv.tank.id, animal.id.num, life.stage, post.mm.sampling, sep = "_")),
         # create column so can plot post.mm.weeks on numeric scale
         post.mm.weeks.num = case_when(post.mm.weeks == "11-13" ~ mean(c(11,13)),
                                       post.mm.weeks == "23-25" ~ mean(c(23,25)),
                                       post.mm.weeks == "51-53" ~ mean(c(51,53)))
         ) %>%
  #add mean developmental time for each juvenile tank, leaving NAs for adult tanks
  left_join(devo.data %>% 
              dplyr::group_by(unique.id.juv) %>%
              summarize(mean.days.forelimb = mean(days.forelimb, na.rm = TRUE)),
            join_by(unique.id == unique.id.juv)
  )

# Create juvenile-only datasets ----------
vib.thresh.juv = vib.thresh %>% filter(life.stage == "juvenile")
hear.thresh.juv = hear.thresh %>% filter(life.stage == "juvenile")

# create clean datasets without NAs to use for dredge
vib.thresh.clean = vib.thresh %>% filter(is.na(thresh.db) == FALSE)
hear.thresh.clean = hear.thresh %>% filter(is.na(thresh.db) == FALSE)


# Summarize datasets-------------------
vib.thresh.summary <- vib.thresh %>%
  group_by(life.stage, num.sampling, freq.hz, .drop = TRUE) %>%
  summarise(n = n(), 
            mean.thresh.ms2 = mean(thresh.ms2, na.rm= TRUE), 
            sd.thresh.ms2 = sd(thresh.ms2, na.rm= TRUE),
            mean.mass.g = mean(mass.g, na.rm = TRUE),
            mean.svl.mm = mean(svl.mm, na.rm = TRUE))

hear.thresh.summary <- hear.thresh %>%
  group_by(life.stage, num.sampling, freq.hz, .drop = TRUE) %>%
  summarise(n = n(), 
            mean.thresh.db = mean(thresh.db, na.rm= TRUE), 
            sd.thresh.db = sd(thresh.db, na.rm= TRUE),
            mean.mass.g = mean(mass.g, na.rm = TRUE),
            mean.svl.mm = mean(svl.mm, na.rm = TRUE))


# Create function for juvenile-only model comparison -----------------------------------------
thresh_model_compare_juv <- function(x){
  #candidate models - manually defined
  #does the effect of age (num.sampling) depend on development time, size, and frequency?
  lmm.full.1 <- lmer(thresh.db ~ svl.mm*mean.days.forelimb*factor(freq.hz)*factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age depend on the effects of development time and size similarly across frequency?
  lmm.3intx.1 <- lmer(thresh.db ~ svl.mm*mean.days.forelimb*factor(num.sampling) + factor(freq.hz) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age function independently from the effects of size, frequency, and age?
  lmm.3intx.2 <- lmer(thresh.db ~ mean.days.forelimb*svl.mm*factor(freq.hz) + factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age depend on development time and frequency?
  lmm.3intx.3 <- lmer(thresh.db ~ mean.days.forelimb*factor(freq.hz)*factor(num.sampling) + svl.mm + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age depend on development time and frequency?
  lmm.3intx.4 <- lmer(thresh.db ~ mean.days.forelimb + svl.mm*factor(freq.hz)*factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age depend on frequency and the effects of size depend on tdevelopment time?
  lmm.2intx.1 <- lmer(thresh.db ~ svl.mm*mean.days.forelimb + factor(freq.hz)*factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of size depend on frequency and the effects of development time depend on age?
  lmm.2intx.2 <- lmer(thresh.db ~ svl.mm*factor(freq.hz) + mean.days.forelimb*factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of size depend on age and the effects of development time depend on frequency?
  lmm.2intx.3 <- lmer(thresh.db ~ svl.mm*factor(num.sampling) + mean.days.forelimb*factor(freq.hz) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of size depend on development time?
  lmm.1intx.1 <- lmer(thresh.db ~ svl.mm*mean.days.forelimb + factor(freq.hz) + factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age depend on development time?
  lmm.1intx.2 <- lmer(thresh.db ~ factor(num.sampling)*mean.days.forelimb + factor(freq.hz) + svl.mm + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of frequency depend on development time?
  lmm.1intx.3 <- lmer(thresh.db ~ factor(freq.hz)*mean.days.forelimb + factor(num.sampling) + svl.mm + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of size depend on age?
  lmm.1intx.4 <- lmer(thresh.db ~ svl.mm*factor(num.sampling) + factor(freq.hz) + mean.days.forelimb + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of size depend on frequency?
  lmm.1intx.5 <- lmer(thresh.db ~ svl.mm*factor(freq.hz) + factor(num.sampling) + mean.days.forelimb + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age depend on frequency?
  lmm.1intx.6 <- lmer(thresh.db ~ svl.mm + factor(num.sampling)*factor(freq.hz) + mean.days.forelimb + (1|combined.id), data = x, na.action = na.omit)
  
  #additive effects of age, development time, and frequency?
  lmm.add <- lmer(thresh.db ~ svl.mm + mean.days.forelimb + factor(freq.hz) + factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit, REML = TRUE)
  
  #null model
  lmm.null <- lmer(thresh.db ~ (1|combined.id), data = x, na.action = na.omit)
  
  #model comparison using AICc
  model.sel = arrange(AICc(lmm.full.1,
                           lmm.3intx.1, lmm.3intx.2, lmm.3intx.3, lmm.3intx.4,
                           lmm.2intx.1, lmm.2intx.2, lmm.2intx.3,
                           lmm.1intx.1, lmm.1intx.2, lmm.1intx.3, lmm.1intx.4, lmm.1intx.5, lmm.1intx.6,
                           lmm.add,
                           lmm.null), AICc) ; print(model.sel)
  final.mod = eval(parse(text = paste(rownames(model.sel)[1]))) #best supported model
  
  # check assumptions of best-fit model
  simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T)
  testDispersion(final.mod)
  testZeroInflation(final.mod)
  testCategorical(final.mod, catPred = x$freq.hz[is.na(x$thresh.db)==FALSE]) 
  testCategorical(final.mod, catPred = x$num.sampling[is.na(x$thresh.db)==FALSE]) 
  
  # estimates from best-supported model
  car::Anova(final.mod, type = "III")
  summary(final.mod)
  
  #export model comparison table to file
  stargazer::stargazer(lmm.full.1,
                       lmm.3intx.1, lmm.3intx.2, lmm.3intx.3, lmm.3intx.4,
                       lmm.2intx.1, lmm.2intx.2, lmm.2intx.3,
                       lmm.1intx.1, lmm.1intx.2, lmm.1intx.3, lmm.1intx.4, lmm.1intx.5, lmm.1intx.6,
                       lmm.add,
                       lmm.null,
                       type = "html", out=paste("modelcomparison", "_", deparse(substitute(x)), ".doc", sep=""), intercept.bottom = F, intercept.top = T, digits = 2)
  
}

# Supp. Table 1 - model comparison -----------------------------------------
thresh_model_compare_juv(vib.thresh.juv)

# Supp. Table 2 - model comparison -----------------------------------------
thresh_model_compare_juv(hear.thresh.juv)












# Analyze: Effect of life stage (adult vs. juvenile three timepoints) and size on VIB threshold ----

lmm.full.1 <- lmer(thresh.db ~ svl.mm*factor(freq.hz)*factor(life.stage.num.sampling) + (1|combined.id), data = vib.thresh.clean, na.action = na.fail)
mod_select <- dredge(lmm.full.1)

nrow(mod_select)

mod_select[1,]

final.mod = summary(get.models(mod_select, delta==0)[[1]])

par(mar = c(3,5,6,4))
plot(mod_select, labAsExpr = TRUE)

model.avg(mod_select, subset = delta < 4)

confset.95p <- get.models(mod_select, cumsum(weight) <= .95)
avgmod.95p <- model.avg(confset.95p)
summary(avgmod.95p)
confint(avgmod.95p)

final.mod$call

final.mod = lmer(thresh.db ~ factor(freq.hz) + factor(life.stage.num.sampling) + 
                   (1 | combined.id) + factor(freq.hz):factor(life.stage.num.sampling), 
                 data = vib.thresh.clean, na.action = na.fail)

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T)
testDispersion(final.mod)
testZeroInflation(final.mod)
testCategorical(final.mod, catPred = vib.thresh.juv.clean$freq.hz) 

car::Anova(final.mod, type = "III")
summary(final.mod)

pairs(emmeans(final.mod, ~ freq.hz), type = "response") 


#candidate models

#does the effect of life stage depend on the effects of size and frequency?
lmm.full.1 <- lmer(thresh.ms2 ~ svl.mm*factor(freq.hz)*life.stage.num.sampling + (1|combined.id), data = vib.thresh, na.action = na.omit)

#does the effect of life stage depend on the effects of frequency?
lmm.1intx.1 <- lmer(thresh.ms2 ~ svl.mm + factor(freq.hz)*life.stage.num.sampling + (1|combined.id), data = vib.thresh, na.action = na.omit)

#does the effect of life stage depend on the effects of size?
lmm.1intx.2 <- lmer(thresh.ms2 ~ svl.mm*life.stage.num.sampling + factor(freq.hz) + (1|combined.id), data = vib.thresh, na.action = na.omit)

#does the effect of size depend on the effects of frequency?
lmm.1intx.3 <- lmer(thresh.ms2 ~ svl.mm*factor(freq.hz) + life.stage.num.sampling + (1|combined.id), data = vib.thresh, na.action = na.omit)

#additive effects of life stage, size, and frequency?
lmm.add <- lmer(thresh.ms2 ~ svl.mm + factor(freq.hz) + life.stage.num.sampling + (1|combined.id), data = vib.thresh, na.action = na.omit)

#null model
lmm.null <- lmer(thresh.ms2 ~ (1|combined.id), data = vib.thresh, na.action = na.omit)

model.sel = arrange(AICc(lmm.full.1, 
                         lmm.1intx.1, lmm.1intx.2, lmm.1intx.3,
                         lmm.add,
                         lmm.null), AICc)
model.sel
final.mod = eval(parse(text = paste(rownames(model.sel)[1])))

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T)
testDispersion(final.mod)
testZeroInflation(final.mod)
testCategorical(final.mod, catPred = vib.thresh$freq.hz[is.na(vib.thresh$thresh.ms2)==FALSE]) 
testCategorical(final.mod, catPred = vib.thresh$life.stage[is.na(vib.thresh$thresh.ms2)==FALSE]) 

#log-transform response variable
lmm.add.log <- lmer(log(thresh.ms2) ~ svl.mm + factor(freq.hz) + life.stage.num.sampling + (1|combined.id), data = vib.thresh, na.action = na.omit)

# check assumptions of best-fit model with nlog-transformed response variable
simulateResiduals(fittedModel = lmm.add.log, quantreg=T, plot = T)
testDispersion(lmm.add.log)
testZeroInflation(lmm.add.log)
testCategorical(lmm.add.log, catPred = vib.thresh$freq.hz[is.na(vib.thresh$thresh.ms2)==FALSE]) 
testCategorical(lmm.add.log, catPred = vib.thresh$life.stage.num.sampling[is.na(vib.thresh$thresh.ms2)==FALSE]) 

# estimates from best-supported model
car::Anova(lmm.add.log, type = "II")
summary(lmm.add.log)

pairs(emmeans(lmm.add.log, ~ life.stage.num.sampling, by = "freq.hz"), type = "response") #backtransformed to the response scale




# Analyze: Effect of life stage and size on HEAR threshold ----

lmm.full.1 <- lmer(thresh.db ~ svl.mm*factor(freq.hz)*factor(life.stage.num.sampling) + (1|combined.id), data = hear.thresh.clean, na.action = na.fail)
mod_select <- dredge(lmm.full.1)

nrow(mod_select)

mod_select[1,]

final.mod = summary(get.models(mod_select, delta==0)[[1]])

par(mar = c(3,5,6,4))
plot(mod_select, labAsExpr = TRUE)

model.avg(mod_select, subset = delta < 4)

confset.95p <- get.models(mod_select, cumsum(weight) <= .95)
avgmod.95p <- model.avg(confset.95p)
summary(avgmod.95p)
confint(avgmod.95p)

final.mod$call

final.mod = lmer(thresh.db ~ factor(freq.hz) + factor(life.stage.num.sampling) + 
                   svl.mm + (1 | combined.id) + factor(freq.hz):factor(life.stage.num.sampling) + 
                   factor(freq.hz):svl.mm + factor(life.stage.num.sampling):svl.mm + 
                   factor(freq.hz):factor(life.stage.num.sampling):svl.mm, data = hear.thresh.clean, 
                 na.action = na.fail)
 

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T)
testDispersion(final.mod)
testZeroInflation(final.mod)
testCategorical(final.mod, catPred = hear.thresh.clean$life.stage.num.sampling) 
testCategorical(final.mod, catPred = hear.thresh.clean$freq.hz) 

car::Anova(final.mod, type = "III")
summary(final.mod)

pairs(emmeans(final.mod, ~ freq.hz), type = "response") #backtransformed to the response scale





#candidate models

#does the effect of life stage depend on the effects of size and frequency?
lmm.full.1 <- lmer(thresh.db ~ svl.mm*factor(freq.hz)*life.stage.num.sampling + (1|combined.id), data = hear.thresh, na.action = na.omit)

#does the effect of life stage depend on the effects of frequency?
lmm.1intx.1 <- lmer(thresh.db ~ svl.mm + factor(freq.hz)*life.stage.num.sampling + (1|combined.id), data = hear.thresh, na.action = na.omit)

#does the effect of life stage depend on the effects of size?
lmm.1intx.2 <- lmer(thresh.db ~ svl.mm*life.stage.num.sampling + factor(freq.hz) + (1|combined.id), data = hear.thresh, na.action = na.omit)

#does the effect of size depend on the effects of frequency?
lmm.1intx.3 <- lmer(thresh.db ~ svl.mm*factor(freq.hz) + life.stage.num.sampling + (1|combined.id), data = hear.thresh, na.action = na.omit)

#additive effects of life stage, size, and frequency?
lmm.add <- lmer(thresh.db ~ svl.mm + factor(freq.hz) + life.stage.num.sampling + (1|combined.id), data = hear.thresh, na.action = na.omit)

#null model
lmm.null <- lmer(thresh.db ~ (1|combined.id), data = hear.thresh, na.action = na.omit)

model.sel = arrange(AICc(lmm.full.1, 
                         lmm.1intx.1, lmm.1intx.2, lmm.1intx.3,
                         lmm.add,
                         lmm.null), AICc)
model.sel
final.mod = eval(parse(text = paste(rownames(model.sel)[1])))

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T)
testDispersion(final.mod)
testZeroInflation(final.mod)
testCategorical(final.mod, catPred = hear.thresh$freq.hz[is.na(hear.thresh$thresh.db)==FALSE]) 
testCategorical(final.mod, catPred = hear.thresh$life.stage.num.sampling[is.na(hear.thresh$thresh.db)==FALSE]) 

#log-transform response variable
lmm.full.1.log <- lmer(log(thresh.db) ~ svl.mm*factor(freq.hz)*life.stage.num.sampling + (1|combined.id), data = hear.thresh, na.action = na.omit)

# check assumptions of best-fit model with nlog-transformed response variable
simulateResiduals(fittedModel = lmm.full.1.log, quantreg=T, plot = T)
testDispersion(lmm.full.1.log)
testZeroInflation(lmm.full.1.log)
testCategorical(lmm.full.1.log, catPred = hear.thresh$freq.hz[is.na(hear.thresh$thresh.db)==FALSE]) 
testCategorical(lmm.full.1.log, catPred = hear.thresh$life.stage.num.sampling[is.na(hear.thresh$thresh.db)==FALSE]) 

# estimates from best-supported model
car::Anova(lmm.full.1, type = "III")
summary(lmm.full.1)

joint_tests(lmm.full.1, by = "freq.hz")
emtrends(lmm.full.1, pairwise ~ life.stage.num.sampling, var = "svl.mm", by = "freq.hz", type = "response")

emmip(lmm.full.1, life.stage.num.sampling ~ svl.mm | freq.hz, mult.name = "variety", cov.reduce = FALSE)



# Figure 1: JUV ONLY threshold by life stage and age within life stage (num.sampling) -------

fig.1a <- ggplot() +
  
  facet_grid(cols = vars(factor(num.sampling)),
             labeller = as_labeller(c('1' = "3 months", '2' = "6 months", '3' = "12 months"))
  )+
  
  geom_point(data = vib.thresh.juv,
             aes(y=clip.thresh.db.plus3.ms2, x = as.numeric(freq.hz), color = factor(num.sampling), group = combined.id),
             size = 3, alpha=0.3, show.legend = FALSE) +
  
  geom_line(data = vib.thresh.juv,
            aes(y=clip.thresh.db.plus3.ms2, x = as.numeric(freq.hz), color = factor(num.sampling), group = combined.id),
            size = 1, alpha=0.3, show.legend = FALSE) +
  
  #add points where threshold above clipping level
  geom_point(data = vib.thresh.juv %>%
               filter(is.na(thresh.ms2) == TRUE), #need to convert units within this line
             fill = "white", pch = 21,
             aes(y=clip.thresh.db.plus3.ms2, x = as.numeric(freq.hz), color = factor(num.sampling), group = combined.id),
             size = 3, stroke = 1, alpha=1.0, show.legend = FALSE) +
  
  # geom_line(data = vib.thresh,
  #           aes(y=thresh.ms2, x = as.numeric(freq.hz), color = life.stage.num.sampling, group = combined.id)) +
  
  stat_summary(data = vib.thresh.juv,
               fun.y=mean, geom="line", size = 1.2, colour="black",
               aes(y=thresh.ms2, x = as.numeric(freq.hz), group = factor(num.sampling)), show.legend = FALSE) +
  
  stat_summary(data = vib.thresh.juv,
               fun.data=function(x){mean_cl_normal(x, conf.int=.683)}, geom="errorbar", 
               width=0.1, size = 0.9, colour="black", alpha=1, 
               aes(y=thresh.ms2, x = as.numeric(freq.hz), group = factor(num.sampling)), show.legend = FALSE) +
  
  stat_summary(data = vib.thresh.juv,
               fun.y=mean, geom="point", colour = "black", pch=21, size=9, 
               aes(y=thresh.ms2, x = as.numeric(freq.hz), fill = factor(num.sampling), group = factor(num.sampling)), show.legend = FALSE) +
  
  scale_colour_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  scale_fill_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=18, color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12, color = "white", face = "bold")) +
  scale_x_continuous(name = "frequency (kHz)", limits = c(100,2000), breaks = c(sort(unique(vib.thresh$freq.hz)), 2000), labels = c("0.1","0.2", "0.3", "0.4", "0.7", "1.1", "1.5", "2.0")) +
  scale_y_continuous(name = "vibration threshold\n(dB re 1 m/s-2)")

# add num.sampling colors to facet grid (code from: https://github.com/tidyverse/ggplot2/issues/2096)
g.1a <- ggplot_gtable(ggplot_build(fig.1a))
strip_t <- which(grepl('strip-t', g.1a$layout$name))
fills <- c("#CC66FF", "#660066", "#330066")
k <- 1
for(i in strip_t){
  j <- which(grepl('rect', g.1a$grobs[[i]]$grobs[[1]]$childrenOrder))
  g.1a$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g.1a) #produces final plot with colored headers


fig.1a.reldB <- ggplot() +
  
  facet_grid(cols = vars(factor(num.sampling)),
             labeller = as_labeller(c('1' = "3 months", '2' = "6 months", '3' = "12 months"))
  )+
  
  geom_point(data = vib.thresh.juv,
             aes(y=rel.dB, x = as.numeric(freq.hz), color = factor(num.sampling), group = combined.id),
             size = 3, alpha=0.3, show.legend = FALSE) +
  
  geom_line(data = vib.thresh.juv,
            aes(y=rel.dB, x = as.numeric(freq.hz), color = factor(num.sampling), group = combined.id),
            size = 1, alpha=0.3, show.legend = FALSE) +
  
  #add points where threshold above clipping level
  geom_point(data = vib.thresh.juv %>%
                filter(is.na(thresh.db) == TRUE), 
              fill = "white", pch = 21,
              aes(y=clip.thresh.db.plus3-120, x = as.numeric(freq.hz), color = factor(num.sampling), group = combined.id), #need to convert units within this line
              size = 3, stroke = 1, alpha=1.0, show.legend = FALSE) +

  geom_line(data = vib.thresh.juv %>%
               filter(is.na(thresh.db) == TRUE), 
             aes(y=clip.thresh.db.plus3-120, x = as.numeric(freq.hz), color = factor(num.sampling), group = combined.id), #need to convert units within this line
             size = 3, alpha=0.3, show.legend = FALSE) +
  
  stat_summary(data = vib.thresh.juv,
               fun.y=mean, geom="line", size = 1.2, colour="black",
               aes(y=rel.dB, x = as.numeric(freq.hz), group = factor(num.sampling)), show.legend = FALSE) +
  
  stat_summary(data = vib.thresh.juv,
               fun.data=function(x){mean_cl_normal(x, conf.int=.683)}, geom="errorbar", 
               width=0.1, size = 0.9, colour="black", alpha=1, 
               aes(y=rel.dB, x = as.numeric(freq.hz), group = factor(num.sampling)), show.legend = FALSE) +
  
  stat_summary(data = vib.thresh.juv,
               fun.y=mean, geom="point", colour = "black", pch=21, size=9, 
               aes(y=rel.dB, x = as.numeric(freq.hz), fill = factor(num.sampling), group = factor(num.sampling)), show.legend = FALSE) +
  
  scale_colour_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  scale_fill_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=18, color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12, color = "white", face = "bold")) +
  scale_x_continuous(name = "frequency (kHz)", limits = c(100,2000), breaks = c(sort(unique(vib.thresh$freq.hz)), 2000), labels = c("0.1","0.2", "0.3", "0.4", "0.7", "1.1", "1.5", "2.0")) +
  scale_y_continuous(name = "vibration threshold\n(dB re hearing)")

# add num.sampling colors to facet grid (code from: https://github.com/tidyverse/ggplot2/issues/2096)
g.1a.reldB <- ggplot_gtable(ggplot_build(fig.1a.reldB))
strip_t <- which(grepl('strip-t', g.1a.reldB$layout$name))
fills <- c("#CC66FF", "#660066", "#330066")
k <- 1
for(i in strip_t){
  j <- which(grepl('rect', g.1a.reldB$grobs[[i]]$grobs[[1]]$childrenOrder))
  g.1a.reldB$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g.1a.reldB) #produces final plot with colored headers


fig.1b <- ggplot() +
  
  facet_grid(cols = vars(factor(num.sampling)),
             labeller = as_labeller(c('1' = "3 months", '2' = "6 months", '3' = "12 months"))
  )+
  
  geom_point(data = hear.thresh.juv,
             aes(y=clip.thresh.db.plus3, x = as.numeric(freq.hz), color = factor(num.sampling), group = combined.id),
             size = 3, alpha=0.3, show.legend = FALSE) +
  
  geom_line(data = hear.thresh.juv,
            aes(y=clip.thresh.db.plus3, x = as.numeric(freq.hz), color = factor(num.sampling), group = combined.id),
            size = 1, alpha=0.3, show.legend = FALSE) +
  
  #add points where threshold above clipping level
  geom_point(data = hear.thresh.juv %>%
               filter(is.na(thresh.db) == TRUE), #need to convert units within this line
             fill = "white", pch = 21,
             aes(y=clip.thresh.db.plus3, x = as.numeric(freq.hz), color = factor(num.sampling), group = combined.id),
             size = 3, stroke = 1, alpha=1.0, width = 50, height = 0, show.legend = FALSE) +

  stat_summary(data = hear.thresh.juv,
               fun.y=mean, geom="line", size = 1.2, colour="black",
               aes(y=thresh.db, x = as.numeric(freq.hz), group = factor(num.sampling)), show.legend = FALSE) +
  
  stat_summary(data = hear.thresh.juv,
               fun.data=function(x){mean_cl_normal(x, conf.int=.683)}, geom="errorbar", 
               width=0.1, size = 0.9, colour="black", alpha=1, 
               aes(y=thresh.db, x = as.numeric(freq.hz), group = factor(num.sampling)), show.legend = FALSE) +
  
  stat_summary(data = hear.thresh.juv,
               fun.y=mean, geom="point", colour = "black", pch=21, size=9, 
               aes(y=thresh.db, x = as.numeric(freq.hz), fill = factor(num.sampling), group = factor(num.sampling)), show.legend = FALSE) +
  
  scale_colour_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  scale_fill_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=18, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y=element_text(size=18, color = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12, color = "white", face = "bold")) +
  scale_x_continuous(name = "frequency (kHz)", limits = c(100,2000), breaks = c(100,200,sort(unique(hear.thresh.juv$freq.hz))), labels = c("0.1","0.2", "0.3", "0.4", "0.7", "1.1", "1.5", "2.0")) +
  scale_y_continuous(name = "hearing threshold\n(dB)")

# add num.sampling colors to facet grid (code from: https://github.com/tidyverse/ggplot2/issues/2096)
g.1b <- ggplot_gtable(ggplot_build(fig.1b))
strip_t <- which(grepl('strip-t', g.1b$layout$name))
fills <- c("#CC66FF", "#660066", "#330066")
k <- 1
for(i in strip_t){
  j <- which(grepl('rect', g.1b$grobs[[i]]$grobs[[1]]$childrenOrder))
  g.1b$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g.1b) #produces final plot with colored headers


png("~/Desktop/R Working Directory/Plots/Figure1.png", units = "in", res = 300, width = 16, height = 12)
ggarrange(g.1a,g.1b,
          ncol = 1, nrow = 2)
dev.off()


png("~/Desktop/R Working Directory/Plots/Figure1.png", units = "in", res = 300, width = 16, height = 12)
ggarrange(g.1a.reldB,g.1b,
          ncol = 1, nrow = 2,
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))
dev.off()

# FIGURE 2:  heatmap: threshold by size and age (num.sampling) -------
#helpful code to make panel plots: https://oscarperpinan.github.io/rastervis/FAQ.html

fig.2a.3 <- levelplot(rel.dB ~ svl.mm * freq.hz, 
                    data = vib.thresh.juv %>% filter(
                      num.sampling == 3),
                    panel = panel.levelplot.points, cex = 1.2,
                    col.regions=colorRampPalette(brewer.pal(9, "Blues"))(25),
                    colorkey=list(labels=list(cex=1.5), title = "vibrational threshold (dB re hearing)"),
                    xlab=list(label = 'snout-vent length (mm)', cex=1.5),
                    ylab=list(label = 'frequency (kHz)', cex = 1.5),
                    ylim = c(50,1550),
                    #auto.key = list(title = "vib threshold", cex.title = 1.2),
                    scales=list(y=list(at=c(100, 200, 300, 400, 700, 1100, 1500), 
                                       labels = c("0.1", "0.2", "0.3", "0.4", "0.7", "1.1", "1.5"),
                                       cex=1),
                                x=list(cex=1)
                    ),
                    main="12 months"
) +
  layer_(panel.2dsmoother(..., n=200))
#grid.text("vib threshold", x=0.92, y=0.2, hjust=0.5, vjust=1)

png("~/Desktop/R Working Directory/Plots/Figure2a.png", units = "in", res = 300, width = 24, height = 12)
gridExtra::grid.arrange(fig.2a.1, fig.2a.2, fig.2a.3,
                        ncol = 3, nrow = 1)
dev.off()


fig.2b.3 <- levelplot(thresh.db ~ svl.mm * freq.hz, 
                    data = hear.thresh.juv %>% filter(num.sampling == 3),
                    panel = panel.levelplot.points, cex = 1.2,
                    col.regions=colorRampPalette(brewer.pal(9, "Blues"))(25),
                    colorkey=list(labels=list(cex=1.5), title = "hearing threshold (dB)"),
                    xlab=list(label = 'snout-vent length (mm)', cex=1.5),
                    ylab=list(label = 'frequency (kHz)', cex = 1.5),
                    ylim = c(250,2050),
                    scales=list(y=list(at=c(300, 500, 700, 1100, 1500, 2000), 
                                       labels = c("0.3", "0.5", "0.7", "1.1", "1.5", "2.0"),
                                       cex=1),
                                x = list(cex=1)
                    ),
                    main="12 months"
) +
  layer_(panel.2dsmoother(..., n=200))

png("~/Desktop/R Working Directory/Plots/Figure2b.png", units = "in", res = 300, width = 24, height = 12)
gridExtra::grid.arrange(fig.2b.1, fig.2b.2, fig.2b.3,
                        ncol = 3, nrow = 1)
dev.off()


# FIGURE 2:  heatmap: threshold by development time and age (num.sampling) -------
#helpful code to make panel plots: https://oscarperpinan.github.io/rastervis/FAQ.html

fig.2c.3 <- levelplot(rel.dB ~ mean.days.forelimb * freq.hz, 
                    data = vib.thresh.juv %>% filter(
                      num.sampling == 3
                    ),
                    panel = panel.levelplot.points, cex = 1.2,
                    col.regions=colorRampPalette(brewer.pal(9, "Blues"))(25),
                    colorkey=list(labels=list(cex=1.5), title = "vibrational threshold (dB re hearing)"),
                    xlab=list(label = 'mean development time (days)', cex=1.5),
                    ylab=list(label = 'frequency (kHz)', cex = 1.5),
                    ylim = c(50,1550),
                    #auto.key = list(title = "vib threshold", cex.title = 1.2),
                    scales=list(y=list(at=c(100, 200, 300, 400, 700, 1100, 1500), 
                                       labels = c("0.1", "0.2", "0.3", "0.4", "0.7", "1.1", "1.5"),
                                       cex=1),
                                x=list(cex=1)
                    ),
                    main="12 months"
) +
  layer_(panel.2dsmoother(..., n=200))
#grid.text("vib threshold", x=0.92, y=0.2, hjust=0.5, vjust=1)

png("~/Desktop/R Working Directory/Plots/Figure2c.png", units = "in", res = 300, width = 24, height = 12)
gridExtra::grid.arrange(fig.2c.1, fig.2c.2, fig.2c.3,
                        ncol = 3, nrow = 1)
dev.off()



fig.2d.3 <- levelplot(thresh.db ~ mean.days.forelimb * freq.hz, 
                      data = hear.thresh.juv %>% filter(num.sampling == 3),
                    panel = panel.levelplot.points, cex = 1.2,
                    col.regions=colorRampPalette(brewer.pal(9, "Blues"))(25),
                    colorkey=list(labels=list(cex=1.5), title = "hearing threshold (dB)"),
                    xlab=list(label = 'mean development time (days)', cex=1.5),
                    ylab=list(label = 'frequency (kHz)', cex = 1.5),
                    ylim = c(250,2050),
                    scales=list(y=list(at=c(300, 500, 700, 1100, 1500, 2000), 
                                       labels = c("0.3", "0.5", "0.7", "1.1", "1.5", "2.0"),
                                       cex=1),
                                x = list(cex=1)
                    ),
                    main = "12 months"
) +
  layer_(panel.2dsmoother(..., n=200))


png("~/Desktop/R Working Directory/Plots/Figure2d.png", units = "in", res = 300, width = 24, height = 12)
gridExtra::grid.arrange(fig.2d.1, fig.2d.2, fig.2d.3,
                        ncol = 3, nrow = 1)
dev.off()


# FIGURE 2?? heatmap plot size by devo time with thresholds
for(i in 1:length(unique(hear.thresh.juv$freq.hz))){
  print(levelplot(thresh.db ~ mean.days.forelimb*svl.mm, data = hear.thresh.juv %>% filter(freq.hz == unique(hear.thresh.juv$freq.hz)[i]),
          panel = panel.levelplot.points, cex = 1.2,
          col.regions=colorRampPalette(brewer.pal(9, "Blues"))(25),
          colorkey=list(labels=list(cex=1.5), title = "hearing threshold (dB)"),
          xlab=list(label = 'mean development time (days)', cex=1.5),
          ylab=list(label = 'svl (mm)', cex = 1.5),
          #ylim = c(250,2050),
          # scales=list(y=list(at=c(300, 500, 700, 1100, 1500, 2000), 
          #                    labels = c("0.3", "0.5", "0.7", "1.1", "1.5", "2.0"),
          #                    cex=1),
                      # x = list(cex=1)
          #)
) +
  layer_(panel.2dsmoother(..., n=200))
)
}

# Figure 3: Juvenile-Adult comparison threshold by life stage and age within life stage (num.sampling) -------

vib.thresh.plot <- ggplot() +
  
  geom_jitter(data = vib.thresh,
             aes(y=thresh.ms2, x = as.numeric(freq.hz), color = life.stage.num.sampling, group = combined.id),
             size = 5, stroke = 1, alpha=0.3, width = 50, height = 0) +
  
  #add points where threshold above clipping level
  geom_jitter(data = vib.thresh %>%
                filter(is.na(thresh.ms2) == TRUE), #need to convert units within this line
              fill = "white", pch = 21,
              aes(y=clip.thresh.db.plus3.ms2, x = as.numeric(freq.hz), color = life.stage.num.sampling, group = combined.id),
              size = 5, stroke = 1, alpha=0.3, width = 50, height = 0) +
  
  # geom_line(data = vib.thresh,
  #           aes(y=thresh.ms2, x = as.numeric(freq.hz), color = life.stage.num.sampling, group = combined.id)) +

  stat_summary(data = vib.thresh,
               fun.y=mean, geom="line", size = 1.2, colour="black",
               aes(y=thresh.ms2, x = as.numeric(freq.hz), group = life.stage.num.sampling)) +
    
  stat_summary(data = vib.thresh,
               fun.data=function(x){mean_cl_normal(x, conf.int=.683)}, geom="errorbar", 
               width=0.1, size = 0.9, colour="black", alpha=1, 
               aes(y=thresh.ms2, x = as.numeric(freq.hz), group = life.stage.num.sampling)) +
  
  stat_summary(data = vib.thresh,
               fun.y=mean, geom="point", colour = "black", pch=21, size=9, 
               aes(y=thresh.ms2, x = as.numeric(freq.hz), fill = life.stage.num.sampling, group = life.stage.num.sampling)) +
  
  scale_colour_manual(values = c("#666666","#CC66FF", "#660066", "#330066"), labels = c("adult", "3 months", "6 months", "12 months")) +
  scale_fill_manual(values = c("#666666","#CC66FF", "#660066", "#330066"), labels = c("adult", "3 months", "6 months", "12 months")) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=24, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y=element_text(size=24, color = "black"),
        axis.title.y = element_text(size=28),
        axis.title.x = element_text(size=28),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "frequency (kHz)", limits = c(50,1550), breaks = sort(unique(vib.thresh$freq.hz)), labels = c("0.1","0.2", "0.3", "0.4", "0.7", "1.1", "1.5")) +
  scale_y_continuous(name = "vibration threshold\n(dB re 1 m/s-2)")


vib.thresh.plot.reldB <- ggplot() +
  
  geom_jitter(data = vib.thresh,
              aes(y=rel.dB, x = as.numeric(freq.hz), color = life.stage.num.sampling, group = combined.id),
              size = 5, stroke = 1, alpha=0.3, width = 50, height = 0, show.legend = FALSE) +
  
  #add points where threshold above clipping level
  geom_jitter(data = vib.thresh %>%
                filter(is.na(thresh.db) == TRUE), #need to convert units within this line
              fill = "white", pch = 21,
              aes(y=clip.thresh.db.plus3-120, x = as.numeric(freq.hz), color = life.stage.num.sampling, group = combined.id),
              size = 5, stroke = 1, alpha=1, show.legend = FALSE) +
  
  # geom_line(data = vib.thresh,
  #           aes(y=thresh.ms2, x = as.numeric(freq.hz), color = life.stage.num.sampling, group = combined.id)) +
  
  stat_summary(data = vib.thresh,
               fun.y=mean, geom="line", size = 1.2, 
               aes(y=rel.dB, x = as.numeric(freq.hz), color = life.stage.num.sampling, group = life.stage.num.sampling)) +
  
  stat_summary(data = vib.thresh,
               fun.data=function(x){mean_cl_normal(x, conf.int=.683)}, geom="errorbar", 
               width=0.1, size = 0.9, colour="black", alpha=1, 
               aes(y=rel.dB, x = as.numeric(freq.hz), group = life.stage.num.sampling)) +
  
  stat_summary(data = vib.thresh,
               fun.y=mean, geom="point", colour = "black", pch=21, size=9, 
               aes(y=rel.dB, x = as.numeric(freq.hz), fill = life.stage.num.sampling, group = life.stage.num.sampling)) +
  
  scale_colour_manual(values = c("#666666","#CC66FF", "#660066", "#330066"), labels = c("adult", "3 months", "6 months", "12 months")) +
  scale_fill_manual(values = c("#666666","#CC66FF", "#660066", "#330066"), labels = c("adult", "3 months", "6 months", "12 months")) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=24, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y=element_text(size=24, color = "black"),
        axis.title.y = element_text(size=28),
        axis.title.x = element_text(size=28),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "frequency (kHz)", limits = c(50,1550), breaks = sort(unique(vib.thresh$freq.hz)), labels = c("0.1","0.2", "0.3", "0.4", "0.7", "1.1", "1.5")) +
  scale_y_continuous(name = "vibration threshold\n(dB re hearing)")




hear.thresh.plot <- ggplot() +
  
  #add vertical lines at dominant frequency of breeding calls (from Yasumiba et al. Why do male and female cane toads, Rhinella marina, respond differently to advertisement calls?)
  geom_rect(data = hear.thresh, fill = "gray95",
            aes(xmin = 560.20, xmax = 765.59,   ymin = -Inf,
                ymax = Inf)) +
  
  geom_jitter(data = hear.thresh,
             aes(y=thresh.db, x = as.numeric(freq.hz), color = life.stage.num.sampling, group = combined.id),
             size = 5, alpha=0.3, width = 50, height = 2) +
  
  #add points where threshold above clipping level
  geom_jitter(data = hear.thresh %>%
               filter(is.na(thresh.db) == TRUE),
             fill = "white", pch = 21,
             aes(y=clip.thresh.db.plus3, x = as.numeric(freq.hz), color = life.stage.num.sampling, group = combined.id),
             size = 5, stroke = 1, alpha=1, width = 50, height = 2) +
  
  stat_summary(data = hear.thresh,
               fun.y=mean, geom="line", size = 1.2, 
               aes(y=thresh.db, x = as.numeric(freq.hz), group = life.stage.num.sampling, color = life.stage.num.sampling)) +
  
  stat_summary(data = hear.thresh,
               fun.data=function(x){mean_cl_normal(x, conf.int=.683)}, geom="errorbar", 
               width=0.1, size = 0.9, colour="black", alpha=1, 
               aes(y=thresh.db, x = as.numeric(freq.hz), group = life.stage.num.sampling)) +
  
  stat_summary(data = hear.thresh,
               fun.y=mean, geom="point", colour = "black", pch=21, size=9, 
               aes(y=thresh.db, x = as.numeric(freq.hz), fill = life.stage.num.sampling, group = life.stage.num.sampling)) +
  
  scale_colour_manual(values = c("#666666","#CC66FF", "#660066", "#330066"), labels = c("adult", "3 months", "6 months", "12 months")) +
  scale_fill_manual(values = c("#666666","#CC66FF", "#660066", "#330066"), labels = c("adult", "3 months", "6 months", "12 months")) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=24, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y=element_text(size=24, color = "black"),
        axis.title.y = element_text(size=28),
        axis.title.x = element_text(size=28),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "frequency (kHz)", limits = c(250,2050), breaks = sort(unique(hear.thresh$freq.hz)), labels = c("0.3", "0.5", "0.7", "1.1", "1.5", "2.0")) +
  scale_y_continuous(name = "hearing threshold (dB)")


png("~/Desktop/R Working Directory/Plots/Figure3.png", units = "in", res = 300, width = 16, height = 12)
ggarrange(vib.thresh.plot.reldB, hear.thresh.plot, 
          ncol = 1,
          nrow = 2,
          common.legend = TRUE,
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))
dev.off()




# PLOT: size as a function of age -------
ggplot() +
  geom_point(data = morph.data.juv,
             aes(x = post.mm.weeks.num, y = svl.mm)) +
  
  geom_line(data = morph.data.juv,
             aes(x = post.mm.weeks.num, y = svl.mm, group = combined.id)) +
  
  theme_bw() +
  scale_x_continuous(name = "age (weeks post metamorphosis") +
  scale_y_continuous(name = "snout-vent length (mm)")





# PLOT: threshold by development time ----
vib.thresh.plot.larvdur <- ggplot() +
  
  facet_grid(cols = vars(factor(freq.hz))) +

  geom_point(data = vib.thresh[vib.thresh$life.stage == "juvenile" & vib.thresh$treatment != "overflow",], pch = 21, color = "black",
             aes(y=thresh.ms2, x = mean.days.forelimb, fill = factor(num.sampling)), ## SMS YOU ARE HERE, color by freq,
             size = 4, alpha=1) +
  
  geom_smooth(method = "lm",
              data = vib.thresh[vib.thresh$life.stage == "juvenile" & vib.thresh$treatment != "overflow",],
              aes(y=thresh.ms2, x = mean.days.forelimb, color = factor(num.sampling), group = factor(num.sampling))) +
  
  scale_color_manual(values=c("gray90", "gray40", "gray20")) +
  scale_fill_manual(values=c("gray90", "gray40", "gray20")) +
  
  # geom_smooth(method = "lm",
  #             data = vib.thresh[vib.thresh$life.stage == "juvenile" & vib.thresh$treatment != "overflow",],
  #             color = "black",
  #             aes(y=thresh.ms2, x = mean.days.forelimb)) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=18, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=18, color = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "larval duration (days)", limits = c(38,70)) +
  scale_y_continuous(name = "vibrational threshold (dB re 1 m/s-2)")


hear.thresh.plot.larvdur <- ggplot() +
  
  facet_grid(cols = vars(factor(freq.hz))) +
  
  geom_point(data = hear.thresh[hear.thresh$life.stage == "juvenile" & hear.thresh$treatment != "overflow",], pch = 21, color = "black",
             aes(y=clip.thresh.db.plus3, x = mean.days.forelimb, fill = factor(num.sampling)), ## SMS YOU ARE HERE, color by freq,
             size = 4, alpha=1) +
  
  geom_smooth(method = "lm",
              data = hear.thresh[hear.thresh$life.stage == "juvenile" & hear.thresh$treatment != "overflow",],
              aes(y=clip.thresh.db.plus3, x = mean.days.forelimb, color = factor(num.sampling), group = factor(num.sampling))) +
  
  scale_color_manual(values=c("gray90", "gray40", "gray20")) +
  scale_fill_manual(values=c("gray90", "gray40", "gray20")) +
  
  # geom_smooth(method = "lm",
  #             data = vib.thresh[vib.thresh$life.stage == "juvenile" & vib.thresh$treatment != "overflow",],
  #             color = "black",
  #             aes(y=thresh.ms2, x = mean.days.forelimb)) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=18, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=18, color = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "larval duration (days)", limits = c(38,70)) +
  scale_y_continuous(name = "hearing threshold (dB)", limits = c(80,140), breaks = seq(80,140,15), labels = seq(80,140,15))

ggarrange(vib.thresh.plot.larvdur, hear.thresh.plot.larvdur,
          nrow = 2,
          common.legend = TRUE)





# PLOT: plotting by mass -------------------------------------------------------------

vib.thresh.plot.mass <- ggplot() +
  
  facet_grid(cols = vars(factor(freq.hz))) +
  
  geom_point(data = vib.thresh %>%
               filter(life.stage == "juvenile",
                      treatment != "overflow",
                      freq.hz != "NA",), 
             pch = 21, color = "black",
             aes(y=thresh.ms2, x = mass.g, fill = factor(num.sampling)), ## SMS YOU ARE HERE, color by freq,
             size = 4, alpha=1) +
  
  geom_smooth(method = "lm",
              data = vib.thresh %>%
                filter(life.stage == "juvenile",
                       treatment != "overflow",
                       freq.hz != "NA",),
              aes(y=thresh.ms2, x = mass.g, color = factor(num.sampling), group = factor(num.sampling))) +
  
  scale_colour_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  scale_fill_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  
  # geom_smooth(method = "lm",
  #             data = vib.thresh[vib.thresh$life.stage == "juvenile" & vib.thresh$treatment != "overflow",],
  #             color = "black",
  #             aes(y=thresh.ms2, x = mean.days.forelimb)) +
  
  theme_bw() +
  theme(strip.text = element_text(size = 20),
    legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=12, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=12, color = "black"),
        axis.title.y = element_text(size=28),
        axis.title.x = element_text(size=28),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "mass (g)", limits = c(0.3,4), breaks = seq(0.2,3.8,0.6), labels = seq(0.2,3.8,0.6)) +
  scale_y_continuous(name = "vib. threshold (dB re 1 m/s-2)")


hear.thresh.plot.mass <- ggplot() +
  
  facet_grid(cols = vars(factor(freq.hz))) +
  
  geom_point(data = hear.thresh %>%
               filter(life.stage == "juvenile",
                      treatment != "overflow",
                      freq.hz != "NA",),
             pch = 21, color = "black",
             aes(y=clip.thresh.db.plus3, x = mass.g, fill = factor(num.sampling)), ## SMS YOU ARE HERE, color by freq,
             size = 4, alpha=1) +
  
  #threshold above clipping level
  geom_point(data = hear.thresh %>%
               filter(life.stage == "juvenile",
                      treatment != "overflow",
                      is.na(thresh.db) == TRUE,
                      freq.hz != "NA",),
             pch = 21, fill = "white",
             aes(y=clip.thresh.db.plus3, x = mass.g, color = factor(num.sampling)), ## SMS YOU ARE HERE, color by freq,
             size = 4, alpha=1) +
  
  geom_smooth(method = "lm",
              data = hear.thresh %>%
                filter(life.stage == "juvenile",
                       treatment != "overflow",
                       is.na(thresh.db) == FALSE,
                       freq.hz != "NA",),
              aes(y=clip.thresh.db.plus3, x = mass.g, color = factor(num.sampling), group = factor(num.sampling))) +
  
  scale_colour_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  scale_fill_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  
  # geom_smooth(method = "lm",
  #             data = vib.thresh[vib.thresh$life.stage == "juvenile" & vib.thresh$treatment != "overflow",],
  #             color = "black",
  #             aes(y=thresh.ms2, x = mean.days.forelimb)) +
  
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=12, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=12, color = "black"),
        axis.title.y = element_text(size=28),
        axis.title.x = element_text(size=28),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "mass (g)", limits = c(0.3,4), breaks = seq(0.2,3.8,0.6), labels = seq(0.2,3.8,0.6)) +
  scale_y_continuous(name = "hearing threshold (dB)", limits = c(80,140), breaks = seq(80,140,15), labels = seq(80,140,15))


png("~/Desktop/R Working Directory/Plots/Juv-Mass.png", units = "in", res = 300, width = 16, height = 12)
ggarrange(vib.thresh.plot.mass, hear.thresh.plot.mass, 
          ncol = 1,
          nrow = 2,
          common.legend = TRUE,
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))
dev.off()


# PLOT: plotting by svl -------------------------------------------------------------

vib.thresh.plot.svl <- ggplot() +
  
  facet_grid(cols = vars(factor(freq.hz))) +
  
  geom_point(data = vib.thresh %>%
               filter(life.stage == "juvenile",
                      treatment != "overflow",
                      freq.hz != "NA",), 
             pch = 21, color = "black",
             aes(y=thresh.ms2, x = svl.mm, fill = factor(num.sampling)), ## SMS YOU ARE HERE, color by freq,
             size = 4, alpha=1) +
  
  geom_smooth(method = "lm",
              data = vib.thresh %>%
                filter(life.stage == "juvenile",
                       treatment != "overflow",
                       freq.hz != "NA",),
              aes(y=thresh.ms2, x = svl.mm, color = factor(num.sampling), group = factor(num.sampling))) +
  
  scale_colour_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  scale_fill_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  
  # geom_smooth(method = "lm",
  #             data = vib.thresh[vib.thresh$life.stage == "juvenile" & vib.thresh$treatment != "overflow",],
  #             color = "black",
  #             aes(y=thresh.ms2, x = mean.days.forelimb)) +
  
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=14, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=14, color = "black"),
        axis.title.y = element_text(size=28),
        axis.title.x = element_text(size=28),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "svl (mm)", limits = c(15,38), breaks = seq(15,38,5), labels = seq(15,38,5)) +
  scale_y_continuous(name = "vib. threshold (dB re 1 m/s-2)")


hear.thresh.plot.svl <- ggplot() +
  
  facet_grid(cols = vars(factor(freq.hz))) +
  
  geom_point(data = hear.thresh %>%
               filter(life.stage == "juvenile",
                      treatment != "overflow",
                      freq.hz != "NA",),
             pch = 21, color = "black",
             aes(y=clip.thresh.db.plus3, x = svl.mm, fill = factor(num.sampling)), ## SMS YOU ARE HERE, color by freq,
             size = 4, alpha=1) +
  
  #threshold above clipping level
  geom_point(data = hear.thresh %>%
               filter(life.stage == "juvenile",
                      treatment != "overflow",
                      is.na(thresh.db) == TRUE,
                      freq.hz != "NA",),
             pch = 21, fill = "white",
             aes(y=clip.thresh.db.plus3, x = svl.mm, color = factor(num.sampling)), ## SMS YOU ARE HERE, color by freq,
             size = 4, alpha=1) +
  
  geom_smooth(method = "lm",
              data = hear.thresh %>%
                filter(life.stage == "juvenile",
                       treatment != "overflow",
                       is.na(thresh.db) == FALSE,
                       freq.hz != "NA",),
              aes(y=clip.thresh.db.plus3, x = svl.mm, color = factor(num.sampling), group = factor(num.sampling))) +
  
  scale_colour_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  scale_fill_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  
  # geom_smooth(method = "lm",
  #             data = vib.thresh[vib.thresh$life.stage == "juvenile" & vib.thresh$treatment != "overflow",],
  #             color = "black",
  #             aes(y=thresh.ms2, x = mean.days.forelimb)) +
  
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=14, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=14, color = "black"),
        axis.title.y = element_text(size=28),
        axis.title.x = element_text(size=28),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "svl (mm)", limits = c(15,38), breaks = seq(15,38,5), labels = seq(15,38,5)) +
  scale_y_continuous(name = "hearing threshold (dB)", limits = c(80,140), breaks = seq(80,140,15), labels = seq(80,140,15))


png("~/Desktop/R Working Directory/Plots/Juv-SVL.png", units = "in", res = 300, width = 16, height = 12)
ggarrange(vib.thresh.plot.svl, hear.thresh.plot.svl, 
          ncol = 1,
          nrow = 2,
          common.legend = TRUE,
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))
dev.off()


# PLOT: plotting by individual svl -------------------------------------------------------------

vib.thresh.plot.svl.indiv <- ggplot() +
  
  facet_grid(cols = vars(factor(freq.hz))) +
  
  geom_point(data = vib.thresh %>%
               filter(life.stage == "juvenile",
                      treatment != "overflow",
                      freq.hz != "NA",), 
             pch = 21, color = "black",
             aes(y=thresh.ms2, x = svl.mm, fill = combined.id), 
             size = 3, alpha=1) +
  
  geom_line(data = vib.thresh %>%
               filter(life.stage == "juvenile",
                      treatment != "overflow",
                      freq.hz != "NA",), 
             pch = 21, 
             aes(y=thresh.ms2, x = svl.mm, color = combined.id), 
             size = 0.5, alpha=1) +
  
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=14, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=14, color = "black"),
        axis.title.y = element_text(size=28),
        axis.title.x = element_text(size=28),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "svl (mm)", limits = c(15,38), breaks = seq(15,38,5), labels = seq(15,38,5)) +
  scale_y_continuous(name = "vib. threshold (dB re 1 m/s-2)")


hear.thresh.plot.svl.indiv <- ggplot() +
  
  facet_grid(cols = vars(factor(freq.hz))) +
  
  geom_point(data = hear.thresh %>%
               filter(life.stage == "juvenile",
                      treatment != "overflow",
                      freq.hz != "NA",), 
             pch = 21, color = "black",
             aes(y=thresh.db, x = svl.mm, fill = combined.id), ## SMS YOU ARE HERE, color by freq,
             size = 3, alpha=1) +
  
  geom_line(data = hear.thresh %>%
              filter(life.stage == "juvenile",
                     treatment != "overflow",
                     freq.hz != "NA",), 
            pch = 21, aes(y=thresh.db, x = svl.mm, color = combined.id), ## SMS YOU ARE HERE, color by freq,
            size = 0.5, alpha=1) +
  
  
  geom_point(data = hear.thresh[hear.thresh$life.stage == "juvenile" & hear.thresh$treatment != "overflow",], pch = 21,
             aes(y=clip.thresh.db.plus3, x = svl.mm, fill = "white", color = combined.id), ## SMS YOU ARE HERE, color by freq,
             size = 3, alpha=1) +
  
  geom_line(data = hear.thresh[hear.thresh$life.stage == "juvenile" & hear.thresh$treatment != "overflow",], pch = 21, 
             aes(y=clip.thresh.db.plus3, x = svl.mm, color = combined.id), ## SMS YOU ARE HERE, color by freq,
             size = 0.5, alpha=1) +
  
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=14, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=14, color = "black"),
        axis.title.y = element_text(size=28),
        axis.title.x = element_text(size=28),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "svl (mm)", limits = c(15,38), breaks = seq(15,38,5), labels = seq(15,38,5)) +
  scale_y_continuous(name = "hearing threshold (dB)", limits = c(80,140), breaks = seq(80,140,15), labels = seq(80,140,15))


png("~/Desktop/R Working Directory/Plots/Juv-SVL-Indiv.png", units = "in", res = 300, width = 16, height = 12)
ggarrange(vib.thresh.plot.svl.indiv, hear.thresh.plot.svl.indiv, 
          ncol = 1,
          nrow = 2,
          common.legend = FALSE,
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))
dev.off()


# PLOT: plotting by individual -------------------------------------------------------------
complete = c("RM_J001_10", "RM_J004_9", "RM_J005_7", "RM_J007_9", "RM_J010_6", "RM_J011_12", "RM_J020_7", "RM_J023_8")

vib.thresh.plot.indiv <- ggplot() +
  
  facet_grid(cols = vars(factor(combined.id))) +
  
  geom_point(data = vib.thresh %>%
               filter(life.stage == "juvenile",
                      treatment != "overflow",
                      freq.hz != "NA",
                      combined.id %in% complete), 
             pch = 21, color = "black",
             aes(y=thresh.ms2, x = freq.hz, fill = factor(num.sampling)), 
             size = 3, alpha=1) +
  
  geom_line(data = vib.thresh %>%
              filter(life.stage == "juvenile",
                     treatment != "overflow",
                     freq.hz != "NA",
                     combined.id %in% complete), 
            pch = 21, 
            aes(y=thresh.ms2, x = freq.hz, color = factor(num.sampling)), 
            size = 0.5, alpha=1) +
  
  scale_colour_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  scale_fill_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  
  theme_bw() +
  theme(strip.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=10, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=10, color = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "frequency (kHz)", limits = c(100,1500), breaks = sort(unique(vib.thresh$freq.hz)), labels = c("0.1","0.2", "0.3", "0.4", "0.7", "1.1", "1.5")) +
  scale_y_continuous(name = "vib. threshold (dB re 1 m/s-2)")


hear.thresh.plot.indiv <- ggplot() +
  
  facet_grid(cols = vars(factor(combined.id))) +
  
  geom_point(data = hear.thresh %>%
               filter(life.stage == "juvenile",
                      treatment != "overflow",
                      freq.hz != "NA",
                      combined.id %in% complete), 
             pch = 21, color = "black",
             aes(y=clip.thresh.db.plus3, x = freq.hz, fill = factor(num.sampling)), ## SMS YOU ARE HERE, color by freq,
             size = 3, alpha=1) +
  
  #threshold above clipping level
  geom_point(data = hear.thresh %>%
               filter(life.stage == "juvenile",
                      treatment != "overflow",
                      is.na(thresh.db) == TRUE,
                      freq.hz != "NA",
                      combined.id %in% complete),
             pch = 21, fill = "white",
             aes(y=clip.thresh.db.plus3, x = freq.hz, color = factor(num.sampling)), ## SMS YOU ARE HERE, color by freq,
             size = 3, alpha=1) +
  
  geom_line(data = hear.thresh %>%
              filter(life.stage == "juvenile",
                     treatment != "overflow",
                     freq.hz != "NA",
                     combined.id %in% complete),
            aes(y=clip.thresh.db.plus3, x = freq.hz, color = factor(num.sampling)), ## SMS YOU ARE HERE, color by freq,
            size = 0.5, alpha=1) +
  
  # extrapolated points because threshold above clipping level
  geom_point(data = hear.thresh %>%
               filter(life.stage == "juvenile",
                      treatment != "overflow",
                      freq.hz != "NA",
                      is.na(thresh.db) == TRUE,
                      combined.id %in% complete), 
             pch = 21, fill = "white", 
             aes(y=clip.thresh.db.plus3, x = freq.hz, color = factor(num.sampling)),
             size = 3, alpha=1) +
  
  scale_colour_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  scale_fill_manual(values = c("#CC66FF", "#660066", "#330066"), labels = c("3 months", "6 months", "12 months")) +
  
  theme_bw() +
  theme(strip.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=10, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=10, color = "black"),
        axis.title.y = element_text(size=24),
        axis.title.x = element_text(size=24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "frequency (kHz)", limits = c(100,1500), breaks = sort(unique(vib.thresh$freq.hz)), labels = c("0.1","0.2", "0.3", "0.4", "0.7", "1.1", "1.5")) +
  scale_y_continuous(name = "hearing threshold (dB)", limits = c(80,140), breaks = seq(80,140,15), labels = seq(80,140,15))


png("~/Desktop/R Working Directory/Plots/Juv-Indiv.png", units = "in", res = 300, width = 16, height = 12)
ggarrange(vib.thresh.plot.indiv, hear.thresh.plot.indiv, 
          ncol = 1,
          nrow = 2,
          common.legend = TRUE,
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))
dev.off()






# PLOT: Body morphology by mass ----

plot.morph.corr.svl <- ggplot() +
  
  geom_point(data = morph.data.juv[morph.data.juv$data.type == "ephys" & morph.data.juv$gs.code == "RM" & morph.data.juv$treatment != "overflow",],
             aes(x = mass.g, y = svl.mm, color = num.sampling.six.cat),
             size = 3) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=18, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=18, color = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "mass (g)", limits = c(0.3,2)) +
  scale_y_continuous(name = "snout-vent length (mm)")

plot.morph.corr.interear <- ggplot() +
  
  geom_point(data = morph.data.juv[morph.data.juv$data.type == "ephys" & morph.data.juv$gs.code == "RM" & morph.data.juv$treatment != "overflow",],
             aes(x = mass.g, y = inter.ear.mm, color = num.sampling.six.cat),
             size = 3) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=18, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=18, color = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "mass (g)", limits = c(0.3,2)) +
  scale_y_continuous(name = "inter-ear width (mm)")


ggarrange(plot.morph.corr.svl, plot.morph.corr.interear,
          nrow = 2,
          common.legend =  TRUE)


# PLOT: Body morphology by larval duration ----
plot.larvdur.corr.mass <- ggplot() +
  
  geom_point(data = morph.data.juv[morph.data.juv$data.type == "ephys" & morph.data.juv$gs.code == "RM" & morph.data.juv$treatment != "overflow",],
             aes(x = mean.days.forelimb, y = mass.g, color = num.sampling.six.cat),
             size = 3) +
  
  geom_smooth(method = "lm",
              data = morph.data.juv[morph.data.juv$data.type == "ephys" & morph.data.juv$gs.code == "RM" & morph.data.juv$treatment != "overflow",],
              aes(y=mass.g, x = mean.days.forelimb, color = factor(num.sampling.six.cat), group = factor(num.sampling.six.cat))) +
  
  
  geom_smooth(method = "lm",
              data = morph.data.juv[morph.data.juv$data.type == "ephys" & morph.data.juv$gs.code == "RM" & morph.data.juv$treatment != "overflow",],
              color = "black",
              aes(y=mass.g, x = mean.days.forelimb)) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=18, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=18, color = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "mean larval duration (days)", limits = c(38,70)) +
  scale_y_continuous(name = "mass (g)")


plot.larvdur.corr.svl <- ggplot() +
  
  geom_point(data = morph.data.juv[morph.data.juv$data.type == "ephys" & morph.data.juv$gs.code == "RM" & morph.data.juv$treatment != "overflow",],
             aes(x = mean.days.forelimb, y = svl.mm, color = num.sampling.six.cat),
             size = 3) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=18, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=18, color = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "mean larval duration (days)", limits = c(38,70)) +
  scale_y_continuous(name = "snout-vent length (mm)")

plot.larvdur.corr.interear <- ggplot() +
  
  geom_point(data = morph.data.juv[morph.data.juv$data.type == "ephys" & morph.data.juv$gs.code == "RM" & morph.data.juv$treatment != "overflow",],
             aes(x = mean.days.forelimb, y = inter.ear.mm, color = num.sampling.six.cat),
             size = 3) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=18, color = "black", angle = 0, hjust = 0.5), 
        axis.text.y=element_text(size=18, color = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "mean larval duration (days)", limits = c(38,70)) +
  scale_y_continuous(name = "inter-ear width (mm)")


ggarrange(plot.larvdur.corr.mass, plot.larvdur.corr.svl, plot.larvdur.corr.interear,
          ncol = 3,
          common.legend =  TRUE)








# Redoing plots after ranid devo paper -------------------------

#plot cm/s2 individual lines colored by treatment and showing means

ggplot(vib.thresh[vib.thresh$Genus.species.code == "RS",]) + 
  
  facet_grid(cols = vars(num.sampling), drop = FALSE) +
  
  #individual tanks
  geom_point(size = 1, alpha = 0.7, pch = 21, aes(y=thresh.cms2, x = Freq..Hz., color = Treatment, fill = Treatment)) +
  geom_line(size = 0.5, alpha = 0.7, aes(y=thresh.cms2, x = Freq..Hz., color = Treatment, group = combined.id)) +
  
  #treatment means
  stat_summary(fun=mean, geom="line", size = 0.8, color = "black", aes(x = Freq..Hz., y = thresh.cms2, group = Treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=20, size = 0.8, colour="black", alpha=1, aes(x = Freq..Hz., y = thresh.cms2, group = Treatment)) +
  stat_summary(fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = Freq..Hz., y = thresh.cms2, fill=Treatment), show.legend = TRUE) + 
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "threshold (cm/s2)") +
  scale_x_continuous(name = "frequency (Hz)", breaks = unique(vib.thresh$Freq..Hz.))


#plot cm/s2 individual lines colored by individual and showing means
ggplot(vib.thresh[vib.thresh$Genus.species.code == "RS",]) + 
  
  facet_grid(cols = vars(num.sampling), drop = FALSE) +
  
  #individual tanks
  geom_point(size = 1, alpha = 0.7, pch = 21, aes(y=thresh.mms2, x = Freq..Hz., color = unique.id.juv, fill = unique.id.juv)) +
  geom_line(size = 0.8, alpha = 0.7, aes(y=thresh.mms2, x = Freq..Hz., color = unique.id.juv, group = unique.id.juv)) +
  
  #treatment means
  stat_summary(fun=mean, geom="line", size = 0.8, color = "black", aes(x = Freq..Hz., y = thresh.mms2), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=20, size = 0.8, colour="black", alpha=1, aes(x = Freq..Hz., y = thresh.mms2), show.legend = FALSE) +
  stat_summary(fun = mean, geom="point", color = "black", fill = "black", pch=21, size=5, stroke = 1, aes(x = Freq..Hz., y = thresh.mms2), show.legend = FALSE) + 
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "bottom",
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "threshold (mm/s2)") +
  scale_x_continuous(name = "frequency (Hz)", breaks = unique(vib.thresh$Freq..Hz.))


#plot faceted by individuals with each timepoint as different color
ggplot(data = vib.thresh) + 
  
  facet_grid(cols = vars(unique.id.juv)) +
  
  #individuals
  geom_point(size = 4, stroke = 1.5, alpha = 0.7, pch = 21, aes(y=Threshold..db., x = as.numeric(Freq..Hz.), color = num.sampling, fill = num.sampling, group = num.sampling),
             position=position_jitterdodge(jitter.width = 20, jitter.height = 0, seed = 0)) +
  
  geom_line(size = 0.5, alpha = 0.7, aes(y=Threshold..db., x = as.numeric(Freq..Hz.), color = num.sampling, group = num.sampling),
            position=position_jitterdodge(jitter.width = 5, jitter.height = 0, seed = 0)) +
  
  scale_colour_gradient(low="lightgray", high = "black") +
  scale_fill_gradient(low="lightgray", high = "black") +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "bottom",
        axis.text.x=element_text(size=20, color = "black", angle = 90, vjust = 0.5, hjust=1), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "threshold (cm/s2)") +
  scale_x_continuous(name = "frequency (kHz)", breaks = sort(unique(vib.thresh$Freq..Hz.)), labels = c("0.1", "0.2"," 0.3", "0.4", "0.5", "0.7", "1.1", "1.5"))



#plot faceted by individuals (LOW DENSITY)
indiv.ld <- ggplot(data = vib.thresh[vib.thresh$combined.id != "RSJ00321" & vib.thresh$Treatment != "low density",]) + 
  
  facet_grid(cols = vars(combined.id)) +
  
  #individuals
  geom_point(size = 4, stroke = 1.5, alpha = 0.7, pch = 21, color = c(natparks.pals("Arches")[4]), aes(y=Threshold..db., x = as.numeric(Freq..Hz.), fill = as.factor(num.sampling), group = as.factor(num.sampling)),
             position=position_jitterdodge(jitter.width = 20, jitter.height = 0, seed = 0)) +
  scale_fill_manual(values=c("darkgray", "black")) +
  
  geom_line(size = 0.5, alpha = 0.7, aes(y=Threshold..db., x = as.numeric(Freq..Hz.), color = as.factor(num.sampling), group = as.factor(num.sampling)),
            position=position_jitterdodge(jitter.width = 5, jitter.height = 0, seed = 0)) +
  
  scale_color_manual(values=c("darkgray", "black")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.text.x=element_text(size=20, color = "black", angle = 90, vjust = 0.5, hjust=1), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "threshold", limits = c(45,92)) +
  scale_x_continuous(name = "frequency (kHz)", breaks = sort(unique(vib.thresh$Freq..Hz.)), labels = c("0.1", "0.2"," 0.3", "0.4"," 0.7", "1.1", "1.5"))


#plot faceted by individuals (HIGH DENSITY)
indiv.hd <- ggplot(data = vib.thresh[vib.thresh$combined.id != "RSJ00321" & vib.thresh$Treatment != "high density",]) + 
  
  facet_grid(cols = vars(combined.id)) +
  
  #individuals
  geom_point(size = 4, stroke = 1.5, alpha = 0.7, pch = 21, color = c(natparks.pals("Arches")[1]), aes(y=Threshold..db., x = as.numeric(Freq..Hz.), fill = as.factor(num.sampling), group = as.factor(num.sampling)),
             position=position_jitterdodge(jitter.width = 20, jitter.height = 0, seed = 0)) +
  scale_fill_manual(values=c("darkgray", "black")) +
  
  geom_line(size = 0.5, alpha = 0.7, aes(y=Threshold..db., x = as.numeric(Freq..Hz.), color = as.factor(num.sampling), group = as.factor(num.sampling)),
            position=position_jitterdodge(jitter.width = 5, jitter.height = 0, seed = 0)) +
  
  scale_color_manual(values=c("darkgray", "black")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.text.x=element_text(size=20, color = "black", angle = 90, vjust = 0.5, hjust=1), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "threshold", limits = c(45,92)) +
  scale_x_continuous(name = "frequency (kHz)", breaks = sort(unique(vib.thresh$Freq..Hz.)), labels = c("0.1", "0.2"," 0.3", "0.4"," 0.7", "1.1", "1.5"))

#panel plot
ggarrange(indiv.ld, indiv.hd,
          nrow = 2)


#plot faceted by individuals (LOW DENSITY) + MASS
deltamass.ld <- ggplot(data = vib.thresh[vib.thresh$combined.id != "RSJ00321" & vib.thresh$Treatment != "low density",]) + 
  
  facet_grid(cols = vars(combined.id)) +
  
  #individuals
  geom_line(size = 1, alpha = 1, color = c(natparks.pals("Arches")[4]), aes(y=Weight..g., x = as.numeric(num.sampling))) +
  geom_point(size = 4, stroke = 1.5, alpha = 1, pch = 21, color = c(natparks.pals("Arches")[4]), aes(y=Weight..g., x = as.numeric(num.sampling), fill = as.factor(num.sampling))) +

  scale_color_manual(values=c("darkgray", "black")) +
  scale_fill_manual(values=c("darkgray", "black")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "mass (g)", limits = c(0.3, 1.8), breaks = seq(0.4,1.8,0.2)) +
  scale_x_continuous(name = "post-metamorphic age (weeks)", limits = c(1,2), breaks = c(1,2), labels = c("11-12", "31-31"))


#plot faceted by individuals (HIGH DENSITY) + MASS
deltamass.hd <- ggplot(data = vib.thresh[vib.thresh$combined.id != "RSJ00321" & vib.thresh$Treatment != "high density",]) + 
  
  facet_grid(cols = vars(combined.id)) +
  
  #individuals
  geom_line(size = 1, alpha = 1, color = c(natparks.pals("Arches")[1]), aes(y=Weight..g., x = as.numeric(num.sampling))) +
  geom_point(size = 4, stroke = 1.5, alpha = 1, pch = 21, color = c(natparks.pals("Arches")[1]), aes(y=Weight..g., x = as.numeric(num.sampling), fill = as.factor(num.sampling))) +
  
  scale_color_manual(values=c("darkgray", "black")) +
  scale_fill_manual(values=c("darkgray", "black")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "mass (g)", limits = c(0.3, 1.8), breaks = seq(0.4,1.8,0.2)) +
  scale_x_continuous(name = "post-metamorphic age (weeks)", limits = c(1,2), breaks = c(1,2), labels = c("11-12", "31-31"))

#panel plot
ggarrange(deltamass.ld, deltamass.hd,
          nrow = 2)


#plot threshold by mass and treatment
ggplot() + 
  
  facet_grid(rows = vars(num.sampling), cols = vars(Freq..Hz.)) +
  
  geom_point(position=position_jitterdodge(jitter.width = 0, jitter.height = 0, seed = 0), pch = 21, size = 2, alpha = 0.7, show.legend = FALSE,
             data = vib.thresh,
             aes(y=Threshold..db., x = Weight..g., color = Treatment, fill = Treatment)) +
  #geom_ribbon(data = predicted.df.morph.mm,
              #alpha = 0.4,
              #mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), show.legend = FALSE) +
  #geom_line(data = predicted.df.morph.mm, size = 1, color = "black", aes(x = x, y = predicted, group = group)) +
  
  geom_smooth(data = vib.thresh, method = "lm",
              mapping = aes(y=Threshold..db., x = Weight..g., color = Treatment, fill = Treatment),
              se=FALSE) +
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "threshold") +
  scale_x_continuous(name = "mass (g)")


#plot threshold by mass
ggplot() + 
  
  facet_grid(rows = vars(num.sampling), cols = vars(Freq..Hz.)) +
  
  geom_point(pch = 21, size = 2, alpha = 0.7, show.legend = FALSE,
             data = vib.thresh,
             aes(y=Threshold..db., x = Weight..g.)) +
  #geom_ribbon(data = predicted.df.morph.mm,
  #alpha = 0.4,
  #mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), show.legend = FALSE) +
  #geom_line(data = predicted.df.morph.mm, size = 1, color = "black", aes(x = x, y = predicted, group = group)) +
  
  geom_smooth(data = vib.thresh, method = "lm",
              mapping = aes(y=Threshold..db., x = Weight..g.),
              se=FALSE) +
  
  #scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  #scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "threshold") +
  scale_x_continuous(name = "mass (g)")


#plots individual changes
unique(vib.thresh$combined.id[vib.thresh$num.sampling == 2])

vib.thresh.change <- data.frame(
  combined.id = sort(rep(unique(vib.thresh$combined.id[vib.thresh$num.sampling == 2]), length(unique(vib.thresh$Freq..Hz.)))),
  freq.hz = rep(sort(unique(vib.thresh$Freq..Hz.)), length(unique(vib.thresh$combined.id[vib.thresh$num.sampling == 2]))),
  treatment = factor(NA, levels = c("high density", "low density", "overflow")) 
)

for(i in 1:length(unique(vib.thresh.change$combined.id))){
  
  for(j in 1:length(unique(vib.thresh.change$freq.hz))){
    
    vib.thresh.change$treatment[vib.thresh.change$combined.id == unique(vib.thresh.change$combined.id)[i]] = factor(vib.thresh$Treatment[vib.thresh$combined.id == unique(vib.thresh.change$combined.id)[i]][1], levels = c("high density", "low density", "overflow")) 
    
    vib.thresh.change$thresh.1[vib.thresh.change$combined.id == unique(vib.thresh.change$combined.id)[i] & vib.thresh.change$freq.hz == unique(vib.thresh.change$freq.hz)[j]] = vib.thresh$Threshold..db.[vib.thresh$combined.id == unique(vib.thresh.change$combined.id)[i] & vib.thresh$Freq..Hz. == unique(vib.thresh.change$freq.hz)[j] & vib.thresh$num.sampling == 1] 
   
    
    vib.thresh.change$thresh.2[vib.thresh.change$combined.id == unique(vib.thresh.change$combined.id)[i] & vib.thresh.change$freq.hz == unique(vib.thresh.change$freq.hz)[j]] = vib.thresh$Threshold..db.[vib.thresh$combined.id == unique(vib.thresh.change$combined.id)[i] & vib.thresh$Freq..Hz. == unique(vib.thresh.change$freq.hz)[j] & vib.thresh$num.sampling == 2] 
    
    }
  }

vib.thresh.change$delta.thresh = vib.thresh.change$thresh.2 - vib.thresh.change$thresh.1

vib.thresh.change$treatment = factor(vib.thresh.change$treatment, levels = c("low density", "high density", "overflow"))

ggplot() + 
  
  geom_point(position=position_jitterdodge(jitter.width = 5, jitter.height = 0.5, seed = 0),
             size = 2, alpha = 0.7, show.legend = TRUE,
             data = vib.thresh.change,
             aes(y=delta.thresh, x = freq.hz, color = treatment, group = combined.id)) +

  geom_line(position=position_jitterdodge(jitter.width = 5, jitter.height = 0.5, seed = 0),
             size = 0.5, alpha = 0.7, show.legend = TRUE,
             data = vib.thresh.change,
             aes(y=delta.thresh, x = freq.hz, color = treatment, group = combined.id)) +


  #add horizontal zero line
  geom_hline(yintercept = 0, colour = "black", lty = 3) +

  #treatment means
  stat_summary(data = vib.thresh.change, fun=mean, geom="line", size = 0.8, color = "black", aes(x = freq.hz, y = delta.thresh, group = treatment), show.legend = TRUE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(data = vib.thresh.change, fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=20, size = 0.8, colour="black", alpha=1, aes(x = freq.hz, y = delta.thresh, group = treatment)) +
  stat_summary(data = vib.thresh.change, fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = freq.hz, y = delta.thresh, fill=treatment), show.legend = TRUE) + 
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "delta threshold", breaks = seq(-18,10,3)) +
  scale_x_continuous(name = "frequency (Hz)", breaks = unique(vib.thresh.change$freq.hz))




#  PLOT RELATIVE THRESHOLD DB plot individual lines colored by treatment and showing means
ggplot(vib.thresh[vib.thresh$Genus.species.code == "RS",]) + 
  
  facet_grid(cols = vars(postmm.sampling), drop = FALSE) +
  
  #individual tanks
  geom_point(size = 1, alpha = 0.7, pch = 21, aes(y=rel.dB, x = Freq..Hz., color = Treatment, fill = Treatment)) +
  geom_line(size = 0.5, alpha = 0.7, aes(y=rel.dB, x = Freq..Hz., color = Treatment, group = combined.id)) +
  
  #treatment means
  stat_summary(fun=mean, geom="line", size = 0.8, color = "black", aes(x = Freq..Hz., y = rel.dB, group = Treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=20, size = 0.8, colour="black", alpha=1, aes(x = Freq..Hz., y = rel.dB, group = Treatment)) +
  stat_summary(fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = Freq..Hz., y = rel.dB, fill=Treatment), show.legend = TRUE) + 
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "relative threshold (dB)") +
  scale_x_continuous(name = "frequency (Hz)", breaks = unique(vib.thresh$Freq..Hz.))







