# Script To-Do List --------------------------------------------------------------
#1. Update Fig. 3 to be a function to create all iterations of figure in ggarrange call
#2. Finalize figures

# Load required packages and set working directory --------------------------------------------------------------

rm(list=ls())
required_pckg = c("dplyr", "readr", "NatParksPalettes", "ggpubr", "RColorBrewer",
                  "lme4", "MuMIn", "DHARMa", "car", "stats", "latticeExtra", "cAIC4", "stargazer", "sjPlot")
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
library("stargazer")


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
         combined.id.freq = factor(paste(combined.id, freq.hz, sep = "_")),
         # convert vibration threshold to m/s2 [[acceleration (m/s2) = (10^((rawdB-120)/20))*9.81]]
         thresh.ms2 = (10^((thresh.db-120)/20))*9.81,
         #thresh.cms2 = 100*((10^((thresh.db-120)/20))*9.81),
         #thresh.mms2 = 1000*((10^((thresh.db-120)/20))*9.81),
         clip.thresh.db.plus3 = if_else(is.na(clip.thresh.db.plus3) == TRUE & is.na(thresh.db) == FALSE, 
                                        thresh.db,
                                        clip.thresh.db.plus3), #allows better plotting because individual can have connected lines across frequencies, but then code overlays white points where clipping occurred
         rel.db = (thresh.db-120),
         clip.thresh.rel.db.plus3 = (clip.thresh.db.plus3-120),
         clip.thresh.db.plus3.ms2 = (10^((clip.thresh.db.plus3-120)/20))*9.81,
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
    !(date == "2025-05-06" & combined.id == "RM_J041_8"),
    !(date == "2025-05-12" & combined.id == "RM_J041_8"),
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
         combined.id.freq = factor(paste(combined.id, freq.hz, sep = "_")),
         clip.thresh.db.plus3 = if_else(is.na(clip.thresh.db.plus3) == TRUE & is.na(thresh.db) == FALSE, 
                                        thresh.db,
                                        clip.thresh.db.plus3), #allows better plotting because individual can have connected lines across frequencies, but then code overlays white points where clipping occurred
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
    !(date == "2025-05-06" & combined.id == "RM_J041_8"), #remove problematic individual J041_08
    !(date == "2025-05-12" & combined.id == "RM_J041_8"), #remove problematic individual J041_08
    !(date == "2025-05-05" & combined.id == "RM_J041_9"), #remove problematic date for one individual J041_09
        combined.id != "RM_NA_1", #odd adult
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

# create clean datasets without NAs to use, if needed
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

vib.thresh.summary.reldB <- vib.thresh %>%
  group_by(life.stage, num.sampling, freq.hz, .drop = TRUE) %>%
  summarise(n = n(), 
            mean.thresh.reldB = mean(rel.dB, na.rm= TRUE), 
            sd.thresh.reldB = sd(rel.dB, na.rm= TRUE),
            mean.mass.g = mean(mass.g, na.rm = TRUE),
            mean.svl.mm = mean(svl.mm, na.rm = TRUE))


hear.thresh.summary <- hear.thresh %>%
  group_by(life.stage, num.sampling, freq.hz, .drop = TRUE) %>%
  summarise(n = n(), 
            mean.thresh.db = mean(thresh.db, na.rm= TRUE), 
            sd.thresh.db = sd(thresh.db, na.rm= TRUE),
            mean.mass.g = mean(mass.g, na.rm = TRUE),
            mean.svl.mm = mean(svl.mm, na.rm = TRUE))

# Create difference-based databases --------------------------

#function to create database that sets threshold to 0 at initial testing then calculates threshold difference between timepoints
thresh_diff_db <- function(db){
  db.temp <- bind_rows(
    # join first sampling point with next two sampling points, filling in the thresh.diff column as the sampling points are joined
    # using thresh.diff.clip when threshold above clipping
    db %>%
      mutate(thresh.db.diff = case_when(is.na(thresh.db) == FALSE ~ 0,
                                        is.na(thresh.db) == TRUE ~ NA),
             clip.thresh.db.diff = case_when(is.na(thresh.db) == FALSE ~ NA,
                                             is.na(thresh.db) == TRUE ~ 0)
      ) %>%
      filter(num.sampling == 1),
    
    db %>% #calculate difference in threshold between first and second sampling, based on joining to id and frequency
      filter(num.sampling == 2) %>%
      left_join( 
        db %>%
          filter(num.sampling == 1) %>% 
          mutate(thresh.db.1 = thresh.db,
                 clip.thresh.db.1 = clip.thresh.db.plus3) %>%
          select(combined.id.freq, thresh.db.1, clip.thresh.db.1) %>%
          left_join(
            db %>%
              filter(num.sampling == 2) %>%
              mutate(thresh.db.2 = thresh.db,
                     clip.thresh.db.2 = clip.thresh.db.plus3) %>%
              select(combined.id.freq, thresh.db.2, clip.thresh.db.2),
            by = "combined.id.freq"
          )  %>%
          mutate(thresh.db.diff = thresh.db.2-thresh.db.1,
                 clip.thresh.db.diff = clip.thresh.db.2 - clip.thresh.db.1) %>%
          select(combined.id.freq, thresh.db.diff, clip.thresh.db.diff),
        by = "combined.id.freq"
      ),
    
    db %>% #calculate difference in threshold between first and third sampling, based on joining to id and frequency
      filter(num.sampling == 3) %>%
      left_join( 
        db %>%
          filter(num.sampling == 1) %>% 
          mutate(thresh.db.1 = thresh.db,
                 clip.thresh.db.1 = clip.thresh.db.plus3) %>%
          select(combined.id.freq, thresh.db.1, clip.thresh.db.1) %>%
          left_join(
            db %>%
              filter(num.sampling == 3) %>%
              mutate(thresh.db.2 = thresh.db,
                     clip.thresh.db.2 = clip.thresh.db.plus3) %>%
              select(combined.id.freq, thresh.db.2, clip.thresh.db.2),
            by = "combined.id.freq"
          )  %>%
          mutate(thresh.db.diff = thresh.db.2-thresh.db.1,
                 clip.thresh.db.diff = clip.thresh.db.2 - clip.thresh.db.1) %>%
          select(combined.id.freq, thresh.db.diff, clip.thresh.db.diff),
        by = "combined.id.freq"
      )
  )
  assign(paste(deparse(substitute(db)), "diff", sep = "."), db.temp, envir = .GlobalEnv)
}

# run function for each database
thresh_diff_db(vib.thresh.juv)
thresh_diff_db(hear.thresh.juv)

# Can double-check difference calculations with the following function check_diff
check_diff = function(db, indiv) {
  View(db %>% filter(combined.id == indiv) %>% select(date, num.sampling, freq.hz, thresh.db, thresh.db.diff, clip.thresh.db.diff))
}
check_diff(vib.thresh.juv.diff, "RM_J001_10")
check_diff(hear.thresh.juv.diff, "RM_J001_10")

# Create function to determine random effects structure ------------
thresh_reffects_juv = function(x){

  #random effect for intercept
lmm.full.1 <- lmer(thresh.db ~ svl.mm*mean.days.forelimb*poly(freq.hz, 3, raw= FALSE)*factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)

#random effect for intercept and slope
lmm.full.1.slopes.1 <- lmer(thresh.db ~ svl.mm*mean.days.forelimb*poly(freq.hz, 3, raw= FALSE)*factor(num.sampling) + (1+factor(num.sampling)|combined.id), data = x, na.action = na.omit)

#fixed intercept, random effect for slope
lmm.full.1.slopes.2 <- lmer(thresh.db ~ svl.mm*mean.days.forelimb*poly(freq.hz, 3, raw= FALSE)*factor(num.sampling) + (0+factor(num.sampling)|combined.id), data = x, na.action = na.omit)

if((anova(lmm.full.1, lmm.full.1.slopes.1, lmm.full.1.slopes.2)$"Pr(>Chisq)"[2] < 0.05 &
    is.na(anova(lmm.full.1, lmm.full.1.slopes.1, lmm.full.1.slopes.2)$"Pr(>Chisq)"[3]) == TRUE)){
      print("use random effect for intercept and slope")}else{
        if(anova(lmm.full.1, lmm.full.1.slopes.1, lmm.full.1.slopes.2)$"Pr(>Chisq)"[3] < 0.05 &
        is.na(anova(lmm.full.1, lmm.full.1.slopes.1, lmm.full.1.slopes.2)$"Pr(>Chisq)"[3]) == FALSE){
        print("use fixed intercept and random effect for slope")}else{
          print("use random effect for intercept model")}
      }
}

thresh_reffects = function(x){
  
  #random effect for intercept
  lmm.full.1 <- lmer(thresh.db ~ svl.mm*poly(freq.hz, 3, raw= FALSE)*factor(life.stage.num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  #random effect for intercept and slope
  lmm.full.1.slopes.1 <- lmer(thresh.db ~ svl.mm*poly(freq.hz, 3, raw= FALSE)*factor(life.stage.num.sampling) + (1+factor(num.sampling)|combined.id), data = x, na.action = na.omit)
  
  #fixed intercept, random effect for slope
  lmm.full.1.slopes.2 <- lmer(thresh.db ~ svl.mm*poly(freq.hz, 3, raw= FALSE)*factor(life.stage.num.sampling) + (0+factor(num.sampling)|combined.id), data = x, na.action = na.omit)
  
  if((anova(lmm.full.1, lmm.full.1.slopes.1, lmm.full.1.slopes.2)$"Pr(>Chisq)"[2] < 0.05 &
      is.na(anova(lmm.full.1, lmm.full.1.slopes.1, lmm.full.1.slopes.2)$"Pr(>Chisq)"[3]) == TRUE)){
    print("use random effect for intercept and slope")}else{
      if(anova(lmm.full.1, lmm.full.1.slopes.1, lmm.full.1.slopes.2)$"Pr(>Chisq)"[3] < 0.05 &
         is.na(anova(lmm.full.1, lmm.full.1.slopes.1, lmm.full.1.slopes.2)$"Pr(>Chisq)"[3]) == FALSE){
        print("use fixed intercept and random effect for slope")}else{
          print("use random effect for intercept model")}
    }
}

# Create function for juvenile-only model comparison -----------------------------------------
thresh_model_compare_juv <- function(x){ 

  #candidate models - manually defined
  #does the effect of age (num.sampling) depend on development time, size, and frequency?
  lmm.full.1 <- lmer(thresh.db ~ svl.mm*mean.days.forelimb*poly(freq.hz, 3, raw= FALSE)*factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  lmm.full.1.ordinal <- lmer(thresh.db ~ svl.mm*mean.days.forelimb*ordered(freq.hz, levels = unique(x$freq.hz))*factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age depend on the effects of development time and size similarly across frequency?
  lmm.3intx.1 <- lmer(thresh.db ~ svl.mm*mean.days.forelimb*factor(num.sampling) + poly(freq.hz, 3, raw= FALSE) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age function independently from the effects of size, frequency, and age?
  lmm.3intx.2 <- lmer(thresh.db ~ mean.days.forelimb*svl.mm*poly(freq.hz, 3, raw= FALSE) + factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age depend on development time and frequency?
  lmm.3intx.3 <- lmer(thresh.db ~ mean.days.forelimb*poly(freq.hz, 3, raw= FALSE)*factor(num.sampling) + svl.mm + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age depend on development time and frequency?
  lmm.3intx.4 <- lmer(thresh.db ~ mean.days.forelimb + svl.mm*poly(freq.hz, 3, raw= FALSE)*factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age depend on frequency and the effects of size depend on development time?
  lmm.2intx.1 <- lmer(thresh.db ~ svl.mm*mean.days.forelimb + poly(freq.hz, 3, raw= FALSE)*factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of size depend on frequency and the effects of development time depend on age?
  lmm.2intx.2 <- lmer(thresh.db ~ svl.mm*poly(freq.hz, 3, raw= FALSE) + mean.days.forelimb*factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of size depend on age and the effects of development time depend on frequency?
  lmm.2intx.3 <- lmer(thresh.db ~ svl.mm*factor(num.sampling) + mean.days.forelimb*poly(freq.hz, 3, raw= FALSE) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of size depend on development time?
  lmm.1intx.1 <- lmer(thresh.db ~ svl.mm*mean.days.forelimb + poly(freq.hz, 3, raw= FALSE) + factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age depend on development time?
  lmm.1intx.2 <- lmer(thresh.db ~ factor(num.sampling)*mean.days.forelimb + poly(freq.hz, 3, raw= FALSE) + svl.mm + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of frequency depend on development time?
  lmm.1intx.3 <- lmer(thresh.db ~ poly(freq.hz, 3, raw= FALSE)*mean.days.forelimb + factor(num.sampling) + svl.mm + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of size depend on age?
  lmm.1intx.4 <- lmer(thresh.db ~ svl.mm*factor(num.sampling) + poly(freq.hz, 3, raw= FALSE) + mean.days.forelimb + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of size depend on frequency?
  lmm.1intx.5 <- lmer(thresh.db ~ svl.mm*poly(freq.hz, 3, raw= FALSE) + factor(num.sampling) + mean.days.forelimb + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of age depend on frequency?
  lmm.1intx.6 <- lmer(thresh.db ~ svl.mm + factor(num.sampling)*poly(freq.hz, 3, raw= FALSE) + mean.days.forelimb + (1|combined.id), data = x, na.action = na.omit)

  #additive effects of age, size, development time, and frequency?
  lmm.add <- lmer(thresh.db ~ svl.mm + mean.days.forelimb + poly(freq.hz, 3, raw= FALSE) + factor(num.sampling) + (1|combined.id), data = x, na.action = na.omit, REML = TRUE)
  
  #null model
  lmm.null <- lmer(thresh.db ~ (1|combined.id), data = x, na.action = na.omit)
  
  #model comparison using AICc
  model.sel = arrange(AICc(
                  lmm.full.1, lmm.full.1.ordinal,
                  lmm.3intx.1, lmm.3intx.2, lmm.3intx.3, lmm.3intx.4,
                  lmm.2intx.1, lmm.2intx.2, lmm.2intx.3,
                  lmm.1intx.1, lmm.1intx.2, lmm.1intx.3, lmm.1intx.4, lmm.1intx.5, lmm.1intx.6,
                  lmm.add,
                  lmm.null), AICc)
     
     final.mod = eval(parse(text = paste(rownames(model.sel)[1]))) #best supported model
  
  # check assumptions of best-fit model
  simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T)
  testDispersion(final.mod)
  testZeroInflation(final.mod)
  testCategorical(final.mod, catPred = x$num.sampling[is.na(x$thresh.db)==FALSE])
  
  # estimates from best-supported model
  car::Anova(final.mod, type = "III")
  summary(final.mod)
  ranef(final.mod)
  
  #export model comparison table to file
  # stargazer::stargazer(lmm.full.1.slopes, lmm.full.1.ordinal,
  #                      lmm.3intx.1.slopes, lmm.3intx.2.slopes, lmm.3intx.3.slopes, lmm.3intx.4.slopes,
  #                      lmm.2intx.1.slopes, lmm.2intx.2.slopes, lmm.2intx.3.slopes,
  #                      lmm.1intx.1.slopes, lmm.1intx.2.slopes, lmm.1intx.3.slopes, lmm.1intx.4.slopes, lmm.1intx.5.slopes, lmm.1intx.6.slopes,
  #                      lmm.add.slopes,
  #                      lmm.null.slopes,
  #                      lmm.full.1,
  #                      lmm.3intx.1, lmm.3intx.2, lmm.3intx.3, lmm.3intx.4,
  #                      lmm.2intx.1, lmm.2intx.2, lmm.2intx.3,
  #                      lmm.1intx.1, lmm.1intx.2, lmm.1intx.3, lmm.1intx.4, lmm.1intx.5, lmm.1intx.6,
  #                      lmm.add,
  #                      lmm.null,
  #                      type = "html", out=paste("modelcomparison", "_", deparse(substitute(x)), ".doc", sep=""), intercept.bottom = F, intercept.top = T, digits = 2)
  
  #rename final.mod and model.sel to store it in R working directory
  assign(paste("final_model", deparse(substitute(x)), sep = "_"), final.mod, envir = .GlobalEnv)
  assign(paste("model_compare", deparse(substitute(x)), sep = "_"), model.sel, envir = .GlobalEnv)
  
}

# Create function for juvenile-adult model comparison -----------------------------------------
thresh_model_compare_adultjuv <- function(x){
  #candidate models - manually defined
  #does the effect of stage/age (life.stage.num.sampling) depend on size and frequency?
  lmm.full.1 <- lmer(thresh.db ~ svl.mm*life.stage.num.sampling*poly(freq.hz, 2, raw=FALSE) + (1|combined.id), data = x, na.action = na.omit)
  
  lmm.full.1.ordinal <- lmer(thresh.db ~ svl.mm*life.stage.num.sampling*ordered(freq.hz, levels = unique(x$freq.hz)) + (1|combined.id), data = x, na.action = na.omit)
  
  lmm.full.1 <- lmer(thresh.db ~ svl.mm*life.stage.num.sampling*poly(freq.hz, 2, raw=FALSE) + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of stage/age depend on size independent of frequency?
  lmm.1intx.1 <- lmer(thresh.db ~ svl.mm*life.stage.num.sampling + poly(freq.hz, 2, raw=FALSE)  + (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of stage/age depend on frequency independent of size?
  lmm.1intx.2 <- lmer(thresh.db ~ svl.mm + poly(freq.hz, 2, raw=FALSE) *life.stage.num.sampling+ (1|combined.id), data = x, na.action = na.omit)
  
  #does the effect of size depend on frequency independent of stage/age?
  lmm.1intx.3 <- lmer(thresh.db ~ svl.mm*poly(freq.hz, 2, raw=FALSE)  + life.stage.num.sampling + (1|combined.id), data = x, na.action = na.omit)
  
  #additive effects of stage/age, size, and frequency?
  lmm.add <- lmer(thresh.db ~ svl.mm + life.stage.num.sampling + poly(freq.hz, 2, raw=FALSE)  + (1|combined.id), data = x, na.action = na.omit, REML = TRUE)
  
  #null model
  lmm.null <- lmer(thresh.db ~ (1|combined.id), data = x, na.action = na.omit)
  
  #model comparison using AICc
  model.sel = arrange(AICc(lmm.full.1, lmm.full.1.ordinal,
                           lmm.1intx.1, lmm.1intx.2, lmm.1intx.3,
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
  print(car::Anova(final.mod, type = "III"))
  print(summary(final.mod))
  
  #export model comparison table to file
  # stargazer::stargazer(lmm.full.1, lmm.full.1.ordinal,
  #                      lmm.1intx.1, lmm.1intx.2, lmm.1intx.3,
  #                      lmm.add,
  #                      lmm.null,
  #                      type = "html", out=paste("modelcomparison", "_", deparse(substitute(x)), ".doc", sep=""), intercept.bottom = F, intercept.top = T, digits = 2)
  
  #rename final.mod and model.sel to store it in R working directory
  assign(paste("final_model", deparse(substitute(x)), sep = "_"), final.mod, envir = .GlobalEnv)
  assign(paste("model_compare", deparse(substitute(x)), sep = "_"), model.sel, envir = .GlobalEnv)
  
}


# Supp. Table 1 - juv-only model comparison for vibration -----------------------------------------

#explore fits to data to determine if transformation needed, which influences the structure of the candidate models in function
ggarrange(ggplot(vib.thresh.juv, 
       aes(x = freq.hz, y = thresh.db, group = combined.id, color = combined.id)) +
  facet_wrap(~ num.sampling) +
  #geom_smooth(se = F, span = 1.2) +
  stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = combined.id, group = combined.id)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=combined.id)) +
  geom_smooth(method = "lm", aes(x=freq.hz, y=thresh.db), inherit.aes = F, se = F, color="black"),
  
ggplot(vib.thresh.juv, 
       aes(x = freq.hz^2, y = thresh.db, group = combined.id, color = combined.id)) +
  facet_wrap(~ num.sampling) +
  #geom_smooth(se = F, span = 1.2) +
  stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = combined.id, group = combined.id)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=combined.id)) +  
  geom_smooth(aes(x=freq.hz^2, y=thresh.db), inherit.aes = F, se = F, color="black"), #quadratic

ggplot(vib.thresh.juv, 
       aes(x = freq.hz^3, y = thresh.db, group = combined.id, color = combined.id)) +
  facet_wrap(~ num.sampling) +
  #geom_smooth(se = F, span = 1.2) +
  stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = combined.id, group = combined.id)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=combined.id)) +  
  geom_smooth(aes(x=freq.hz^3, y=thresh.db), inherit.aes = F, se = F, color="black"),

nrow = 3, ncol = 1, common.legend = TRUE
)

thresh_model_compare_juv(vib.thresh.juv)

#pairwise comparison for frequency
pairs(emmeans::emmeans(final_model_vib.thresh.juv, ~num.sampling))
car::Anova(final_model_vib.thresh.juv, type = "III")
model_compare_vib.thresh.juv
summary(final_model_vib.thresh.juv)
vcov(final_model_vib.thresh.juv) 
ranef(final_model_vib.thresh.juv)
coef(final_model_vib.thresh.juv)

MuMIn::r.squaredGLMM(final_model_vib.thresh.juv)
sjPlot::tab_model(final_model_vib.thresh.juv)


# Supp. Table 2 - juv-only model comparison for hearing -----------------------------------------

#explore fits to data to determine if transformation needed, which influences the structure of the candidate models in function
ggarrange(ggplot(hear.thresh.juv, 
                 aes(x = freq.hz, y = thresh.db, group = combined.id, color = combined.id)) +
            facet_wrap(~ num.sampling) +
            #geom_smooth(se = F, span = 1.2) +
            stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = combined.id, group = combined.id)) +
            stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=combined.id)) +
            geom_smooth(method = "lm", aes(x=freq.hz, y=thresh.db), inherit.aes = F, se = F, color="black"),
          
          ggplot(hear.thresh.juv, 
                 aes(x = freq.hz^2, y = thresh.db, group = combined.id, color = combined.id)) +
            facet_wrap(~ num.sampling) +
            #geom_smooth(se = F, span = 1.2) +
            stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = combined.id, group = combined.id)) +
            stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=combined.id)) +  
            geom_smooth(aes(x=freq.hz^2, y=thresh.db), inherit.aes = F, se = F, color="black"), #quadratic
          
          ggplot(hear.thresh.juv, 
                 aes(x = freq.hz^3, y = thresh.db, group = combined.id, color = combined.id)) +
            facet_wrap(~ num.sampling) +
            #geom_smooth(se = F, span = 1.2) +
            stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = combined.id, group = combined.id)) +
            stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=combined.id)) +  
            geom_smooth(aes(x=freq.hz^3, y=thresh.db), inherit.aes = F, se = F, color="black"),
          
          nrow = 3, ncol = 1, common.legend = TRUE
)

thresh_reffects_juv(hear.thresh.juv)
thresh_model_compare_juv(hear.thresh.juv)

#view model comparison table & final model
model_compare_hear.thresh.juv
final_model_hear.thresh.juv

#pairwise comparison for frequency
pairs(emmeans::emmeans(final_model_hear.thresh.juv, ~num.sampling))
car::Anova(final_model_hear.thresh.juv, type = "III")
summary(final_model_hear.thresh.juv)

MuMIn::r.squaredGLMM(final_model_hear.thresh.juv)
sjPlot::tab_model(final_model_hear.thresh.juv)

# Supp. Table 3 - juv-adult model comparison for vibration ------------------
#explore fits to data to determine if transformation needed
ggarrange(ggplot(vib.thresh, 
                 aes(x = freq.hz, y = thresh.db, group = combined.id, color = combined.id)) +
            facet_wrap(~ life.stage.num.sampling) +
            #geom_smooth(se = F, span = 1.2) +
            stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = combined.id, group = combined.id)) +
            stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=combined.id)) +
            geom_smooth(method = "lm", aes(x=freq.hz, y=thresh.db), inherit.aes = F, se = F, color="black"),
          
          ggplot(vib.thresh, 
                 aes(x = freq.hz^2, y = thresh.db, group = combined.id, color = combined.id)) +
            facet_wrap(~ life.stage.num.sampling) +
            #geom_smooth(se = F, span = 1.2) +
            stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = combined.id, group = combined.id)) +
            stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=combined.id)) +  
            geom_smooth(aes(x=freq.hz^2, y=thresh.db), inherit.aes = F, se = F, color="black"), #quadratic
          
          ggplot(vib.thresh, 
                 aes(x = freq.hz^3, y = thresh.db, group = combined.id, color = combined.id)) +
            facet_wrap(~ life.stage.num.sampling) +
            #geom_smooth(se = F, span = 1.2) +
            stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = combined.id, group = combined.id)) +
            stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=combined.id)) +  
            geom_smooth(aes(x=freq.hz^3, y=thresh.db), inherit.aes = F, se = F, color="black"),
          
          nrow = 3, ncol = 1, common.legend = TRUE
)

thresh_model_compare_adultjuv(vib.thresh)

#view model comparison table and final model
model_compare_vib.thresh
final_model_vib.thresh

#pairwise comparison for life.stage.num.sampling
pairs(emmeans::emmeans(final_model_vib.thresh, ~freq.hz, by = "life.stage.num.sampling"))
pairs(emmeans::emmeans(final_model_vib.thresh, ~life.stage.num.sampling, by = "freq.hz"))
car::Anova(final_model_vib.thresh, type = "III")
summary(final_model_vib.thresh)

MuMIn::r.squaredGLMM(final_model_vib.thresh)
sjPlot::tab_model(final_model_vib.thresh)


# Supp. Table 4 - juv-adult model comparison for hearing ------------------
#explore fits to data to determine if transformation needed
ggarrange(ggplot(hear.thresh, 
                 aes(x = freq.hz, y = thresh.db, group = combined.id, color = combined.id)) +
            facet_wrap(~ life.stage.num.sampling) +
            #geom_smooth(se = F, span = 1.2) +
            stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = combined.id, group = combined.id)) +
            stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=combined.id)) +
            geom_smooth(method = "lm", aes(x=freq.hz, y=thresh.db), inherit.aes = F, se = F, color="black"),
          
          ggplot(hear.thresh, 
                 aes(x = freq.hz^2, y = thresh.db, group = combined.id, color = combined.id)) +
            facet_wrap(~ life.stage.num.sampling) +
            #geom_smooth(se = F, span = 1.2) +
            stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = combined.id, group = combined.id)) +
            stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=combined.id)) +  
            geom_smooth(aes(x=freq.hz^2, y=thresh.db), inherit.aes = F, se = F, color="black"), #quadratic
          
          ggplot(hear.thresh, 
                 aes(x = freq.hz^3, y = thresh.db, group = combined.id, color = combined.id)) +
            facet_wrap(~ life.stage.num.sampling) +
            #geom_smooth(se = F, span = 1.2) +
            stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = combined.id, group = combined.id)) +
            stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=combined.id)) +  
            geom_smooth(aes(x=freq.hz^3, y=thresh.db), inherit.aes = F, se = F, color="black"),
          
          nrow = 3, ncol = 1, common.legend = TRUE
)


thresh_model_compare_adultjuv(hear.thresh)

#view model comparison and final model
model_compare_hear.thresh
final_model_hear.thresh

#pairwise comparison for life.stage.num.sampling
pairs(emmeans::emmeans(final_model_hear.thresh, ~freq.hz, by = "life.stage.num.sampling"))
pairs(emmeans::emmeans(final_model_hear.thresh, ~life.stage.num.sampling, by = "freq.hz"))
car::Anova(final_model_hear.thresh, type = "III")
summary(final_model_hear.thresh)

MuMIn::r.squaredGLMM(final_model_hear.thresh)
sjPlot::tab_model(final_model_hear.thresh)

# Figure 1: Juvenile-Adult comparison threshold by life stage and age within life stage (num.sampling) -------

#Create function
thresh_plot <- function(d, y_var, y_var_clip, y_title){
  y_var = sym(y_var)
  y_var_clip = sym(y_var_clip)

  ggplot() +
  
  geom_jitter(data = d,
              aes(y=!!y_var, x = freq.hz, color = life.stage.num.sampling, group = combined.id),
              size = 5, stroke = 1, alpha=0.3, width = 50, height = 0) +
  
  #add points where threshold above clipping level
  geom_jitter(data = d %>%
                filter(is.na(!!y_var) == TRUE), #need to convert units within this line
              fill = "white", pch = 21,
              aes(y=!!y_var_clip, x = freq.hz, color = life.stage.num.sampling, group = combined.id),
              size = 5, stroke = 1, alpha=0.6, width = 50, height = 0) +
  
  stat_summary(data = d,
               fun.y=mean, geom="line", size = 1.2, 
               aes(y=!!y_var, x = freq.hz, color = life.stage.num.sampling, group = life.stage.num.sampling)) +
  
  stat_summary(data = d,
               fun.data=function(d){mean_cl_normal(d, conf.int=.683)}, geom="errorbar", 
               width=0.1, size = 0.9, colour="black", alpha=1, 
               aes(y=!!y_var, x = freq.hz, group = life.stage.num.sampling)) +
  
  stat_summary(data = d,
               fun.y=mean, geom="point", colour = "black", pch=21, size=9, 
               aes(y=!!y_var, x =freq.hz, fill = life.stage.num.sampling, group = life.stage.num.sampling)) +
  
  scale_colour_manual(values = c("#666666","#CC66FF", "#660066", "#330066"), labels = c("adult", "3 months", "6 months", "12 months")) +
  scale_fill_manual(values = c("#666666","#CC66FF", "#660066", "#330066"), labels = c("adult", "3 months", "6 months", "12 months")) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=22, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y=element_text(size=22, color = "black"),
        axis.title.y = element_text(size=22),
        axis.title.x = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "frequency (kHz)", limits = c(50,1550), breaks = sort(unique(d$freq.hz)), labels = sort(unique(d$freq.hz))) +
  scale_y_continuous(name = y_title)
}

png("~/Desktop/R Working Directory/Plots/Figure1.png", units = "in", res = 300, width = 16, height = 12)
ggarrange(thresh_plot(vib.thresh, "rel.db", "clip.thresh.rel.db.plus3", "vibration threshold\n(dB re 1 m/s-2)"), 
          thresh_plot(hear.thresh, "thresh.db", "clip.thresh.db.plus3", "hearing threshold\n(dB re SPL)"), 
          ncol = 1,
          nrow = 2,
          common.legend = TRUE,
          labels = c("a", "b"),
          font.label = list(size = 22, color = "black"))
dev.off()

# Figure 2: JUV ONLY threshold by life stage and age within life stage (num.sampling) -------

thresh_plot_juv <- function(d, y_var, y_var_clip, y_title){
  y_var = sym(y_var)
  y_var_clip = sym(y_var_clip)
  
  fig <- ggplot() +
    
    facet_grid(cols = vars(factor(num.sampling)),
               labeller = as_labeller(c('1' = "3 months", '2' = "6 months", '3' = "12 months"))
    )+
    
    geom_point(data = d,
               aes(y=!!y_var_clip, x = freq.hz, color = factor(num.sampling), group = combined.id),
               size = 3, alpha=0.3, show.legend = FALSE) +
    
    geom_line(data = d,
              aes(y=!!y_var_clip, x = freq.hz, color = factor(num.sampling), group = combined.id),
              size = 1, alpha=0.3, show.legend = FALSE) +
    
    #add points where threshold above clipping level
    geom_point(data = d %>%
                 filter(is.na(!!y_var) == TRUE), #need to convert units within this line
               fill = "white", pch = 21,
               aes(y=!!y_var_clip, x = freq.hz, color = factor(num.sampling), group = combined.id),
               size = 3, stroke = 1, alpha=1.0, show.legend = FALSE) +
    
    # mean points and lines based on non-clipped values
    stat_summary(data = d,
                 fun.y=mean, geom="line", size = 1.2, colour="black",
                 aes(y=!!y_var, x = freq.hz, group = factor(num.sampling)), show.legend = FALSE) +
    
    stat_summary(data = d,
                 fun.data=function(x){mean_cl_normal(x, conf.int=.683)}, geom="errorbar", 
                 width=0.1, size = 0.9, colour="black", alpha=1, 
                 aes(y=!!y_var, x = freq.hz, group = factor(num.sampling)), show.legend = FALSE) +
    
    stat_summary(data = d,
                 fun.y=mean, geom="point", colour = "black", pch=21, size=9, 
                 aes(y=!!y_var, x = freq.hz, fill = factor(num.sampling), group = factor(num.sampling)), show.legend = FALSE) +
    
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
    scale_x_continuous(name = "frequency (kHz)", limits = c(100,2000), 
                       breaks = if(max(d$freq.hz) < 2000){
                         c(sort(unique(d$freq.hz)), 2000)}else{
                           c(100,200,sort(unique(d$freq.hz)))}, 
                       labels = if(max(d$freq.hz) < 2000){
                         c(sort(unique(d$freq.hz)), 2000)}else{
                           c(100,200,sort(unique(d$freq.hz)))}) +
    scale_y_continuous(name = y_title)
  
  # add num.sampling colors to facet grid (code from: https://github.com/tidyverse/ggplot2/issues/2096)
  g.1a <- ggplot_gtable(ggplot_build(fig))
  strip_t <- which(grepl('strip-t', g.1a$layout$name))
  fills <- c("#CC66FF", "#660066", "#330066")
  k <- 1
  for(i in strip_t){
    j <- which(grepl('rect', g.1a$grobs[[i]]$grobs[[1]]$childrenOrder))
    g.1a$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  grid::grid.draw(g.1a)
  
  #rename g.1a and store it in R working directory
  assign(paste("fig", deparse(substitute(d)), sep = "_"), g.1a, envir = .GlobalEnv)
}

thresh_plot_juv(vib.thresh.juv, "rel.db", "clip.thresh.rel.db.plus3", "vibration threshold (dB re 1 m/s-2)")
thresh_plot_juv(hear.thresh.juv, "thresh.db", "clip.thresh.db.plus3", "hearing threshold (dB re SPL)")

png("~/Desktop/R Working Directory/Plots/Figure2.png", units = "in", res = 300, width = 16, height = 12)
ggarrange(fig_vib.thresh.juv,
          fig_hear.thresh.juv,
          ncol = 1, nrow = 2)
dev.off()


# Figure 3:  heatmap: threshold by size and age (num.sampling) -------
#helpful code to make panel plots: https://oscarperpinan.github.io/rastervis/FAQ.html

fig.2a.3 <- levelplot(rel.dB ~ svl.mm * freq.hz, 
                    data = vib.thresh.juv %>% filter(
                      num.sampling == 3),
                    panel = panel.levelplot.points, cex = 1.2,
                    col.regions=colorRampPalette(brewer.pal(9, "Blues"))(25),
                    colorkey=list(at=seq(-105,-30,10), labels=list(cex=1.5, at=seq(-110,-20,10)), title = "vibrational threshold (dB re hearing)"),
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
                    colorkey=list(at=seq(75,135,10), labels=list(cex=1.5, at=seq(75,135,10)), title = "hearing threshold (dB)"),
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


stroke = 1, alpha=0.3, width = 50, height = 0, show.legend = FALSE) +
  
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
  scale_y_continuous(name = "vibration threshold\n(dB re SPL)")




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

# Figure 3 alternate: individual threshold relative to own curve at 3 months (- or + dB from no change at 0) -- easy to see at a glance if individuals get better or worse across frequencies for 6 months and 12 months -------------------
diff_thresh_plot_juv <- function(d, y_var, y_var_clip, y_title){
  y_var = sym(y_var)
  y_var_clip = sym(y_var_clip)
  
  fig <- ggplot() +
    
    geom_hline(yintercept = 0, color = "gray15", size = 1, linetype = 1) + #add horizontal line at 0 change
    
    geom_point(data = d %>%
                 filter(num.sampling > 1),
               aes(y=!!y_var_clip, x = freq.hz, color = factor(num.sampling), group = id.life.stage.num.sampling),
               size = 3, alpha=0.3, show.legend = FALSE) +
    
    geom_line(data = d %>%
                filter(num.sampling > 1),
              aes(y=!!y_var_clip, x = freq.hz, color = factor(num.sampling), group = id.life.stage.num.sampling),
              size = 1, alpha=0.3, show.legend = FALSE) +
    
    #add points where threshold above clipping level
    geom_point(data = d %>%
                 filter(num.sampling > 1) %>%
                 filter(is.na(!!y_var) == TRUE), #need to convert units within this line
               fill = "white", pch = 21,
               aes(y=!!y_var_clip, x = freq.hz, color = factor(num.sampling), group = id.life.stage.num.sampling),
               size = 3, stroke = 1, alpha=1.0, show.legend = FALSE) +
    
    # mean points and lines based on non-clipped values
    stat_summary(data = d %>%
                   filter(num.sampling > 1),
                 fun.y=mean, geom="line", size = 1.2, colour="black",
                 aes(y=!!y_var, x = freq.hz, group = factor(num.sampling)), show.legend = FALSE) +

    stat_summary(data = d %>%
                   filter(num.sampling > 1),
                 fun.data=function(x){mean_cl_normal(x, conf.int=.683)}, geom="errorbar",
                 width=0.1, size = 0.9, colour="black", alpha=1,
                 aes(y=!!y_var, x = freq.hz, group = factor(num.sampling)), show.legend = FALSE) +

    stat_summary(data = d %>%
                   filter(num.sampling > 1),
                 fun.y=mean, geom="point", colour = "black", pch=21, size=9,
                 aes(y=!!y_var, x = freq.hz, fill = factor(num.sampling), group = factor(num.sampling)), show.legend = TRUE) +

    scale_colour_manual(values = c("#660066", "#330066"), labels = c("6 months", "12 months")) +
    scale_fill_manual(values = c("#660066", "#330066"), labels = c("6 months", "12 months")) +
    
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
    scale_x_continuous(name = "frequency (Hz)", limits = c(100,2000), 
                       breaks = if(max(d$freq.hz) < 2000){
                         c(sort(unique(d$freq.hz)), 2000)}else{
                           c(100,200,sort(unique(d$freq.hz)))}, 
                       labels = if(max(d$freq.hz) < 2000){
                         c(sort(unique(d$freq.hz)), 2000)}else{
                           c(100,200,sort(unique(d$freq.hz)))}) +
    scale_y_continuous(name = y_title, limits = c(-35,25), breaks = seq(-35,25,5), labels = seq(-35,25,5))
}

png("~/Desktop/R Working Directory/Plots/Figure3.png", units = "in", res = 300, width = 16, height = 12)
ggarrange(
  diff_thresh_plot_juv(vib.thresh.juv.diff, "thresh.db.diff", "clip.thresh.db.diff", "vibration threshold relative to 3 months"),
  diff_thresh_plot_juv(hear.thresh.juv.diff, "thresh.db.diff", "clip.thresh.db.diff", "hearing threshold relative to 3 months"),
  common.legend = TRUE,
  labels = c("a", "b"),
  font.label = list(size = 20, color = "black")
)
dev.off()
