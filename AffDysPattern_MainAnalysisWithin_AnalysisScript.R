wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(wd, "/Data/Preprocessed_Data"))

df <- read.csv("EmoPattern_neuralData.csv")



#################################################
# load packages and functions

library("psych")
library("ggplot2")
library("effectsize")
library("lme4")
library("plyr")
library("meta")
library("metafor")
library("BayesFactor")
library("brms")
library("nlme")
library("afex")
library("corrplot")
library("reshape")
library("caret")

setwd(wd)
source("AffDysPattern_functions.R")


#################################################
# prepare variables

# calculate contrasts "negative vs neutral"
df$Pines_Diff <- df$NegPines - df$NeutPines
df$Kragel_Diff_Amused <- df$NegViewKragelAmused - df$NeutViewKragelAmused
df$Kragel_Diff_Angry <- df$NegViewKragelAngry - df$NeutViewKragelAngry
df$Kragel_Diff_Content <- df$NegViewKragelContent - df$NeutViewKragelContent
df$Kragel_Diff_Fearful <- df$NegViewKragelFearful - df$NeutViewKragelFearful
df$Kragel_Diff_Neutral <- df$NegViewKragelNeutral - df$NeutViewKragelNeutral
df$Kragel_Diff_Sad <- df$NegViewKragelSad - df$NeutViewKragelSad
df$Kragel_Diff_Surprised <- df$NegViewKragelSurprised - df$NeutViewKragelSurprised

# get variable names of difference scores for all patterns
pattern_Diff_Names <- names(df)[which(names(df) == "Pines_Diff"):which(names(df) == "Kragel_Diff_Surprised")]

# ID as factor
df$ID <- as.factor(df$ID)

# remove remitted and trauma control groups and recode group
df <- df[df$group != "TC" & df$group != "BPDremit", ]
df$group <- ifelse(df$group == "HC", "Healthy", "Clinical")

# emotion labels
emotion_label <- c("Negative Affect", "Amused", "Angry", "Content", "Fearful", "Neutral", "Sad", "Surprised")



#################################################
#  within mega-analysis (negative vs neutral)


#########################
# prepare variables

# scale pattern difference scores to SD (without centering!)
for(i in 1:length(pattern_Diff_Names)){
  df[, paste0(pattern_Diff_Names[i], "_scaled")] <- scaleWithinRuns(pattern_Diff_Names[i], design = "within")
}


#########################
# Conduct frequentist mega-analysis

# prepare results df for mega analysis
resultsMEGAwithin <- as.data.frame(
  matrix(
    nrow = length(pattern_Diff_Names), 
    ncol = 4))
names(resultsMEGAwithin) <- c("pattern", "estimate", "lb", "ub")


# conduct mega analysis on all patterns
for(i in 1:length(pattern_Diff_Names)){
  
  dv <- paste0(pattern_Diff_Names[i], "_scaled")
  f <- as.formula(paste(dv, "1 + (1|study:ID)", sep = " ~ "))
  
  MEGAwithinModel <- lmer(
    data = df,
    formula = f
  )
  
  cat("__________________________________________________________________________________\n\n", 
      "PATTERN: ", pattern_Diff_Names[i], "\n\n\n")
  print(summary(MEGAwithinModel))
  print(round(confint(MEGAwithinModel), 2))
  
  est <- round(summary(MEGAwithinModel)$coefficients[1], 2)
  lb <- round(confint(MEGAwithinModel)[3, 1], 2)
  ub <- round(confint(MEGAwithinModel)[3, 2], 2)
  
  resultsMEGAwithin[i, ] <- list(pattern_Diff_Names[i], est, lb, ub)
  
}

# show  results
resultsMEGAwithin
mean(abs(resultsMEGAwithin$estimate))

# plot results
setwd(paste0(wd, "/Figures"))

dfplotMEGAwithin <- resultsMEGAwithin
dfplotMEGAwithin$emotion_label <- factor(emotion_label, levels = c("Negative Affect", "Fearful", "Angry", "Sad", "Amused", "Content", "Surprised", "Neutral"))

ggplot(data = dfplotMEGAwithin, aes(y = estimate, x = emotion_label)) +
  geom_point(size = 3) + geom_errorbar(aes(ymin = lb, ymax = ub)) +
  
  ylab(expression(atop("Pattern Expression", paste("[Cohen's ", italic("d"), ", Negative - Neutral]")))) + xlab(NULL) +
  
  theme_classic() +
  theme(axis.line.x=element_blank(), axis.ticks = element_blank(), axis.title.x=element_blank()) + 
  geom_segment(aes(x=0,y=0,xend=9,yend=0)) 

ggsave("PlotMEGAwithin.png", device = "png")
ggsave("PlotMEGAwithin.pdf", device = "pdf")






#########################
# calculate Bayes Factors

# [the following code is commented out, as it takes ~20 minutes to run and sometimes breaks (unsystematically!) the R session.
# Therefore, models were saved as external objects.
# Also, unfortunately, seeds could not be set due to errors when running the H1 model with "seed = ...", limiting exact reproducibility]

# directions of pattern-specific hypotheses
# [higher: H1 = Negative > Neutral]
# patternHypothesesWithin <- c("higher", "lower", "higher", "lower", "higher", "lower", "higher", "higher")
# 
# setwd(paste0(wd, "/BayesModels"))
# for(i in 1:length(pattern_Diff_Names)){
# 
#    # set one-sided priors
#    if(patternHypothesesWithin[i] == "higher"){
#      cauchyPrior <- prior("cauchy(0, 0.707)", class = b, lb = 0)
#    }else{
#      cauchyPrior <- prior("cauchy(0, 0.707)", class = b, ub = 0)
#    }
# 
#    # H1 model
#    dv <- paste0(pattern_Diff_Names[i], "_scaled")
#    f1 <- as.formula(paste(dv, "0 + intercept + (0+intercept|study:ID)", sep = " ~ "))
# 
#    MEGAwithinModelBayes <- brms::brm(
#      data = df,
#      formula = f1,
#      save_all_pars = T, prior = cauchyPrior, sample_prior = "yes", cores = 3,
#      control = list(adapt_delta = 0.90, max_treedepth = 15),
#      file = paste0(pattern_Diff_Names[i], "_H1")
#    )
# 
#    # H0 model
#    dv <- paste0(pattern_Diff_Names[i], "_scaled")
#    f0 <- as.formula(paste(dv, "0 + (1|study:ID)", sep = " ~ "))
# 
#    MEGAwithinModelBayes_null <- brms::brm(
#      data=df,
#      formula = f0,
#      save_all_pars = T, sample_prior = "yes", cores = 3,
#      control = list(adapt_delta = 0.90, max_treedepth = 15),
#      file = paste0(pattern_Diff_Names[i], "_H0"))
# 
# }





# compute Bayesfactors from Bayesian models

# prepare loop
setwd(paste0(wd, "/BayesModels"))
BFresultsWithin <- as.data.frame(
  matrix(nrow = length(pattern_Diff_Names),
         ncol = 2,
         dimnames = list(NULL, c("pattern", "BF01"))),
)

# print diagnostics for H1 models
getH1modelDiagnostics <- T

# loop through emotion patterns and compute BF
for(i in 1:length(pattern_Diff_Names)){
  
  #load models
  H0 <- readRDS(paste0(pattern_Diff_Names[i], "_H0.rds"))
  H1 <- readRDS(paste0(pattern_Diff_Names[i], "_H1.rds"))
  
  #compute and save BF01
  set.seed(i)
  bf01 <- round(brms::bayes_factor(H0, H1)$bf, 2)
  BFresultsWithin[i, ] <- list(pattern_Diff_Names[i], bf01)
  
  # [optional] plot diagnostics for H1 model
  if(getH1modelDiagnostics){
    print(plot(H1))
    div_errors <- nuts_params(H1)
    print(pairs(H1, np = div_errors))
    print(pp_check(H1))
  }
}

# show results
BFresultsWithin$BF10 <- round(1/BFresultsWithin$BF01, 2)
BFresultsWithin$emotion_label <- emotion_label
BFresultsWithin




#########################
# calculate for single studies/runs

# make dataframe with data for study 1 averaged over runs
dfemoregAgg <- ddply(df[df$study == "emoreg", ], .(ID, group), numcolwise(mean))
dfemoregAgg$run <- rep(1, nrow(dfemoregAgg))
dfemoregAgg$rundescr <- rep("average", nrow(dfemoregAgg))
dfemoregAgg$study <- rep("emoreg", nrow(dfemoregAgg))
dfsingleStudies <- rbind(df[df$study != "emoreg", ], dfemoregAgg[, names(df)])

# hypotheses directions
patternHypothesesWithin <- c("higher", "lower", "higher", "lower", "higher", "lower", "higher", "higher")

# prepare loop and results df
studiesRuns <- unique(dfsingleStudies[, c("study", "run")])
studiesRunsResults <- as.data.frame(
  matrix(nrow = nrow(studiesRuns)*length(pattern_Diff_Names), 
         ncol = 7,
         dimnames = list(NULL, c("pattern", "study", "run", "cohens_d", "lb", "ub", "BF10"))
  )
)
counter = 0

# loop through studies/runs
for(i in 1:nrow(studiesRuns)){
  
  # subset study/run
  dfsingleStudiesSingle <- dfsingleStudies[dfsingleStudies$study == studiesRuns$study[i] & dfsingleStudies$run == studiesRuns$run[i], ]
  
  # loop through emotion patterns
  for(j in 1:length(pattern_Diff_Names)){
    
    counter <- counter + 1
    
    if(patternHypothesesWithin[j] == "higher"){
      alternative = c(0, Inf)
    }else if(patternHypothesesWithin[j] == "lower"){
      alternative = c(-Inf, 0)
    }else{
      alternative = c(-Inf, Inf)
    }
    
    cd <- round(cohens_d(dfsingleStudiesSingle[, pattern_Diff_Names[j]])$Cohens_d, 2)
    lb <- round(cohens_d(dfsingleStudiesSingle[, pattern_Diff_Names[j]])$CI_low, 2)
    ub <- round(cohens_d(dfsingleStudiesSingle[, pattern_Diff_Names[j]])$CI_high, 2)
    BF <- round(extractBF(ttestBF(dfsingleStudiesSingle[, pattern_Diff_Names[j]], 
                                  nullInterval=alternative))[1, "bf"], 2)
    
    results <- c(pattern_Diff_Names[j], 
                 studiesRuns$study[i], 
                 studiesRuns$run[i],
                 cd,
                 lb,
                 ub,
                 BF)
    
    studiesRunsResults[counter, ] <- results
    
  }
  
}

studiesRunsResults$BF10 <- ifelse(as.numeric(studiesRunsResults$BF10) > 100, ">100", studiesRunsResults$BF10)



# average results for study 2 with mixed model

# short df
dfreact <- df[df$study == "react", ]

# prepare loop and results df
studiesRunsResultsReact <- as.data.frame(
  matrix(nrow = length(pattern_Diff_Names), 
         ncol = 7,
         dimnames = list(NULL, c("pattern", "study", "run", "cohens_d", "lb", "ub", "BF10"))
  )
)

for(i in 1:length(pattern_Diff_Names)){
  
  dv <- paste0(pattern_Diff_Names[i], "_scaled")
  f <- as.formula(paste(dv, "1 + (1|ID)", sep = " ~ "))
  
  model.react <- lmer(formula = f, data = dfreact)
  cd <- round(summary(model.react)$coefficients[1], 2)
  lb <- round(confint(model.react)[3, 1], 2)
  ub <- round(confint(model.react)[3, 2], 2)
  
  
  results <- c(pattern_Diff_Names[i], 
               "react", 
               1,
               cd,
               lb,
               ub,
               NA)
  
  studiesRunsResultsReact[i, ] <- results
}

studiesRunsResultsReact

# merge results with results from the two other studies
dfMainResultsWithin <- rbind(studiesRunsResults[studiesRunsResults$study != "react", ], studiesRunsResultsReact)

dfMainResultsWithin$cohens_d <- as.numeric(dfMainResultsWithin$cohens_d)
dfMainResultsWithin$lb <- as.numeric(dfMainResultsWithin$lb)
dfMainResultsWithin$ub <- as.numeric(dfMainResultsWithin$ub)
dfMainResultsWithin$emotion_label <- factor(rep(c("Negative Affect",  "Amused", "Angry", "Content", "Fearful", "Neutral", "Sad", "Surprised"), 3), 
                                             levels = c("Negative Affect", "Fearful", "Angry", "Sad", "Amused", "Content", "Surprised", "Neutral"))
dfMainResultsWithin <- dfMainResultsWithin[, -c(3, 7)]
dfMainResultsWithin$study <- ifelse(dfMainResultsWithin$study == "emoreg", "Study 1",
                                     ifelse(dfMainResultsWithin$study == "react", "Study 2",
                                            "Study 3"))


# prepare mega-analysis df
resultsMEGAwithin$emotion_label <- factor(c("Negative Affect",  "Amused", "Angry", "Content", "Fearful", "Neutral", "Sad", "Surprised"), 
                                           levels = c("Negative Affect", "Fearful", "Angry", "Sad", "Amused", "Content", "Surprised", "Neutral"))
names(resultsMEGAwithin)[2] <- "cohens_d"


setwd(paste0(wd, "/Figures"))

ggplot() +
  
  geom_point(data = resultsMEGAwithin, aes(y = cohens_d, x = emotion_label), size = 4, colour = "white") +
  
  geom_segment(aes(x=0,y=0,xend=9,yend=0)) +
  
  geom_errorbar(data = dfMainResultsWithin, aes(y = cohens_d, x = emotion_label, colour = study, ymin = lb, ymax = ub), position=position_dodge(width=0.7), alpha = 0.20) + 
  geom_point(data = dfMainResultsWithin, aes(y = cohens_d, x = emotion_label, colour = study), position=position_dodge(width=0.7), alpha = 0.20) + 
  geom_errorbar(data = resultsMEGAwithin, aes(y = cohens_d, x = emotion_label, ymin = lb, ymax = ub), size = 1.5, colour = "white") +
  geom_point(data = resultsMEGAwithin, aes(y = cohens_d, x = emotion_label), size = 4, colour = "white") +
  geom_errorbar(data = resultsMEGAwithin, aes(y = cohens_d, x = emotion_label, ymin = lb, ymax = ub)) +
  geom_point(data = resultsMEGAwithin, aes(y = cohens_d, x = emotion_label), size = 3) +
  
  ylab(expression(atop("Pattern Expression", paste("[Cohen's ", italic("d"), ", Negative  - Neutral]")))) + xlab(NULL) +
  ylim(-3.5, 4.5) +
  
  theme_classic() +
  theme(axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle = 20), legend.title = element_blank()) +
  scale_y_continuous(breaks=c(seq(from=-3, to = 4, by = 1)))


ggsave("plotMEGAWithin.png", device = "png")
ggsave("plotMEGAWithin.pdf", device = "pdf")



#################################################
# Accuracy


PinesAcc <- function(data, study, run){
  
  # prepare classification df
  dfsubset <- data[data$study == study & data$run == run, ]
  dfsubsetMelt <- melt(dfsubset, id.vars = "ID", measure.vars = c("NegPines", "NeutPines"))
  
  # calculate logistic regression
  logModel <- glm(variable ~ value, data = dfsubsetMelt, family = "binomial")
  
  # accuracy metrics
  logModel_pred <- ifelse(predict(logModel, type = "link") > 0, "NeutPines", "NegPines")
  train_tab <- table(predicted = logModel_pred, actual = dfsubsetMelt$variable)
  confusionMatrix(train_tab, positive = "NegPines")
}

PinesAcc(dfsingleStudies, "emoreg", 1)
PinesAcc(dfsingleStudies, "react", 1)
PinesAcc(dfsingleStudies, "react", 2)
PinesAcc(dfsingleStudies, "react", 3)
PinesAcc(dfsingleStudies, "ewmt", 1)



# logistic regression
ewmtModel <- glm(variable ~ value, data = dfewmtMelt, family = "binomial")
summary(ewmtModel)

# accuracy
ewmtModel_pred <- ifelse(predict(ewmtModel, type = "link") > 0, "NeutPines", "NegPines")
train_tab <- table(predicted = ewmtModel_pred, actual = dfewmtMelt$variable)
confusionMatrix(train_tab, positive = "NegPines")

