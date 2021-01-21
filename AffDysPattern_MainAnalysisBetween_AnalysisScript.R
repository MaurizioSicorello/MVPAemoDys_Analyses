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
library("stringr")

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


# calculate sample sizes
Nvars <- c("ID", "group", "study")
nrow(unique(df[, Nvars]))
nrow(unique(df[df$group == "Clinical", Nvars]))
nrow(unique(df[df$study == "emoreg", Nvars]))
nrow(unique(df[df$study == "emoreg" & df$group == "Clinical", Nvars]))
nrow(unique(df[df$study == "react", Nvars]))
nrow(unique(df[df$study == "react" & df$group == "Clinical", Nvars]))
nrow(unique(df[df$study == "ewmt", Nvars]))
nrow(unique(df[df$study == "ewmt" & df$group == "Clinical", Nvars]))


# correlation of responsivity between different emotions INVALID BECAUSE OF UNMODELLED DEPENDENCIES!
dfCorr <- df[, which(names(df) == "Pines_Diff"):which(names(df) == "Kragel_Diff_Surprised")]
signatureCorrs <- cor(dfCorr)
corrplot(signatureCorrs, order = "hclust", addrect = 3)

pca.results <- prcomp(signatureCorrs)
varimax.results <- varimax(signatureCorrs)
screeplot(pca.results)


#################################################
#  between mega-analysis (psychopathology vs healthy)


#########################
# prepare variables


# standardize pattern difference scores
for(i in 1:length(pattern_Diff_Names)){
  df[, paste0(pattern_Diff_Names[i], "_scalePool")] <- scaleWithinRuns(pattern_Diff_Names[i], design = "between", group_var = "group")
}

# weighted regression effect coding for group 
df <- ddply(df, .(study, run), transform, group_weightEff = weightedEffectCoding(group))



#########################
# Conduct frequentist mega-analysis

# validate that fixed intercept and intercept variance between studies equal zero
summary(
  lmer(
    data = df,
    Pines_Diff_scalePool ~ 1 + (1|study:ID) + (1|study)
  )
)


# prepare results df for mega analysis
resultsMEGAbetween <- as.data.frame(
  matrix(
    nrow = length(pattern_Diff_Names),
    ncol = 4))
names(resultsMEGAbetween) <- c("pattern", "estimate", "lb", "ub")


# conduct mega analysis on all patterns
for(i in 1:length(pattern_Diff_Names)){

  dv <- paste0(pattern_Diff_Names[i], "_scalePool")
  f <- as.formula(paste(dv, "0 + group_weightEff + (1|study:ID)", sep = " ~ "))

  MEGAbetweenModel <- lmer(
    data = df,
    formula = f
  )

  cat("__________________________________________________________________________________\n\n",
      "PATTERN: ", pattern_Diff_Names[i], "\n\n\n")
  print(summary(MEGAbetweenModel))
  print(round(confint(MEGAbetweenModel), 2))

  est <- round(summary(MEGAbetweenModel)$coefficients[1], 2)
  lb <- round(confint(MEGAbetweenModel)[3, 1], 2)
  ub <- round(confint(MEGAbetweenModel)[3, 2], 2)

  resultsMEGAbetween[i, ] <- list(pattern_Diff_Names[i], est, lb, ub)

}

# show  results
resultsMEGAbetween

# plot results
setwd(paste0(wd, "/Figures"))

dfplotMEGAbetween <- resultsMEGAbetween
dfplotMEGAbetween$emotion_label <- factor(c("Negative Affect",  "Amused", "Angry", "Content", "Fearful", "Neutral", "Sad", "Surprised"), 
                                          levels = c("Negative Affect", "Fearful", "Angry", "Sad", "Amused", "Content", "Surprised", "Neutral"))

ggplot(data = dfplotMEGAbetween, aes(y = estimate, x = emotion_label)) +
  geom_point(size = 3) + geom_errorbar(aes(ymin = lb, ymax = ub)) +
  
  ylab(expression(atop("Pattern Expression", paste("[Cohen's ", italic("d"), ", Emotion Dysregulation - Healthy Controls]")))) + xlab(NULL) +
  ylim(-1, 1) +
  
  theme_classic() +
  theme(axis.line.x=element_blank(), axis.ticks = element_blank(), axis.title.x=element_blank()) + 
  geom_segment(aes(x=0,y=0,xend=9,yend=0)) 

ggsave("PlotMEGAbetween.png", device = "png")
ggsave("PlotMEGAbetween.pdf", device = "pdf")




#########################
# calculate Bayes Factors

# [the following code is commented out, as it takes ~20 minutes to run and sometimes breaks (unsystematically!) the R session.
# Therefore, models were saved as external objects.
# Also, unfortunately, seeds could not be set due to errors when running the H1 model with "seed = ...", limiting exact reproducibility]

# directions of pattern-specific hypotheses
# [higher: H1 = Negative > Neutral]
# patternHypothesesbetween <- c("higher", "lower", "higher", "lower", "higher", "lower", "higher", "unequal")
#  
# setwd(paste0(wd, "/BayesModels"))
# for(i in 1:length(pattern_Diff_Names)){
# 
#   # set one-sided priors
#   if(patternHypothesesbetween[i] == "higher"){
#     cauchyPrior <- prior("cauchy(0, 0.707)", class = b, lb = 0)
#   }else if(patternHypothesesbetween[i] == "lower"){
#     cauchyPrior <- prior("cauchy(0, 0.707)", class = b, ub = 0)
#   }else{
#     cauchyPrior <- prior("cauchy(0, 0.707)", class = b)
#   }
# 
#   # H1 model
#   dv <- paste0(pattern_Diff_Names[i], "_scalePool")
#   f1 <- as.formula(paste(dv, "0 + group_weightEff + (1|study:ID)", sep = " ~ "))
# 
#   MEGAbetweenModelBayes <- brms::brm(
#     data = df,
#     formula = f1,
#     save_all_pars = T, prior = cauchyPrior, sample_prior = "yes", cores = 3,
#     control = list(adapt_delta = 0.99, max_treedepth = 15),
#     file = paste0(pattern_Diff_Names[i], "_between_H1")
#   )
# 
#   # H0 model
#   dv <- paste0(pattern_Diff_Names[i], "_scalePool")
#   f0 <- as.formula(paste(dv, "0 + (1|study:ID)", sep = " ~ "))
# 
#   MEGAbetweenModelBayes_null <- brms::brm(
#     data=df,
#     formula = f0,
#     save_all_pars = T, sample_prior = "yes", cores = 3,
#     control = list(adapt_delta = 0.90, max_treedepth = 15),
#     file = paste0(pattern_Diff_Names[i], "_between_H0"))
# 
# }



# compute Bayesfactors from Bayesian models
setwd(paste0(wd, "/BayesModels"))
BFresultsBetween <- as.data.frame(
  matrix(nrow = length(pattern_Diff_Names),
         ncol = 2,
         dimnames = list(NULL, c("pattern", "BF01"))),
)

getH1modelDiagnostics <- T

for(i in 1:length(pattern_Diff_Names)){
  
  #load models
  H0 <- readRDS(paste0(pattern_Diff_Names[i], "_between_H0.rds"))
  H1 <- readRDS(paste0(pattern_Diff_Names[i], "_between_H1.rds"))
  
  #compute and save BF01
  set.seed(i)
  bf01 <- round(brms::bayes_factor(H0, H1)$bf, 2)
  BFresultsBetween[i, ] <- list(pattern_Diff_Names[i], bf01)
  
  # [optional] plot diagnostics for H1 model
  if(getH1modelDiagnostics){
    print(plot(H1))
    div_errors <- nuts_params(H1)
    print(pairs(H1, np = div_errors))
    print(pp_check(H1))
  }
}

# show results
BFresultsBetween



#########################
# calculate for single studies/runs

# make dataframe with data for study 1 averaged over runs
dfemoregAgg <- ddply(df[df$study == "emoreg", ], .(ID, group), numcolwise(mean))
dfemoregAgg$run <- rep(1, nrow(dfemoregAgg))
dfemoregAgg$rundescr <- rep("average", nrow(dfemoregAgg))
dfemoregAgg$study <- rep("emoreg", nrow(dfemoregAgg))
dfsingleStudies <- rbind(df[df$study != "emoreg", ], dfemoregAgg[, names(df)])


# prepare loop and results df
studiesRuns <- unique(dfsingleStudies[, c("study", "run")])
studiesRunsResults <- as.data.frame(
  matrix(nrow = nrow(studiesRuns)*length(pattern_Diff_Names), 
         ncol = 7,
         dimnames = list(NULL, c("pattern", "study", "run", "cohens_d", "lb", "ub", "BF01"))
  )
)

counter = 0

patternHypothesesbetween <- c("higher", "lower", "higher", "lower", "higher", "lower", "higher", "unequal")

# loop through studies/runs
for(i in 1:nrow(studiesRuns)){
  
  # subset study/run
  dfsingleStudiesSingle <- dfsingleStudies[dfsingleStudies$study == studiesRuns$study[i] & dfsingleStudies$run == studiesRuns$run[i], ]
  
  # loop through emotion patterns
  for(j in 1:length(pattern_Diff_Names)){
    
    counter <- counter + 1
    
    if(patternHypothesesbetween[j] == "higher"){
      alternative = c(0, Inf)
    }else if(patternHypothesesbetween[j] == "lower"){
        alternative = c(-Inf, 0)
      }else{
        alternative = c(-Inf, Inf)
      }
    
    
    BF <- round(extractBF(ttestBF(dfsingleStudiesSingle[dfsingleStudiesSingle$group == "Clinical", pattern_Diff_Names[j]], 
                                  dfsingleStudiesSingle[dfsingleStudiesSingle$group == "Healthy", pattern_Diff_Names[j]],
                                  nullInterval=alternative))[1, "bf"], 2)
    BF <- round(1/BF, 2)
    cd_model <- cohens_d(dfsingleStudiesSingle[dfsingleStudiesSingle$group == "Clinical", pattern_Diff_Names[j]], 
                         dfsingleStudiesSingle[dfsingleStudiesSingle$group == "Healthy", pattern_Diff_Names[j]])
    cd <- as.numeric(round(cd_model$Cohens_d, 2))
    lb <- round(cd_model$CI_low, 2)
    ub <- round(cd_model$CI_high, 2)
    
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

studiesRunsResults


# average results for study 2 with mixed model

# short df
dfreact <- df[df$study == "react", ]

# prepare loop and results df
studiesRunsResultsReact <- as.data.frame(
  matrix(nrow = length(pattern_Diff_Names), 
         ncol = 7,
         dimnames = list(NULL, c("pattern", "study", "run", "cohens_d", "lb", "ub", "BF01"))
  )
)

for(i in 1:length(pattern_Diff_Names)){
  
  dv <- paste0(pattern_Diff_Names[i], "_scalePool")
  f <- as.formula(paste(dv, "0 + group_weightEff + (1|ID)", sep = " ~ "))
  
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
dfMainResultsBetween <- rbind(studiesRunsResults[studiesRunsResults$study != "react", ], studiesRunsResultsReact)

dfMainResultsBetween$cohens_d <- as.numeric(dfMainResultsBetween$cohens_d)
dfMainResultsBetween$lb <- as.numeric(dfMainResultsBetween$lb)
dfMainResultsBetween$ub <- as.numeric(dfMainResultsBetween$ub)
dfMainResultsBetween$emotion_label <- factor(rep(c("Negative Affect",  "Amused", "Angry", "Content", "Fearful", "Neutral", "Sad", "Surprised"), 3), 
                                          levels = c("Negative Affect", "Fearful", "Angry", "Sad", "Amused", "Content", "Surprised", "Neutral"))
dfMainResultsBetween <- dfMainResultsBetween[, -c(3, 7)]
dfMainResultsBetween$study <- ifelse(dfMainResultsBetween$study == "emoreg", "Study 1",
                                     ifelse(dfMainResultsBetween$study == "react", "Study 2",
                                            "Study 3"))


# prepare mega-analysis df
resultsMEGAbetween$emotion_label <- factor(c("Negative Affect",  "Amused", "Angry", "Content", "Fearful", "Neutral", "Sad", "Surprised"), 
                                             levels = c("Negative Affect", "Fearful", "Angry", "Sad", "Amused", "Content", "Surprised", "Neutral"))
names(resultsMEGAbetween)[2] <- "cohens_d"


setwd(paste0(wd, "/Figures"))

ggplot() +
  
  geom_point(data = resultsMEGAbetween, aes(y = cohens_d, x = emotion_label), size = 4, colour = "white") +
  
  geom_segment(aes(x=0,y=0,xend=9,yend=0)) +
  
  geom_errorbar(data = dfMainResultsBetween, aes(y = cohens_d, x = emotion_label, colour = study, ymin = lb, ymax = ub), position=position_dodge(width=0.7), alpha = 0.20) + 
  geom_point(data = dfMainResultsBetween, aes(y = cohens_d, x = emotion_label, colour = study), position=position_dodge(width=0.7), alpha = 0.20) + 
  geom_errorbar(data = resultsMEGAbetween, aes(y = cohens_d, x = emotion_label, ymin = lb, ymax = ub), size = 1.5, colour = "white") +
  geom_point(data = resultsMEGAbetween, aes(y = cohens_d, x = emotion_label), size = 4, colour = "white") +
  geom_errorbar(data = resultsMEGAbetween, aes(y = cohens_d, x = emotion_label, ymin = lb, ymax = ub)) +
  geom_point(data = resultsMEGAbetween, aes(y = cohens_d, x = emotion_label), size = 3) +
  
  ylab(expression(atop("Pattern Expression", paste("[Cohen's ", italic("d"), ", Emotion Dysregulation - Healthy Controls]")))) + xlab(NULL) +
  
  theme_classic() +
  theme(axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle = 20), legend.title = element_blank()) +
  scale_y_continuous(limits = c(-1.5, 1.5), breaks=c(seq(from=-1.5, to = 1.5, by = 0.5)))

ggsave("plotMEGAbetween.png", device = "png")
ggsave("plotMEGAbetween.pdf", device = "pdf")


#################################################
#  between mega-analysis with medication as covariate

# recode medication status

# prepare results df for mega analysis
resultsMEGAbetweenMed <- as.data.frame(
  matrix(
    nrow = length(pattern_Diff_Names),
    ncol = 4))
names(resultsMEGAbetweenMed) <- c("pattern", "estimate", "lb", "ub")


# conduct mega analysis on all patterns
for(i in 1:length(pattern_Diff_Names)){
  
  dv <- paste0(pattern_Diff_Names[i], "_scalePool")
  f <- as.formula(paste(dv, "1 + group_weightEff + med_bin + (1|study:ID)", sep = " ~ "))
  
  MEGAbetweenMedModel <- lmer(
    data = df,
    formula = f
  )
  
  cat("__________________________________________________________________________________\n\n",
      "PATTERN: ", pattern_Diff_Names[i], "\n\n\n")
  print(summary(MEGAbetweenMedModel))
  print(round(confint(MEGAbetweenMedModel), 2))
  
  est <- round(summary(MEGAbetweenMedModel)$coefficients[2], 2)
  lb <- round(confint(MEGAbetweenMedModel)[4, 1], 2)
  ub <- round(confint(MEGAbetweenMedModel)[4, 2], 2)
  
  resultsMEGAbetweenMed[i, ] <- list(pattern_Diff_Names[i], est, lb, ub)
  
}

# show  results
resultsMEGAbetweenMed

# plot results
setwd(paste0(wd, "/Figures"))

dfplotMEGAbetweenMed <- resultsMEGAbetweenMed
dfplotMEGAbetweenMed$emotion_label <- factor(c("Negative Affect",  "Amused", "Angry", "Content", "Fearful", "Neutral", "Sad", "Surprised"), 
                                          levels = c("Negative Affect", "Fearful", "Angry", "Sad", "Amused", "Content", "Surprised", "Neutral"))

ggplot(data = dfplotMEGAbetweenMed, aes(y = estimate, x = emotion_label)) +
  geom_point(size = 3) + geom_errorbar(aes(ymin = lb, ymax = ub)) +
  
  ylab(expression(atop("Pattern Expression (medication controlled)", paste("[Cohen's ", italic("d"), ", Emotion Dysregulation - Healthy Controls]")))) + xlab(NULL) +
  ylim(-1, 1) +
  
  theme_classic() +
  theme(axis.line.x=element_blank(), axis.ticks = element_blank(), axis.title.x=element_blank()) + 
  geom_segment(aes(x=0,y=0,xend=9,yend=0)) 

ggsave("PlotMEGAbetweenMed.png", device = "png")
ggsave("PlotMEGAbetweenMed.pdf", device = "pdf")



#################################################
#  between mega-analysis on neutral baseline


#########################
# prepare variables

# get names of variables for pattern expression in neutral condition
pattern_Neut_Names <- c("NeutPines", names(df)[str_detect(names(df), "NeutViewKragel")][-1]) # recode -1 to more safely exclude the kragel negEmotion pattern

# standardize pattern scores 
for(i in 1:length(pattern_Neut_Names)){
  df[, paste0(pattern_Neut_Names[i], "_scalePool")] <- scaleWithinRuns(pattern_Neut_Names[i], design = "between", group_var = "group")
}


#########################
# conduct frequentist mega-analysis

# prepare results df for mega.analysis
resultsMEGAbetweenNeut <- as.data.frame(
  matrix(
    nrow = length(pattern_Neut_Names),
    ncol = 4))
names(resultsMEGAbetweenNeut) <- c("pattern", "estimate", "lb", "ub")


# conduct mega analysis on all patterns
for(i in 1:length(pattern_Neut_Names)){
  
  dv <- paste0(pattern_Neut_Names[i], "_scalePool")
  f <- as.formula(paste(dv, "0 + group_weightEff + (1|study:ID)", sep = " ~ "))
  
  MEGAbetweenModel <- lmer(
    data = df,
    formula = f
  )
  
  cat("__________________________________________________________________________________\n\n",
      "PATTERN: ", pattern_Neut_Names[i], "\n\n\n")
  print(summary(MEGAbetweenModel))
  print(round(confint(MEGAbetweenModel), 2))
  
  est <- round(summary(MEGAbetweenModel)$coefficients[1], 2)
  lb <- round(confint(MEGAbetweenModel)[3, 1], 2)
  ub <- round(confint(MEGAbetweenModel)[3, 2], 2)
  
  resultsMEGAbetweenNeut[i, ] <- list(pattern_Neut_Names[i], est, lb, ub)
  
}

# show  results
resultsMEGAbetweenNeut


# plot results
setwd(paste0(wd, "/Figures"))

dfplotMEGAbetweenNeut <- resultsMEGAbetweenNeut
dfplotMEGAbetweenNeut$emotion_label <- factor(c("Negative Affect",  "Amused", "Angry", "Content", "Fearful", "Neutral", "Sad", "Surprised"), 
                                              levels = c("Negative Affect", "Fearful", "Angry", "Sad", "Amused", "Content", "Surprised", "Neutral"))

ggplot(data = dfplotMEGAbetweenNeut, aes(y = estimate, x = emotion_label)) +
  geom_point(size = 3) + geom_errorbar(aes(ymin = lb, ymax = ub)) +
  
  ylab(expression(atop("Pattern Expression (neutral condition)", paste("[Cohen's ", italic("d"), ", Emotion Dysregulation - Healthy Controls]")))) + xlab(NULL) +
  ylim(-1, 1) +
  
  theme_classic() +
  theme(axis.line.x=element_blank(), axis.ticks = element_blank(), axis.title.x=element_blank()) + 
  geom_segment(aes(x=0,y=0,xend=9,yend=0)) 

ggsave("PlotMEGAbetweenNeut.jpg", device = "jpeg")
ggsave("PlotMEGAbetweenNeut.pdf", device = "pdf")


#################################################
#  between mega-analysis on negative  baseline


#########################
# prepare variables

# get names of variables for pattern expression in negative  condition
pattern_Neg_Names <- c("NegPines", names(df)[str_detect(names(df), "NegViewKragel")][-1])

# standardize pattern scores 
for(i in 1:length(pattern_Neg_Names)){
  df[, paste0(pattern_Neg_Names[i], "_scalePool")] <- scaleWithinRuns(pattern_Neg_Names[i], design = "between", group_var = "group")
}


#########################
# conduct frequentist mega-analysis

# prepare results df for mega.analysis
resultsMEGAbetweenNeg <- as.data.frame(
  matrix(
    nrow = length(pattern_Neg_Names),
    ncol = 4))
names(resultsMEGAbetweenNeg) <- c("pattern", "estimate", "lb", "ub")


# conduct mega analysis on all patterns
for(i in 1:length(pattern_Neg_Names)){
  
  dv <- paste0(pattern_Neg_Names[i], "_scalePool")
  f <- as.formula(paste(dv, "0 + group_weightEff + (1|study:ID)", sep = " ~ "))
  
  MEGAbetweenModel <- lmer(
    data = df,
    formula = f
  )
  
  cat("__________________________________________________________________________________\n\n",
      "PATTERN: ", pattern_Neg_Names[i], "\n\n\n")
  print(summary(MEGAbetweenModel))
  print(round(confint(MEGAbetweenModel), 2))
  
  est <- round(summary(MEGAbetweenModel)$coefficients[1], 2)
  lb <- round(confint(MEGAbetweenModel)[3, 1], 2)
  ub <- round(confint(MEGAbetweenModel)[3, 2], 2)
  
  resultsMEGAbetweenNeg[i, ] <- list(pattern_Neg_Names[i], est, lb, ub)
  
}

# show  results
resultsMEGAbetweenNeg

# compare CIs for negative-neutral with negative
mean(resultsMEGAbetweenNeg$ub - resultsMEGAbetweenNeg$lb)/mean(resultsMEGAbetween$ub - resultsMEGAbetween$lb)


# plot results
setwd(paste0(wd, "/Figures"))

dfplotMEGAbetweenNeg <- resultsMEGAbetweenNeg
dfplotMEGAbetweenNeg$emotion_label <- factor(c("Negative Affect",  "Amused", "Angry", "Content", "Fearful", "Neutral", "Sad", "Surprised"), 
                                             levels = c("Negative Affect", "Fearful", "Angry", "Sad", "Amused", "Content", "Surprised", "Neutral"))

ggplot(data = dfplotMEGAbetweenNeg, aes(y = estimate, x = emotion_label)) +
  geom_point(size = 3) + geom_errorbar(aes(ymin = lb, ymax = ub)) +
  
  ylab(expression(atop("Pattern Expression (negative condition)", paste("[Cohen's ", italic("d"), ", Emotion Dysregulation - Healthy Controls]")))) + xlab(NULL) +
  ylim(-1, 1) +
  
  theme_classic() +
  theme(axis.line.x=element_blank(), axis.ticks = element_blank(), axis.title.x=element_blank()) + 
  geom_segment(aes(x=0,y=0,xend=9,yend=0)) 

ggsave("PlotMEGAbetweenNeg.jpg", device = "jpeg")
ggsave("PlotMEGAbetweenNeg.pdf", device = "pdf")

