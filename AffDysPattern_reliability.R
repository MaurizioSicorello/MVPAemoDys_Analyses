wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(wd, "/Data/Preprocessed_Data"))

df <- read.csv("EmoPattern_neuralData.csv")



#################################################
# load packages and functions

library("ltm")
library("plyr")
library("tidyr")
library("corrplot")
library("performance")
library("lavaan")


#################################################
# prepare variables

# calculate contrasts "negative vs neutral"
df$Pines_Diff_NegAffect <- df$NegPines - df$NeutPines
df$Kragel_Diff_Amused <- df$NegViewKragelAmused - df$NeutViewKragelAmused
df$Kragel_Diff_Angry <- df$NegViewKragelAngry - df$NeutViewKragelAngry
df$Kragel_Diff_Content <- df$NegViewKragelContent - df$NeutViewKragelContent
df$Kragel_Diff_Fearful <- df$NegViewKragelFearful - df$NeutViewKragelFearful
df$Kragel_Diff_Neutral <- df$NegViewKragelNeutral - df$NeutViewKragelNeutral
df$Kragel_Diff_Sad <- df$NegViewKragelSad - df$NeutViewKragelSad
df$Kragel_Diff_Surprised <- df$NegViewKragelSurprised - df$NeutViewKragelSurprised
df$Schulze_Diff_AmyHippL <- df$NegViewAmyHippL - df$NeutViewAmyHippL

# get variable names of difference scores for all patterns
pattern_Diff_Names <- names(df)[which(names(df) == "Pines_Diff_NegAffect"):which(names(df) == "Schulze_Diff_AmyHippL")]

# ID as factor
df$ID <- as.factor(df$ID)



#################################################
# reliability study 1

# make study 1 df
dfemoreg <- df[df$study == "emoreg", ]


study2rel <- as.data.frame(
  matrix(
    nrow = length(pattern_Diff_Names),
    ncol = 2,
    dimnames = list(NULL, c("pattern", "alpha"))
  )
)

for(i in 1:length(pattern_Diff_Names)){
  alphaModel <- cronbach.alpha(pivot_wider(dfemoreg, id_cols = "ID", names_from = "run", values_from = pattern_Diff_Names[i])[, -1], CI = TRUE)
  study2rel[i, ] <- list(pattern_Diff_Names[i], round(alphaModel$alpha, 2))
}


# results
study2rel

# reliability of patterns
study2relPattern <- study2rel[-nrow(study2rel), ]
round(mean(study2relPattern$alpha), 2)
range(study2relPattern$alpha)


# correlation matrix for activation between runs, averaged for all patterns

# create array of between-run correlation matrices, separately for patterns
pattern_Diff_Names_noAmy <- pattern_Diff_Names[-which(pattern_Diff_Names == "Schulze_Diff_AmyHippL")]
MAT <- vector(mode = "list", length = length(pattern_Diff_Names_noAmy))

for(i in 1:length(pattern_Diff_Names_noAmy)){
  MAT[[i]] <- cor(pivot_wider(dfemoreg, id_cols = "ID", names_from = "run", values_from = pattern_Diff_Names_noAmy[i])[, -1])
}
C <- do.call(cbind, MAT)
C <- array(C, dim=c(dim(MAT[[1]]), length(MAT)))

# calculate mean correlation matrix
apply(C, c(1,2), mean)



#################################################
# reliability study 2

# make study 2 df
dfreact <- df[df$study == "react" & df$group != "BPDremit", ]

study2rel <- as.data.frame(
  matrix(
    nrow = length(pattern_Diff_Names),
    ncol = 2,
    dimnames = list(NULL, c("pattern", "alpha"))
  )
)

for(i in 1:length(pattern_Diff_Names)){
  alphaModel <- cronbach.alpha(pivot_wider(dfreact, id_cols = "ID", names_from = "run", values_from = pattern_Diff_Names[i])[, -1], CI = TRUE)
  study2rel[i, ] <- list(pattern_Diff_Names[i], round(alphaModel$alpha, 2))
}


# results
study2rel

# reliability of patterns
study2relPattern <- study2rel[-nrow(study2rel), ]
round(mean(study2relPattern$alpha), 2)
range(study2relPattern$alpha)



# correlation matrix for activation between runs, averaged for all patterns

# create array of between-run correlation matrices, separately for patterns
pattern_Diff_Names_noAmy <- pattern_Diff_Names[-which(pattern_Diff_Names == "Schulze_Diff_AmyHippL")]
MAT <- vector(mode = "list", length = length(pattern_Diff_Names_noAmy))

for(i in 1:length(pattern_Diff_Names_noAmy)){
  MAT[[i]] <- cor(pivot_wider(dfreact, id_cols = "ID", names_from = "run", values_from = pattern_Diff_Names_noAmy[i])[, -1],
                  use = "pairwise.complete.obs")
}
C <- do.call(cbind, MAT)
C <- array(C, dim=c(dim(MAT[[1]]), length(MAT)))

# calculate mean correlation matrix
apply(C, c(1,2), mean)





