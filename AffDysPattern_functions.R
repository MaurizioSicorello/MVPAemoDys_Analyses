# pooled SD of x by groups "Clinical" and "Healthy"
sdPool <- function(x, group){
  d1 = x[group == "Clinical"]
  d2 = x[group == "Healthy"]
  
  n1 = length(d1)
  n2 = length(d2)
  
  v1 = var(d1)
  v2 = var(d2)
  
  sdPool = sqrt(((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2))
  return(sdPool)
}


# center on mean and standardize on pooled SD
scalePooled <- function(x, group){
  
  sdPooled <- sdPool(x, group)
  out <- (x-mean(x))/sdPooled
  
  return(out)
  
}


# scale within runs 
scaleWithinRuns <- function(x, data = df, design = "within", group_var = NULL){
  
  # within design
  if(design == "within"){
    
    res <- ddply(data, .(study, run), here(transform), out = get(x)/sd(get(x))) # no centering
    
    # between design  
  }else if(design == "between"){
    
    res <- ddply(data, .(study, run), here(transform), out = scalePooled(get(x), get(group_var))) # centering
    
    # warnings
  }else
  {warning("invalid design argument")
  }
  
  # order output variable according to input dataframe
  res_ordered <- res[
    order(
      match(
        paste(res[,"study"],res[,"run"],res[,"ID"]),
        paste(data[,"study"],data[,"run"],data[,"ID"])
      )
    )
    ,]
  
  # return scaled variable
  return(res_ordered$out)
  
}


# weighted effect coding (regression)
weightedEffectCoding <- function(group){
  #weighted codes
  code1 <- 1
  code2 <- -(sum(group == "Clinical")/sum(group == "Healthy"))
  #scale group difference to 1
  codeDiff <- code1-code2
  code1_stan <- code1/codeDiff
  code2_stan <- code2/codeDiff
  #output codes
  outputCodes <- ifelse(group == "Clinical", code1_stan, code2_stan)
  outputCodes
}


# FUNCTION to prepare data for three-level meta-analysis (cohens d and SE)
# formulas for independent cohens d's sampling variance of d based on Borenstein 2009.
# for dependent cohen's d's sampling variance: http://wvbauer.com/lib/exe/fetch.php/talks:2019_viechtbauer_lsp_ma_longitudinal.pdf
prepMeta <- function(x, data = df, betweenDesign = T){
  
  # empty df
  dfMAData <- as.data.frame(
    matrix(
      nrow=nrow(unique(df[, c("study", "run")])), 
      ncol = 4)
  )
  names(dfMAData) <- c("study", "run", "effectSize", "Variance")
  
  studyNames <- unique(df$study)
  counter <- 0
  
  
  # outer loop 
  #[ enter single study. get number of runs]
  for(i in 1:length(studyNames)){
    
    dfSub <- df[df$study == studyNames[i], ]
    studyName <- i
    nruns <- length(unique(dfSub$run))
    
    
    # inner loop
    # [calculate cohen's d and SE for each run]
    for(j in 1:nruns){
      
      counter <- counter + 1
      dfSubRun <- dfSub[dfSub$run == j, ]
      run <- j
      
      if(betweenDesign == T){
        
        cd <- cohens_d(dfSubRun[dfSubRun$group == "Clinical", x], 
                       dfSubRun[dfSubRun$group == "Healthy", x])$Cohens_d
        n1 <- nrow(dfSubRun[dfSubRun$group == "Clinical", ])
        n2 <- nrow(dfSubRun[dfSubRun$group == "Healthy", ])
        SEsquared <- (n1+n2)/(n1*n2)+cd^2/(2*(n1+n2))
      }
      else{
        
        cd <- cohens_d(dfSubRun[, x])$Cohens_d
        n <- length(dfSubRun[, x])
        SEsquared <- 1/n + cd^2/(2*n)
      }
      
      results <- c(studyName, run, cd, SEsquared)
      
      dfMAData[counter, ] <- results
      
      
    }
  }
  return(dfMAData)
}