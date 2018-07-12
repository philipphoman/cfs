#-----------------------------------------------------------------------
# Reversal learning CFS
#
# Functions
#
#-----------------------------------------------------------------------
# PH, 3/23/18
#-----------------------------------------------------------------------
rm(list = ls())
#library(plyr)
#library(dplyr)
#library(ggplot2)
#library(cowplot)
#library(QuantPsyc)
#library(car)

if (!require("pacman")) install.packages("pacman")
library("pacman")
pacman::p_load(
          "cowplot",
          "png",
          "astsa",
          "magick",
          "here",
          "png",
          "lme4",
          "lsmeans",
          "grid",
          "gridExtra",
          "flexmix",
          "DescTools",
          "grid",
          "plyr",
          "car",
          "QuantPsyc",
          "ggplot2",
          "tidyr",
          "dplyr"
        )




#-----------------------------------------------------------------------
readdat <- function(fn="../data/cfs.csv", df=NULL, minur=NULL,
                    subids=NULL, resptype="crur") {
  # Reads the csv source data and applies criteria as specified
  #
  # Args:
  #   fn      : csv data base file name
  #   usmin   : minimum mean US response for subjects to be included
  #   subids  : subject ids to be included
  #   resptype: response type (conditioned, unconditioned)
  #
  # Returns:
  #   df      : data frame of scanned csv values 
  if (is.null(df)) df <- read.csv(fn)

  # Add early and late phase
  df$phase[df$ctrial %in% c( 1: 4)]  <- "Early"
  df$phase[df$ctrial %in% c( 5: 8)]  <- "Late"
  df$phase[df$ctrial %in% c( 9:12)]  <- "Early"
  df$phase[df$ctrial %in% c(13:16)]  <- "Late"

  # Add safety marker
  df$safety[df$stim=="CSplus"&df$stage=="Acquisition"] <- "Unsafe"
  df$safety[df$stim=="CSplus"&df$stage=="Reversal"] <- "Safe"
  df$safety[df$stim=="CSminus"&df$stage=="Acquisition"] <- "Safe"
  df$safety[df$stim=="CSminus"&df$stage=="Reversal"] <- "Unsafe"
  df$safety <- factor(df$safety)

  dftmp <- df

  if (!is.null(minur)) {
	  oldsize     <- nrow(subset(dftmp, !(stim  == "CSplusUS"  |
                                        stim  == "CSminusUS" |
                                        group == ""          |
                                        stim  == "") & trial==2))
    #oldsize     <- nrow(subset(dftmp, trial==1))
    dftmpm      <- subset(dftmp, stim  == "CSplusUS"  |
                                 stim  == "CSminusUS" |
                                 group == ""          |
                                 stim  == "")
    dftmpmus    <- ddply(dftmpm, c("id"),
                         summarise, meanur=mean(scr, na.rm=TRUE))
    #p           <- print(dftmpmus$meanus)
    dftmp       <- merge(dftmp, dftmpmus, all=TRUE)
    dftmp       <- subset(dftmp, meanur > minur)
    newsize     <- nrow(subset(dftmp, trial==2))
    str         <- paste("Caution: removed", oldsize - newsize,
                         "subjects that did not pass the US criterion")
    print(str)
    #p           <- print(unique(dftmp$meanus))
  }

  if (!is.null(subids[1])) {
    tmp         <- print(paste("Caution: Restricting ids to the",
                                "list provided!"))
    dftmp       <- subset(dftmp, id %in% subids)
  }

  # restrict to non reinforced trials and complete cases; note that this
  # will remove the first trial as both orders start with a shock; thus,
  # to make analyses in wide format we will need to restrict the data
  # frame to the second trial (as this is is a non-reinforced trial in
  # both trial orders)
  switch(resptype,
         # restrict to CR
         cr   = {
           df  <- subset(dftmp, !(stim  == "CSplusUS"  |
                                  stim  == "CSminusUS" |
                                  group == ""          |
                                  stim  == ""))
         },
         # restrict to trials following a shock
         crds   = {
           df  <- subset(dftmp, (safety     == "Safe"  &
                                 precshock  == 1) |
                                 (safety     == "Unsafe" &
                                 precshock  == 1))
         },
         # restrict to CR and acquisition stage
         cra   = {
           df  <- subset(dftmp, !(stim  == "CSplusUS"  |
                                  stim  == "CSminusUS" |
                                  stage == "Reversal"  |
                                  group == ""          |
                                  stim  == ""))
         },
         # restrict to CR and acquisition early stage
         crae   = {
           df  <- subset(dftmp, !(stim  == "CSplusUS"  |
                                  stim  == "CSminusUS" |
                                  stage == "Reversal"  |
                                  group == ""          |
                                  stim  == "") &
                                  phase == "Early")
         },
         # restrict to CR and acquisition late stage
         cral   = {
           df  <- subset(dftmp, !(stim  == "CSplusUS"  |
                                  stim  == "CSminusUS" |
                                  stage == "Reversal"  |
                                  group == ""          |
                                  stim  == "") &
                                  phase == "Late")
         },
         # restrict to CR and acquisition late stage, CS+
         cralrlp   = {
           df  <- subset(dftmp, !(stim  == "CSplusUS"  |
                                  stim  == "CSminusUS" |
                                  group == ""          |
                                  stim  == "") &
                                  phase == "Late" &
                                  stim  == "CSplus")
         },
         # restrict to CR and reversal late stage, CS-
         cralrlm   = {
           df  <- subset(dftmp, !(stim  == "CSplusUS"  |
                                  stim  == "CSminusUS" |
                                  group == ""          |
                                  stim  == "") &
                                  phase == "Late" &
                                  stim  == "CSminus")
         },
         # restrict to CR and reversal stage
         crr   = {
           df  <- subset(dftmp, !(stim  == "CSplusUS"    |
                                  stim  == "CSminusUS"   |
                                  stage == "Acquisition" |
                                  group == ""            |
                                  stim  == ""))
         },
         # restrict to CR and reversal early stage
         crre   = {
           df  <- subset(dftmp, !(stim  == "CSplusUS"    |
                                  stim  == "CSminusUS"   |
                                  stage == "Acquisition" |
                                  group == ""            |
                                  stim  == "") &
                                  phase == "Early")
         },
         # restrict to CR and reversal late stage
         crrl   = {
           df  <- subset(dftmp, !(stim  == "CSplusUS"    |
                                  stim  == "CSminusUS"   |
                                  stage == "Acquisition" |
                                  group == ""            |
                                  stim  == "") &
                                  phase == "Late")
         },
         # restrict to UR
         ur   = {
           df  <- subset(dftmp, !(stim  == "CSplusUS"  |
                                  stim  == "CSminusUS" |
                                  group == ""          |
                                  stim  == ""))
         },
         # restrict to UR and CR
         crur = {
           df  <- subset(dftmp, !(group == "" |
                                  stim  == ""))
         },
         # restrict to UR acquisition
         urp  = {
           df  <- subset(dftmp,  (stim  == "CSplusUS") & 
                                !(stim  == "" |
                                  group == "" )) 
         },
         # restrict to UR reversal 
         urm  = {
           df  <- subset(dftmp,  (stim  == "CSminusUS") & 
                                !(stim  == "" |
                                  group == "" )) 
         },

         # restrict to subjects with non-exclude flag 
         nonexcluded = {
           df  <- subset(dftmp,  is.na(exclude))
         },

         # restrict to subjects with correct mri 
         withgoodmri = {
           df  <- subset(dftmp,  (hasrightnumofmrivols=="yes"))
         },

         # restrict to subjects with incorrect mri 
         withbadmri = {
           df  <- subset(dftmp,  (hasrightnumofmrivols=="no"))
         },

         # restrict to subjects with SCR data
         withscr = {
           df  <- subset(dftmp,  !is.na(scr))
         },

         # restrict to subjects with DCM data
         withdcm = {
           df  <- subset(dftmp,  !is.na(dcm2) & ntrials==32)
         },

         only32trials = {
           df  <- subset(dftmp, trial<=32)
         })
	
	df$group  <- factor(df$group)
	df$id     <- factor(df$id)
	df$stim   <- factor(df$stim)
  return(df)
}

createfulldf <- function(df=NULL, resptype="cr",
                         oc="scr", ocsur="scr") {
  # Creates full data frame by summarizing acq, rev, and revlearn (which
  # corresponds to the interaction of stage x stim) for each subject and
  # adding these values as to the data frame
  #
  # Args:
  #   df       : data frame
  #   resptype : type of SCR (conditioned=cr, unconditioned=ur, etc.)
  #   oc       : outcome measure
  #              (skin conductance peak2peak or dcm model based)
  #   ocsur    : outcome measure for unconditioned response
  #              (response to shocks)
  # Returns:
  #   dff      : merged data frame

  #df              <- subset(df, select=c(id, group, stage, stim, ctrial,
  #                                       dcm2, dcm3))

  dfus            <- calcmeanur(df=df, ocsur=ocsur)


  dfa             <- readdat(df=df, resptype="crds")
  dfa$outcome     <- dfa[, oc]
  out             <- calcstrengthint(df=dfa,
                                     fcts=c("id", "safety"), invf=1)
  dfads            <- out[[2]]
  dfads[, "meanrespds"]   <- dfads$m
  dfads$m          <- NULL
  dfads$stage      <- NULL
  dfads$phase      <- NULL
  
  
  dfa             <- readdat(df=df, resptype="cra")
  dfa$outcome     <- dfa[, oc]
  out             <- calcstrengthint(df=dfa,
                                     fcts=c("id", "stim"), invf=1)
  dfap            <- out[[2]]
  dfap[, "acqt"]   <- dfap$m
  dfap$m          <- NULL
  dfap$stage      <- NULL
  dfap$phase      <- NULL

  dfa             <- readdat(df=df, resptype="crae")
  dfa$outcome     <- dfa[, oc]
  out             <- calcstrengthint(df=dfa,
                                     fcts=c("id", "stim"), invf=1)
  dfaep            <- out[[2]]
  dfaep[, "acqearly"]   <- dfaep$m
  dfaep$m          <- NULL
  dfaep$stage      <- NULL
  dfaep$phase      <- NULL

  dfa             <- readdat(df=df, resptype="cral")
  dfa$outcome     <- dfa[, oc]
  out             <- calcstrengthint(df=dfa,
                                     fcts=c("id", "stim"), invf=1)
  dfalp            <- out[[2]]
  dfalp[, "acqlate"]   <- dfalp$m
  dfalp$m          <- NULL
  dfalp$stage      <- NULL
  dfalp$phase      <- NULL


  dfr             <- readdat(df=df, resptype="crr")
  dfr$outcome     <- dfr[, oc]
  out             <- calcstrengthint(df=dfr,
                                     fcts=c("id", "stim"), invf=-1)
  dfrp            <- out[[2]]
  dfrp[, "revt"]   <- dfrp$m
  dfrp$m          <- NULL
  dfrp$stage      <- NULL
  dfrp$phase      <- NULL


  dfr             <- readdat(df=df, resptype="crre")
  dfr$outcome     <- dfr[, oc]
  out             <- calcstrengthint(df=dfr,
                                     fcts=c("id", "stim"), invf=-1)
  dfrep            <- out[[2]]
  dfrep[, "revearly"]   <- dfrep$m
  dfrep$m          <- NULL
  dfrep$stage      <- NULL
  dfrep$phase      <- NULL

  dfr             <- readdat(df=df, resptype="crrl")
  dfr$outcome     <- dfr[, oc]
  out             <- calcstrengthint(df=dfr,
                                     fcts=c("id", "stim"), invf=-1)
  dfrlp            <- out[[2]]
  dfrlp[, "revlate"]   <- dfrlp$m
  dfrlp$m          <- NULL
  dfrlp$stage      <- NULL
  dfrlp$phase      <- NULL

  # "Inhibition" according to R01
  # Mean CS+ in acqlate minus mean CS+ in revlate
  dfalrlp             <- readdat(df=df, resptype="cralrlp")
  dfalrlp$outcome     <- dfalrlp[, oc]
  out             <- calcstrengthint(df=dfalrlp,
                                     fcts=c("id", "stage"), invf=-1)
  dfalrlplp            <- out[[2]]
  dfalrlplp[, "inhib"]   <- dfalrlplp$m
  dfalrlplp$m          <- NULL
  dfalrlplp$stage      <- NULL
  dfalrlplp$phase      <- NULL
  
  # "Update" according to R01
  # Mean CS- in acqlate minus mean CS- in revlate
  dfalrlm             <- readdat(df=df, resptype="cralrlm")
  dfalrlm$outcome     <- dfalrlm[, oc]
  out             <- calcstrengthint(df=dfalrlm,
                                     fcts=c("id", "stage"), invf=1)
  dfalrlmlp            <- out[[2]]
  dfalrlmlp[, "update"]   <- dfalrlmlp$m
  dfalrlmlp$m          <- NULL
  dfalrlmlp$stage      <- NULL
  dfalrlmlp$phase      <- NULL


  dfrl            <- readdat(df=df, resptype="cr")
  dfrl$outcome    <- dfrl[, oc]
  out             <- calcstrengthint(df=dfrl,
                                     fcts=c("id", "stage", "stim"),
                                     invf=-1)
  dfrelp           <- out[[2]]
  dfrelp[, "revlearn"]  <- dfrelp$m
  dfrelp$m         <- NULL
  dfrelp$stage     <- NULL

  dff             <- Reduce(function(x, y) merge(x, y, all=TRUE),
                            list(dfap, dfaep, dfalp, dfrp, dfrep,
                                 dfrlp, dfrelp, dfalrlplp,
                                 dfalrlmlp, dfus, dfads))

  dff$acq      <- (dff$acqearly + dff$acqlate)/2
  dff$rev      <- (dff$revearly + dff$revlate)/2

  # Apergis-Schoute2017, PNAS
  # "Generalization score": 1 - (sign(rev) * rev)
  dff$genscore<- 1 - (sign(dff$rev) * abs(dff$rev))


  # Stage difference: (face a - face b) - (face b - face a)
  dff$stagediff <- dff$acq - dff$rev
  
  #print(dff$id)
  #print(df$id)
  #dfff              <- merge(df, dff)
  
  return(dff)
}


calcmeanur <- function(df, ocsur="scr") {
  # Calculate the mean unconditioned response per subject
  #
  #   Args:
  #     df   : data frame
  #     ocsur: outcome type for shock response
  #
  #   Returns:
  #     df: data frame with new covariate meanur

  switch(ocsur,
         scr = {
           dfm <- ddply(df, c("id"), summarise,
                        meanur=mean(scr, na.rm=TRUE))
         },
         scrsqrt = {
           dfm <- ddply(df, c("id"), summarise,
                        meanur=mean(scrsqrt, na.rm=TRUE))
         },
         dcm3 = {
           dfm <- ddply(df, c("id"), summarise,
                        meanur=mean(dcm3, na.rm=TRUE))
         },
         dcm4 = {
           dfm <- ddply(df, c("id"), summarise,
                        meanur=mean(dcm4, na.rm=TRUE))
         })
  return(dfm)
}


preparedataset <-  function(df=NULL) {
  # Prepare the dataset for analysis
  #
  # Args:
  #  df: data frame
  # Returns:
  #  dfp: processed data frame

  r                                      <- df
  r$group                                <- factor(r$group)
  r$manipulation                         <- factor(r$manipulation)
  r$spider                               <- factor(r$spider)
  r$stage                                <- factor(r$stage)
  r$id                                   <- factor(r$id)
  r$ordertype <- NA
  r$ordertype[r$order=="A"|r$order=="C"] <- "RevSpAFirst"
  r$ordertype[r$order=="B"|r$order=="D"] <- "RevSpBFirst"

  r$spidernum <- NA
  r$spidernum[r$spider=="Spider A"]      <- 0
  r$spidernum[r$spider=="Spider B"]      <- 1

  r$visrespnum <- NA
  r$visrespnum[r$visresp==1]             <- 0
  r$visrespnum[r$visresp==2]             <- 1

  # set binary marker for shocks
  r$shockapplied <- NA
  r$shockapplied[r$spider=="Spider A"&r$stage=="Acquisition"]  <- 1
  r$shockapplied[r$spider=="Spider B"&r$stage=="Acquisition"]  <- 0
  r$shockapplied[r$spider=="Spider A"&r$stage=="Reversal"]  <- 1
  r$shockapplied[r$spider=="Spider B"&r$stage=="Reversal"]  <- 0 
  #rshockapplied  <- factor(r$shockapplied)

  dfp <-  r
  return(dfp)
}

plotscr <- function(df, l=TRUE, lty=1, lwd=1, p=TRUE, pch=1, pcex=1.0,
                    eb=TRUE, lcol="black", pcol="black", stard=0.5) {
  # Plots points, lines, and error bars of a scr-dataframe
  #
  # Args:
  #   df   : should include: ctrial, m, (eb=error bar)
  #   l    : plot lines (default=TRUE)
  #   lwd  : line width (default=1)
  #   p    : plot points (default=TRUE)
  #   pch  : format of points (default=1)
  #   eb   : plot error bars (default=TRUE)
  #   lcol : line color (default=black) 
  #   pcol : point color (default=black)
  #
  # Returns:
      
  if (l == TRUE) {
    l       <- lines(df$ctrial, df$m, cex=1.0, lwd=lwd, col=lcol,
                     lty=lty)
  }

  if (p == TRUE) {
    plt     <- points(df$ctrial, df$m,
                     type="p", pch=pch,
                     cex=pcex, col=pcol) 
  }
        
  if (eb == TRUE) {
    mses    <- df$m - df$eb
    pses    <- df$m + df$eb
    a       <- arrows(df$ctrial, mses, df$ctrial, pses, length=0.00,
                      angle=90, code=3, lwd=lwd)
  }
  if (!is.null(df$pval)) {
    st <- starsfromp(df$pval)
    text(df$ctrial, pses+stard, st)
    print(st)
  }
}

starsfromp <- function(pval, c1="~", c2="*") {
  # Parses pvalues, returns asterisks

  as <- vector(mode="character", length=length(pval))
  for (i in 1:length(pval)) {
    p <- pval[i]
    if (p < 0.1) as[i] <- c1 
    if (p < 0.05) as[i] <- paste(c2, sep="")
    if (p < 0.01) as[i] <- paste(c2, c2, sep="")
    if (p < 0.001) as[i] <- paste(c2, c2, c2, sep="")
  }
  return(as)
}

transfy <- function(y, transf) {
  # Applies transformation to outcome measure
  #
  # Args:
  #   y:      outcome measure
  #   transf: type of transformation (log, sqrt)
  #
  # Returns:
  #   yt:     y transformed

  switch(transf,
         log  = {
           yt <- log(y + 1)
         },
         sqrt = {
           yt <- sqrt(y) 
         },
         {
           # default
           yt <- y 
         })

  return(yt)
}

calcstrengthint <- function (df, fcts, invf=-1) {
  # Calculates the strength of a 3-way interaction per subject
  #
  # Args
  #   df  :    data frame with input values
  #   fcts:    factors by which data frame should by summarized
  #   invf:    multiply individual interactions with this factor
  #
  # Returns:
  #   output: list of data frames

  if (length(fcts) == 3) {
    # summarize the stimulus diff per subject in three steps
    pdm         <- ddply(df,    fcts,
                         summarise, m1=mean(outcome, na.rm=TRUE))
    pdm.s       <- ddply(pdm,   fcts[1:2],
                         summarise, m2=diff(m1, na.rm=TRUE))
    pdm.ss      <- ddply(pdm.s, fcts[1],
                         summarise, m=diff(m2, na.rm=TRUE))
  } else {
    # summarize the stimulus diff per subject in two steps
    #print(subset(df, select=c(id, group, stage, stim, trial, dcm2)))
    pdm         <- ddply(df,    fcts,
                         summarise, m1=mean(outcome, na.rm=TRUE))
    pdm.ss      <- ddply(pdm,   fcts[1],
                         summarise, m=diff(m1, na.rm=TRUE))
  }
  pdm.ss$m    <- pdm.ss$m * invf  
  pdmf        <- data.frame(m=NA, se=NA)
  #pdmf$m      <- abs(mean(pdm.ss$m, na.rm=TRUE))
  #pdmf$m      <- mean(pdm.ss$m, na.rm=TRUE) * (-1)
  pdmf$m      <- mean(pdm.ss$m, na.rm=TRUE) 
  pdmf$sd     <- sd  (pdm.ss$m, na.rm=TRUE)
  pdmf$N      <- sum(!is.na(pdm.ss$m))
  pdmf$se     <- pdmf$sd/sqrt(pdmf$N)
  output      <- list(pdmf, pdm.ss)
  return (output)
}

plotreg <- function (df, pcorr=TRUE, psig=TRUE,
                     xax=TRUE, yax=TRUE, pch=19,
                     main=NA, xlab=NA, ylab=NA, lbdp=FALSE) {
  # Plots linear regression with individual data points, regression
  # line, and significance level (optional)
  #
  # Args:
  #   df   : data frame
  #  pcorr : currently not implemented
  #  psig  : currently not implemented
  #  xax   : print x axis
  #  yax   : print y axis
  #  pch   : type of point
  #  main  : title for plot
  #  xlab  : x label
  #  ylab  : y label
  #  lbdp  : label individual data points 
  #
  # Returns:
  #rs.m <- subset(rs.m,trial==1&!is.na(revmean)&!is.na(rs.m[,covs[jj]]))
  plot(df$cov, df$m, xaxt="n", yaxt="n", xlab=xlab, ylab=ylab, pch=pch,
       cex=0.7, main=main, cex.main=0.7)
  if (xax  == TRUE) axis(1, cex.axis=0.7)
  if (!is.na(xlab)) mtext(xlab, 1, 2)
  if (yax  == TRUE) axis(2, cex.axis=0.7, las=2)
  if (!is.na(ylab)) mtext(ylab, 2, 3)
  if (lbdp == TRUE) text(df$cov, df$m, labels=df$id, cex=0.5)

  lm                  <- lm(df$m ~ df$cov)
  a                   <- abline(lm, lwd=2)
  rsq                 <- round(summary(lm)$r.squared, 3)
  rr                  <- round(sqrt(rsq), 2)
  cr                  <- cor.test(df$cov, df$m)
  pv                  <- round(cr$p.value, 3)
  ast                 <- ""
  if (pv < 0.1) ast   <- ""
  if (pv < 0.05) ast  <- "*"
  if (pv < 0.01) ast  <- "**"
  if (pv < 0.001) ast <- "***"
  rr                  <- round(cr$estimate, 2)
  legend("topright",paste("r=", rr, ast, sep = ""), bty = "n")
}
          
calcphtd <- function (phparams, lambda, stim) {
  # Prepares associability and prediction error curves
  # NB: Not yet implemented, needs to be checked
  #
  # Update rule:
  #       V = V + kappa * alpha * (lambda - V)
  #   alpha = eta * abs(lambda - V) + (1 - eta) * alpha 
  #
  # Args:
  #   phparams: vector of alpha_0, V_0, eta, kappa
  #   lambda  : vector of shocks
  #       stim: vector of stimuli
  #
  # Returns:
  #   output:   list of data frames
  initalpha <- phparams[1]
  initV     <- phparams[2]
  eta       <- phparams[3]
  kappa     <- phparams[4]
  V[1]      <- initV
  alpha[1]  <- initalpha
  delta[1]  <- NA
  for (i in 1:(length(stim)-1)) {
    delta[i]   <- lambda[i] - V[i]
    V[i+1]     <- V[i] + kappa * alpha[i] * delta[i]
    alpha[i+1] <- eta * abs(delta[i]) + (1 - eta) * alpha[i] 
  }
  df.a      <- data.frame(ctrial=c(1:length(alpha), m=alpha))
  df.p      <- data.frame(ctrial=c(1:length(delta), m=delta))
  output    <- list(df.a, df.p)
  return(output)
}

calcmergecov <- function(df1=NULL, df2=NULL, fcts=NULL,
                         name="mycov", invf=-1) {
  # Calculate covariate and merge with data frame
  #
  #  Args:
  #  Returns:
  #    dfm: merged data frame with new covariate mycov

  out         <- calcstrengthint(df2, fcts=fcts, invf=invf)
  dfi.ss      <- out[[2]]
  # Merge with original data frame
  dfm         <- merge(df1, dfi.ss, all=TRUE) 
  dfm[, name] <- dfm$m
  dfm$m       <- NULL
  return(dfm)
}

calcanova <- function(df=df, y=NULL) {
  # Calculate ANOVA on data frame given factors
  #
  #   Args:
  #     df  : data frame
  #     y   : outcome
  #   Returns:
  #     atab: ANOVA table

  atab <- vector(mode="list", length=length(y))
  dfa  <- subset(df, trial==2)
  for (i in 1:length(y)) {
    #print(y)
    dfa$y <- dfa[, y[i]]
    a <- aov(y ~ group, data=dfa)
    print(summary(a))
    #atab[i] <- a
  }
  #return(atab)
}

calclm <- function(df=df, resptype="cr", y=NULL,
                   grp=NULL, fcts=NULL, subids=NULL,
                   grpterm="* group", plot="n", printvif=NULL) {
  # Calculate linear models on y given factors
  #
  #   Args:
  #     df: data frame
  #      y: outcome
  #    grp: vector indicating if the group interaction is to
  #         be calculated
  #   fcts: vector of factors
  #  Returns:
  #      l: linear model output
  d   <- readdat(df=df, resptype=resptype, subids=subids)
  
  out <- list()
  ct  <- 0
  dfa <- subset(d, trial==2)
  for (i in 1:length(y)) {
    for (j in 1:length(fcts)) {
      # skip if y and factor are equal
      if (y[i] == fcts[j]) next
      ct        <- ct+1
      #dfa$y     <- dfa[, y[i]]
      #dfa$fct   <- dfa[, fcts[j]]
      f         <- paste(y[i], "~", fcts[j], ifelse(fcts[j] == "group",
                                                "", grpterm))
      cat("\n\n")
      print(paste("Model:", f))
      cat("\n\n")
      l         <- lm(f, data=dfa)
      print(summary(l))
      print(anova(l))
      #print(lm.beta(l))
      if(!is.null(printvif)) {
        print(vif(l))
      }
      if(plot=="y") {
        dff <- data.frame(cov=dfa[, fcts[j]], m=dfa[, y[i]])
        plotreg(dff)
      }
      out[[ct]] <- l
    }
  }
  return(out)
}
  

plotidresp <- function (fn=NULL, df=NULL, oc="scr", resptype="ur",
                        transf=NULL, subids=NULL,
                        scy="indiv", ti=NULL, ylim=c(0, 0.1),
                        xlim=c(1, 16)) {
  # Plots US responses for each subject
  #
  #  Args:
  #     fn       : file name of data base
  #     df       : data frame to work on 
  #     oc       : outcome type (SCR, DCM etc.)
  #     resptype : response type (conditioned/unconditioned response)
  #     transf   : transformation to apply on y
  #     subids   : vector of subject ids to restrict analysis
  #       scy    : scale y axis (indiv, group)
  #
  #  Returns:

  #r         <- readdat(fn=fn, df=df, resptype=resptype, subids=subids)
  r         <- readdat(fn=fn, df=df, resptype=resptype)
  # filter again
  filter    <- ifelse(grepl("scr", oc), "withscr", "withdcm")
  r         <- readdat(fn=fn, df=r,  resptype=filter)
  r$outcome <- transfy(r[, oc], transf) 

  #---------------------------------------------------------------------
  # Loop over ids and plot individual time courses
  #---------------------------------------------------------------------
  pr       <- par(mar=c(0.25, 0.25, 0, 0), oma=c(1, 1, 3.5, 1))
  l        <- layout(matrix(c(1:100), byrow=TRUE, ncol=10))
  #mxy      <- 1     
  for (s in 1:nlevels(r$id)) {
    bgc     <- "grey80"
    ri      <- subset(r, id==levels(r$id)[s])
    ri$m    <- ri$outcome
    ri      <- ri[order(ri$id, ri$trial),]

    # How to scale the y axis
    switch(scy,
           indiv = {
             if (is.na(ri$outcome[1])) {
               mxy = 0.1
               ylim     <- c(0, mxy)
             } else { 
               mxy <- max(ri$outcome, na.rm=TRUE) * 1.2
               ylim     <- c(0, mxy)
             }
           },
           group = {
             mxy      <- max(r$outcome, na.rm=TRUE) * 1.2
             ylim     <- c(0, mxy)
           },
           manual = {
             mxy      <- ylim[2]
           })
            
    tmp     <- plot(x=ri$ctrial, y=ri$outcome, xlab="", ylab="",
                    ylim=ylim, xlim=xlim,
                    col="white", axes=FALSE)

    # Flag the plot if id is excluded
    if (!ri$id[1] %in% subids) {
      print(paste("Adding exclusion flag to plot for", ri$id[1]))
      #mtext("EXCL", 3, -0.60, cex=0.3, adj=1, col="red", font=2)
      bgc="grey60"
    }

    rc      <- rect(par("usr")[1], par("usr")[3], par("usr")[2],
                    par("usr")[4], col=bgc, border=0)

    switch(resptype,
           cr   = {
             tmp <- plotscr(df=subset(ri, stim=="CSplus"),
                            p=FALSE, eb=FALSE, lcol="blue")
             tmp <- plotscr(df=subset(ri, stim=="CSminus"),
                            p=FALSE, eb=FALSE, lcol="red")
             lgl <- c("CS+", "CS-")
             lgc <- c("blue", "red")
           },
           ur   = {
             tmp <- plotscr(df=subset(ri, stim=="CSplus"),
                            p=FALSE, eb=FALSE, lcol="darkblue")
             tmp <- plotscr(df=subset(ri, stim=="CSminus"),
                            p=FALSE, eb=FALSE, lcol="darkred")
             lgl <- c("CSplusUS", "CSminusUS")
             lgc <- c("darkblue", "darkred")
           },
           sr   = {
             ri$ctrial <- ri$trial
             ri$m <- ri$spidernum
             tmp <- plotscr(df=ri, p=FALSE, eb=FALSE, lcol="red")
             ri$m <- ri$visrespnum
             tmp <- plotscr(df=ri, p=FALSE, eb=FALSE, lcol="blue")
             lgl <- c("Spider", "Response")
             lgc <- c("red", "blue")

           },
           shr   = {
             ri$ctrial <- ri$trial
             ri$m <- ri$shockapplied
             tmp <- plotscr(df=ri, p=FALSE, eb=FALSE, lcol="red")
             ri$m <- ri$visrespnum
             tmp <- plotscr(df=ri, p=FALSE, eb=FALSE, lcol="blue")
             lgl <- c("Shock", "Response")
             lgc <- c("red", "blue")
           })
    #tmp     <- plotscr(df=ri, p=FALSE, eb=FALSE, lcol="darkred")
    lb      <- paste(ri$id[1], " (", ri$group[1], ")", sep="")
    tmp     <- mtext(lb, 3, -0.60, cex=0.3, adj=0)

    # add also a tag that this is an excluded subject
    if (!ri$id[1] %in% subids) {
    #  print(paste("Adding exclusion flag to plot for", ri$id[1]))
      mtext("EXCL", 3, -0.60, cex=0.3, adj=1, col="red", font=2)
    }

    if (s == 1) {
      legend(x=5*max(ri$ctrial),
             y=mxy*1.6, lty=1, cex=0.75, lwd=1.5, col=lgc,
             legend=lgl, ncol=1, xpd=NA, bty="n")
      if (!is.null(ti)) title(ti, outer=TRUE, cex.main=0.8, adj=1)

    }

  }
}
  
plotgroupresp <- function (fn=NULL, df=NULL, oc="scr", resptype="cr",
                           transf=NULL, subids=NULL,
                           ylab=ylab, ti=NULL, plotind=TRUE) {
  # Plots group mean responses with error bars and strength of learning 
  #
  #  Args:
  #     fn       : file name of data base
  #     df       : data frame to work on 
  #     oc       : outcome type (SCR, DCM etc.)
  #     resptype : response type (conditioned/unconditioned response)
  #     transf   : transformation to apply on y
  #     subids   : vector of subject ids to restrict analysis
  #     ylab     : label for y axis
  #
  #  Returns:

  r         <- readdat(fn=fn, df=df, resptype=resptype, subids=subids)
  r$outcome <- transfy(r[, oc], transf) 

  #---------------------------------------------------------------------
  # Loop over groups and plot group means
  #---------------------------------------------------------------------
  pr       <- par(mar=c(2, 1, 2, 1), oma=c(22, 6, 6, 6))
  l        <- layout(matrix(c(1, 1, 1, 1, 1, 2,
                              3, 3, 3, 3, 3, 4),
                            byrow=TRUE, ncol=6))
  # layout.show(l)
	for (i in 1:nlevels(r$group)) {
    # create means and se
    rm          <- ddply(r, c("group", "stage", "ctrial", "stim"),
                         summarise, n=sum(!is.na(outcome)),
                         m=mean(outcome, na.rm=TRUE),
                         sd=sd(outcome, na.rm=TRUE), eb=sd/sqrt(n))
    mxy         <- max(rm$m)   * 1.2
    mxx         <- max(rm$ctrial) * 1.2
    rmg         <- subset(rm, group==levels(group)[i])
    rmm         <- subset(rmg, stim=='CSminus') 
    rmp         <- subset(rmg, stim=='CSplus') 
    tmp         <- plot(x=rmm$ctrial, y=rmm$m, xlab="", ylab="",
                        xlim=c(1, 16), ylim=c(0, mxy*1.1),
                        xaxt="n", bty="n")
    tmp         <- plotscr(df=rmm, pch=1)
    tmp         <- plotscr(df=rmp, pch=19)
    tmp         <- axis(1, at=c(1:16))

    if (i == 1) legend("top", inset=-0.25, lty=1, cex=1, lwd=1,
                       pch=c(19, 1), legend=c("CS+","CS-"),
                       ncol=1, xpd=NA, bty="n")
    if (i == 2) mtext(ylab, 2, 3, font=2, cex=1.5, outer=TRUE)
    if (i == 2) mtext(      "Trial", 1, 3, font=2, cex=1)

    # ------------------------------------------------------------------
    # Also plot strength of reversal learning with 95% CIs
    # ------------------------------------------------------------------
    #out         <- calcstrengthint(subset(r, group==levels(group)[i]),
    #                                     fcts=c("id", "stage", "stim"))
    #fi         <- out[[1]]


    if (!is.null(plotind)) {
      dfii        <- ddply(subset(r, trial==2), c("group"),
                           .fun = function(x) { 
                             c(m=mean(x[, 'revlearn'], na.rm=TRUE),
                               sd=sd(x[, 'revlearn'], na.rm=TRUE),
                               n=sum(!is.na(x[, 'revlearn']))
                               )})
      dfi         <- subset(dfii, group==levels(r$group)[i])

      ncp         <- abs(qt(0.05/2, dfi$n-1))
      dfi$se      <- dfi$sd/sqrt(dfi$n)
      dfi$ci      <- dfi$se * ncp
      dfi$eb      <- dfi$ci
      dfi$ctrial  <- 0.5
      tmp         <- plot(x=dfi$ctrial, y=dfi$m, xlab="", ylab="",
                          xlim=c(0, 1), ylim=c(0, mxy*1.1),
                          bty="n", axes=FALSE)
      tmp         <- plotscr(df=dfi, l=FALSE, pch=19, lwd=2)

      # label data point
      tmp         <- text(0.8, dfi$m, labels=round(dfi$m, 2), cex=0.5)

      lb          <- paste(levels(rm$group)[i], " (N=", rmg$n[1], ")",
                           sep="")
      tmp         <- abline(h=0, lty=3)
      tmp         <- text(1.75, mxy/2, labels=lb,
                          font=2, cex=1.5, srt = -90, xpd=NA)
      if (i == 3) mtext("Reversal\nLearning\nIndex", 1, 4, font=2,
                        cex=0.75)
    } else {
      plot.new()
      #lb          <- paste(levels(rm$group)[i], " (N=", rmg$n[1], ")",
      #                     sep="")
      #tmp         <- text(1.75, mxy/2, labels=lb,
      #                    font=2, cex=1.5, srt = -90, xpd=NA)
    }
	}
  if (!is.null(ti)) title(ti, outer=TRUE, cex.main=0.8)
}


plotgroupreg <- function (fn=NULL, df=NULL, oc="scr", resptype="cr",
                          transf="", subids=NULL, pch=19,
                          lbdp=FALSE, ylab=ylab, covs=NULL,
                          invf=-1, ti=NULL) {
  # Plots group linear regressions
  #
  #  Args:
  #     fn       : file name of data base
  #     df       : data frame to work on 
  #     oc       : outcome type (SCR, DCM etc.)
  #     resptype : response type (conditioned/unconditioned response)
  #     transf   : transformation to apply on y
  #     subids   : vector of subject ids to restrict analysis
  #     pch      : point format
  #     lbdp     : label data points
  #     covs     : vector of covariates
  #
  #  Returns:
  r         <- readdat(fn=fn, df=df, resptype=resptype, subids=subids)
  #r <- df
  #r         <- readdat(fn=fn, df=df, resptype=resptype)
  r$outcome <- transfy(r[, oc], transf) 

  #---------------------------------------------------------------------
  # Loop over groups and plot regressions
  #---------------------------------------------------------------------

  pr   <- par(mar=c(3, 2, 2, 0), oma=c(13, 4, 9, 0))
  l    <- layout(matrix(c(c( 1: 8),
                          c( 9:16),
                          c(17:24)), byrow=TRUE, ncol=8),
                 widths=c(rep(1, 7), 0.5), heights=c(1, 1, 1))

  #pr   <- par(mar=c(3, 2, 2, 0), oma=c(5, 4, 3, 0))
  #l    <- layout(matrix(c(c( 1: 8),
  #                        c( 9:16),
  #                        c(17:24),
  #                        c(25:32),
  #                        c(33:40)), byrow=TRUE, ncol=8),
  #               widths=c(rep(1, 7), 0.5), heights=c(1, 1, 1))

	for (i in 1:nlevels(r$group)) {
    yax    <- TRUE
    for (c in 1:length(covs)) {
      #if (covs[c] == "") {
      #  plot.new()
      #  next
      #}
      dftmp     <- subset(r, group==levels(group)[i])
      dftmp$cov <- dftmp[, covs[c]]
      dftmp     <- unique(subset(dftmp, select=c(m, cov, id)))
      #print(nrow(dftmp))
      if (c %% 7 != 1)  yax <- FALSE

      #print("Ok")
      tmp    <- plotreg(dftmp,
                        xax=TRUE, yax=yax, xlab=NA, lbdp=lbdp, pch=pch,
                        main=toupper(covs[c]))

      #tmp    <- plotreg(subset(dftmp, trial==2&!is.na(m)&!is.na(cov)),
      #                  xax=TRUE, yax=yax, xlab=NA, lbdp=lbdp, pch=pch,
      #                  main=toupper(levels(covs)[c]))

      #lb     <- "Reversal Learning"
      lb     <- ylab 
      if (c %% 7 == 1 && i == 1) mtext(lb, 2, 2, cex=1.5, font=2,
                                       outer=TRUE)

      if (c %% 7 == 0) {
        plot.new()
        lb          <- levels(r$group)[i]
        t           <- text(-0.1, 0.5, labels=lb,
                            font=2, cex=1.5, srt = -90, xpd=NA)
      } 
    }
  } 
  if (!is.null(ti)) title(ti, outer=TRUE, cex.main=0.8)
}

## plotgroupreg <- function (fn=NULL, df=NULL, oc="scr", resptype="cr",
##                           fcts, transf=NULL, subids=NULL, pch=19,
##                           lbdp=FALSE, ylab=ylab, covs=NULL, invf=-1) {
##   # Plots group linear regressions
##   #
##   # Note: this function should be revised and just take y as argument
##   # instead of calculating it on the fly
##   #
##   #  Args:
##   #     fn       : file name of data base
##   #     df       : data frame to work on 
##   #     oc       : outcome type (SCR, DCM etc.)
##   #     resptype : response type (conditioned/unconditioned response)
##   #     transf   : transformation to apply on y
##   #     subids   : vector of subject ids to restrict analysis
##   #     pch      : point format
##   #     lbdp     : label data points
##   #     covs     : vector of covariates
##   #     ylab     : y label
##   #     invf     : inversion factor (to switch sign of a difference) 
##   #
##   #  Returns:

##   r         <- readdat(fn=fn, df=df, resptype=resptype, subids=subids)
##   r$outcome <- transfy(r[, oc], transf) 

##   #---------------------------------------------------------------------
##   # Loop over groups (hcc, vcc, vptsd) and plot regressions
##   # Covariates are: CAPS, PCL
##   #---------------------------------------------------------------------
##   pr   <- par(mar=c(3, 2, 2, 0), oma=c(13, 4, 9, 0))
##   l    <- layout(matrix(c(c( 1: 8),
##                           c( 9:16),
##                           c(17:24)), byrow=TRUE, ncol=8),
##                  widths=c(rep(1, 7), 0.5), heights=c(1, 1, 1))
##   #layout.show(l)
##   if (is.null(covs)) {
##     covs <- c("pcl",
##               "factor1pcl",
##               "factor2pcl",
##               "factor3pcl",
##               "factor4pcl",
##               "factor5pcl",
##               "caps")
##   }
## 	for (i in 1:nlevels(r$group)) {
##     out    <- calcstrengthint(subset(r, group==levels(group)[i]),
##                               #fcts=c("id", "stage", "stim"))
##                               fcts=fcts, invf=invf)
##     dfi.ss <- out[[2]]
##     yax    <- TRUE
##     for (c in 1:length(covs)) {
##       dftmp     <- merge(dfi.ss, subset(r, group==levels(group)[i],
##                                         all=TRUE))
##       dftmp$cov <- dftmp[, covs[c]]
##       #dftmp$m   <- abs(dftmp$m)
##       #dftmp$m   <- dftmp$m * (-1)
##       dftmp     <- unique(subset(dftmp, select=c(m, cov, id)))
##       #print(paste("There are", nrow(dftmp), "datapoints"))
##       #print(dftmp)
##       #dftmp     <- subset(dftmp, trial==2, select=c(m, cov, id))
##       if (c %% 7 != 1)  yax <- FALSE

##       tmp    <- plotreg(dftmp,
##                         xax=TRUE, yax=yax, xlab=NA, lbdp=lbdp, pch=pch,
##                         main=toupper(covs[c]))

##       #tmp    <- plotreg(subset(dftmp, trial==2&!is.na(m)&!is.na(cov)),
##       #                  xax=TRUE, yax=yax, xlab=NA, lbdp=lbdp, pch=pch,
##       #                  main=toupper(levels(covs)[c]))

##       #lb     <- "Reversal Learning"
##       lb     <- ylab 
##       if (c %% 7 == 1 && i == 2) mtext(lb, 2, 2, cex=1.5, font=2,
##                                        outer=TRUE)

##       if (c %% 7 == 0) {
##         plot.new()
##         lb          <- levels(r$group)[i]
##         t           <- text(-0.1, 0.5, labels=lb,
##                             font=2, cex=1.5, srt = -90, xpd=NA)
##       } 
##     }
##   } 
## }

plotgroupstages <- function (fn=NULL, df=NULL, oc="scr", resptype="cr",
                             transf=NULL, subids=NULL, eb="ci",
                             ylab=NULL, dosc=NULL, yax=NULL,
                             abl=NULL, stages=NULL, ti=NULL,
                             mxyif=0.5, ylims=NULL) {
  # Plots stim diffs per stage
  #
  #  Args:
  #     fn       : file name of data base
  #     df       : data frame to work on 
  #     oc       : outcome type (SCR, DCM etc.)
  #     resptype : response type (conditioned/unconditioned response)
  #     transf   : transformation to apply on y
  #     subids   : vector of subject ids to restrict analysis
  #     eb       : type of error bar
  #     ylab     : label for y axis
  #     dosc     : vector of scaling factors for y axis
  #      yax     : boolean vector for drawing y axis
  #     stages   : vector summarizing stim diff per stage
  #     ti       : title string
  #     mxyif    : mxy inflation factor 
  #
  #  Returns:

  df         <- readdat(fn=fn, df=df, resptype=resptype, subids=subids)
  df$outcome <- transfy(df[, oc], transf) 

  #---------------------------------------------------------------------
  # Loop over groups and plot stim diffs 
  #---------------------------------------------------------------------
  pr         <- par(mar=c(5, 1, 2, 2), oma=c(2, 6, 2, 6))
  l          <- layout(matrix(c(1,  2,  3,  4,
                                5,  6,  7,  8,
                                9, 10, 11, 12),
                            byrow=TRUE, ncol=4))

  if (is.null(stages)) {
    stages     <- c("acqearly", "acqlate", "revearly", "revlate", "acq",
                    "rev", "inhib", "update", "revlearn")
  }
  if (is.null(yax)) {
    yax        <- c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE,
                    TRUE)
  }
  abl        <- c(rep(TRUE, length(stages)))
  #dosc       <- c(2.4, NA, NA, NA, 1.4, NA, NA, NA, 1.4)
  #dfi        <- ddply(subset(df, trial==2), c("group"), summarise,
  #                     m=mean(acq, na.rm=TRUE),
  #                     sd=sd(acq, na.rm=TRUE),
  #                     n=sum(!is.na(acq)),
  #                     se=sd/sqrt(n))
  
  for (i in 1:length(stages)) {
    if (stages[i] == "") {
      plot.new()
      next
    }
    dfi        <- ddply(subset(df, trial==2), c("group"),
                        .fun = function(x) { 
                          c(m=mean(x[, stages[i]], na.rm=TRUE),
                            sd=sd(x[, stages[i]], na.rm=TRUE),
                            n=sum(!is.na(x[, stages[i]]))
                            )})
    dfi$ctrial <- c(1:nlevels(df$group))
    ncp        <- abs(qt(0.05/2, dfi$n-1))
    dfi$se     <- dfi$sd/sqrt(dfi$n)
    dfi$ci     <- dfi$se * ncp
    #dfi$eb     <- ifelse(eb=="ci", dfi$ci, dfi$se)
    if (eb == "ci") {
      dfi$eb <- dfi$ci
    } else {
      dfi$eb <- dfi$se
    }
    
      
    print(dfi)

    # take care of scaling
    # also consider that values might be negative
    if (!is.na(dosc[i])) {
      mxy <- max(abs(dfi$m)) * dosc[i]
      if(mean(dfi$m, na.rm=TRUE) < 0) { 
        ylim <- c(-mxy, mxy*mxyif)
      } else {
        ylim <- c(-mxy*mxyif, mxy)
      }
    } else {
      ylim <- c(-max(abs(dfi$m))*0.5, max(abs(dfi$m)))
    }

    if(!is.null(ylims[i, 1])) {
      ylim <- ylims[i, 1:2]
    }

    print(ylim)

    mxx        <- max(dfi$ctrial) * 1.1

    tmp        <- plot(x=dfi$ctrial, y=dfi$m, xlab="", ylab="",
                       xlim=c(0.8, mxx), ylim=ylim,
                       bty="n", xaxt="n", yaxt="n", main=stages[i])

    tmp        <- axis(1, at=c(1:nlevels(dfi$group)),
                       labels=levels(dfi$group),
                       las=2, font=2)

    tmp        <- plotscr(df=dfi, l=FALSE, pch=19, lwd=2)
    tmp        <- text(dfi$ctrial+0.4, dfi$m,
                       labels=round(dfi$m, 2), cex=0.5, xpd=NA)
    if (yax[i] == TRUE) axis(2)                           
    if (abl[i] == TRUE) abline(h=0, lty=3)
  }
  mtext(ylab, 2, 3, font=2, outer=TRUE, cex=1.5)
  if (!is.null(ti)) title(ti, outer=TRUE, cex.main=0.8)
}

plotgroupfp <- function (fn=NULL, df=NULL, oc="scr", resptype="cr",
                         transf=NULL, subids=NULL, ti=NULL) {
  # Plots group free parameters 
  #
  #  Args:
  #     fn       : file name of data base
  #     df       : data frame to work on 
  #     oc       : outcome type (SCR, DCM etc.)
  #     resptype : response type (conditioned/unconditioned response)
  #     transf   : transformation to apply on y
  #     subids   : vector of subject ids to restrict analysis
  #
  #  Returns:

  r         <- readdat(fn=fn, df=df, resptype=resptype, subids=subids)
  r$outcome <- transfy(r[, oc], transf) 

  #---------------------------------------------------------------------
  # Loop over groups and plot PH-TD means
  #---------------------------------------------------------------------
  pr       <- par(mar=c(7, 1, 2, 2), oma=c(6, 6, 6, 6))
  l        <- layout(matrix(c(1, 2, 3, 4,
                              5, 6, 7, 8),
                            byrow=TRUE, ncol=4))
  ph       <- c("alpha", "v", "eta", "kappa")

  # plot fully constrained und partly constrained estimates
  lbs      <- c("Fully Constrained", "Kappa unconstrained")
  suffix   <- c("ab","")
  for (i in 1:2) {
    phparams        <- c(paste(ph[1], oc, suffix[i], sep=""),
                         paste(ph[2], oc, suffix[i], sep=""),
                         paste(ph[3], oc, suffix[i], sep=""),
                         paste(ph[4], oc, suffix[i], sep=""))
                                     
    phparamsgr      <- c(expression(bold(alpha[0])),
                         expression("V"[0]),
                         expression(bold(eta)),
                         expression(bold(kappa)))

    # each of the four free parameters
    for (j in 1:length(phparams)) {
      
      #print(paste("Working on", phparams[j]))
      #for (j in 1:nlevels(r$group)) {
      php         <- r[, phparams[j]] 
      r$phparam   <- php


      # Summarize
      rm          <- ddply(subset(r, trial==2), c("group"),
                           summarise, n=sum(!is.na(phparam)),
                           m=mean(phparam, na.rm=TRUE),
                           sd=sd(phparam, na.rm=TRUE), eb=sd/sqrt(n))

      # take care of y axis scaling
      mxy         <- ifelse(i==2&j==4, max(rm$m * 1.5), 1)
      
      rm$ctrial   <- c(1:nlevels(r$group))
      #mxy         <- max(rm$m)   * 1.2
      mxx         <- max(rm$ctrial) * 1.1
      #ncp         <- abs(qt(0.05/2, dfi$N-1))
      #dfi$se      <- dfi$se * ncp
      tmp         <- plot(x=rm$ctrial, y=rm$m, xlab="", ylab="",
                          xlim=c(0.8, mxx), ylim=c(0, mxy),
                          bty="n", xaxt="n", yaxt="n")
      tt          <- title(phparamsgr[j], cex.main=1.8, font=2)
      tmp         <- axis(1, at=c(1:nlevels(rm$group)),
                          labels=levels(rm$group),
                          las=2, font=2)
      
      tmp         <- plotscr(df=rm, l=FALSE, pch=19, lwd=2)
      if (i > 0) axis(2)
      if (j == 4)  text(3.75, mxy/2, labels=lbs[i],
                        font=2, cex=1.5, srt = -90, xpd=NA)

    }
  }
  if (!is.null(ti)) title(ti, outer=TRUE, cex.main=0.8)
}

plothist <- function(fn=NULL, df=NULL, oc="scr",
                     resptype="cr", transf=NULL,
                     subids=NULL, covs=NULL, abl=NULL,
                     ti=NULL, xlab=NULL) {
  # Plots histograms for each covariate in a vector
  #
  #   Args:
  #     fn    : file name
  #     df    : data frame
  #     oc    : outcome measure
  #   resptype: type of response
  #    transf : transformations to be applied on outcome measure
  #    subids : vector of included subject ids
  #    covs   : vector of covariates
  #    abl    : vector of vertical ablines to be plotted
  #    ti     : title string
  #    xlab   : x label
  #   Returns:

  r         <- readdat(fn=fn, df=df, resptype=resptype, subids=subids)
  r$outcome <- transfy(r[, oc], transf) 

  #---------------------------------------------------------------------
  # Loop over groups and plot PH-TD means
  #---------------------------------------------------------------------
  pr       <- par(mar=c(7, 1, 2, 2), oma=c(6, 6, 6, 6))
  l        <- layout(matrix(c( 1,  2,  3,  4,  5,
                               6,  7,  8,  9, 10,
                              11, 12, 13, 14, 15),
                            byrow=TRUE, ncol=5))
  for (h in 1:nlevels(r$group)) {
    rs <- subset(r, group==levels(r$group)[h])
    for (i in 1:length(covs)) {
      hi         <- hist(rs[, covs[i]], plot=FALSE)
      hi$density <- hi$counts/sum(hi$counts) * 100
      tmp        <- plot(hi, col="gray52", xlab="", freq=FALSE,
                         main=paste(toupper(covs[i])), ylim=c(0, 100))
      if (!is.na(abl[i])) {
        print(paste("Cov is", covs[i], ", index is", i, ", abl is", abl[i]))
        abline(v=abl[i], col="red", lty=2, lwd=1.2)
      }
      if (i == 4) {
        plot(x=NULL, y=NULL, xlim=c(0, 1), ylim=c(0, 100), axes=FALSE,
             xlab="", ylab="", col="white")
        lb          <- levels(r$group)[h]
        tmp         <- text(-0.2, 50, labels=lb,
                            font=2, cex=1.5, srt = -90, xpd=NA)
      }
      
    }
  }
  tmp <- mtext("Percentage", 2, 3, cex=1.75, font=2, outer=TRUE, las=3)
  tmp <- mtext(xlab, 1, 0, cex=1.75, font=2, outer=TRUE)
  if (!is.null(ti)) title(ti, outer=TRUE, cex.main=1.8)
}

makemat <- function(li=NULL) {
  # Create a matrix of all coefficients and statistics from list of lm
  # output
  #
  # Args:
  #   li : list of lm objects
  # Returns:
  #  mat : matrix

  l   <- lapply(li, function(x) {
    r=rev(append(rownames(summary(x)$coef), format(formula(x))))
    #m=rbind(NA, summary(x)$coef)
    #rownames(m)=r
    })
    
  mat <- do.call("rbind", l)
  return(mat)
}

makematanova <- function(li=NULL) {
  # Create a matrix of all ANOVA statistics from list of lm
  # output
  #
  # Args:
  #   li : list of lm objects
  # Returns:
  #  mat : matrix

  l   <- lapply(li, function(x) anova(x))
  mat <- do.call("rbind", l)
  return(mat)
}

sumupas2017 <- function(dfas=NULL) {
  # Sum up the findings from the paper by Apergis-Shoute et al. 2017,
  # PNAS
  #
  # Note: the SCR average values are read out from the graphs.
  #

  # sample sizes
  nhcc                           <- 35
  nocd                           <- 43

  pdf("../lib/cfs_as2017.pdf")
  # graphics
  pr                             <- par(mar=c(14, 1, 2, 2), oma=c(2, 6,
                                                                  2, 6))
  l                              <- layout(matrix(c( 1, 1),  
                                                  byrow=TRUE, ncol=2))

  # summarize Fig. 1
  dfas <- data.frame(group =c(rep(c("HCC", "HCC", "OCD", "OCD"), 4)),
                     ctrial=c(1:16),
                     stage =c(rep("Acquisition", 8), rep("Reversal", 8)),
                     phase =c(rep(c(rep("Early", 4), rep("Late", 4)), 2)),
                     stim  =c(rep(c("CSplus", "CSminus"), 8)),
                     m     =c(1.46, 0.66,  1.12, 0.51,
                              1.45, 0.508, 1.02, 0.48,
                              1.32, 0.9,   0.87, 0.66,
                              1.5,  0.70,   0.9,  0.69),
                     se    =c(0.14, 0.08, 0.09, 0.09, 0.16, 0.08, 0.1,
                              0.08, 0.13, 0.13, 0.1, 0.08, 0.22, 0.13,
                              0.135, 0.1))
  dfas$eb <- dfas$se

  # plot and compare with Fig. 1
  plot(dfas$ctrial, dfas$m, axes=FALSE, xlab="", ylab="",
       ylim=c(0, 1.7), col="white")
  plotscr(df=subset(dfas, stim=="CSplus"&group=="HCC"&
                          ctrial %in% seq(1, 16, 4)),
                    pch=19, l=FALSE, pcol="black")
  plotscr(df=subset(dfas, stim=="CSminus"&group=="HCC"&
                          ctrial %in% seq(2, 16, 4)),
                    pch=1, l=FALSE, pcol="black")
  plotscr(df=subset(dfas, stim=="CSplus"&group=="OCD"&
                          ctrial %in% seq(3, 16, 4)),
                    pch=19, l=FALSE, pcol="black")
  plotscr(df=subset(dfas, stim=="CSminus"&group=="OCD"&
                          ctrial %in% seq(4, 16, 4)),
                    pch=1, l=FALSE, pcol="black")

  axis(1, at=c(1:16), labels=paste(dfas$group, dfas$stage,
                                   dfas$phase, dfas$stim), las=3)
  axis(2)
  mtext("Mean SCR with SEM", 2, 3, font=2)

  # summarize as differential SCR
  dfasm                          <- ddply(dfas, c("group", "stage",
                                                  "phase"), summarise,
                                          mm= (-1) * diff(m))

  # calculate noncentrality parameters per group
  ncphcc                         <- abs(qt(0.05/2, nhcc-1))
  ncpocd                         <- abs(qt(0.05/2, nocd-1))

  # sample sizes
  dfasm$n[dfasm$group=="HCC"]    <- nhcc
  dfasm$n[dfasm$group=="OCD"]    <- nocd

  # use reported t-values to calculate appropriate standard errors
  dfasm$t                        <- c(5.812, 6.630, 3.274, 3.836, 5.808,
                                      5.196, 1.562, 1.056)
  dfasm$se                       <- dfasm$mm/dfasm$t
  dfasm$pval                     <- round(2 * pt(dfasm$t, dfasm$n-1,
                                                 lower=FALSE), 4)
  dfasm$ci[dfasm$group=="HCC"]   <- dfasm$se[dfasm$group=="HCC"] * ncphcc
  dfasm$ci[dfasm$group=="OCD"]   <- dfasm$se[dfasm$group=="OCD"] * ncpocd
  dfasm$sse                      <- NULL
  dfasm$sd                       <- dfasm$se * sqrt(dfasm$n)
  dfasm$var                      <- dfasm$sd^2
  

  # calculate the possible between group differences and sum up 
  dfasn                          <- ddply(dfasm, c("stage", "phase"),
                                          summarise,
                                          m= (-1) * diff(mm),
                                          sse=sqrt(sum(se^2)))
  dfasn$ncp                      <- abs(qt(0.05/2, (nhcc+nocd-1)))
  dfasn$eb                       <- dfasn$sse * dfasn$ncp
  dfasn$ctrial                   <- c(1:4)

  # plot
  pr                             <- par(mar=c(4, 1, 10, 11), oma=c(2, 6,
                                                                  2, 6))
  l                              <- layout(matrix(c( 1),  
                                                  byrow=TRUE, ncol=1))
  # summarize Fig. 1
  plot(dfasn$ctrial, dfasn$m, axes=FALSE, xlab="", ylab="", bty="n",
       ylim=c(-0.2, 1.2))
  plotscr(df=dfasn, l=FALSE, pch=19, lwd=2)
  axis(1, at=c(1:4), labels=c("Acq Early",
                              "Acq Late", "Rev Early", "Rev Late"),
       font=2, las=3)
  axis(2)
  mtext("Mean differential SCR with 95% CI", 2, 3, font=2)
  mtext("Group differences", 3, 0, font=2)
  abline(h=0, lty=3, lwd=1.2)
  dev.off()
}

createthstr <- function(minur=-Inf, minacq=-Inf, minrev=-Inf,
                        maxperc=Inf, hasprob=NA) {
  # Create a string indicating the applied threshold for inclusion
  #
  #   Args:
  #     minur: minimal mean UR
  #     minacq: minimal mean acq
  #     minrev: minimal mean rev
  #
  #   Returns:
  #     str: string
  str <- NULL
  if (!minur == -Inf) {
    str <- paste("Inclusion criterion: minimal mean UR >", minur)
  } else if (!minacq == -Inf) {
    str <- paste("Inclusion criterion: minimal acq >", minacq)
  } else if (!minrev == -Inf) {
    str <- paste("Inclusion criterion: minimal rev >", minrev)
  }
  return(str)
}

regroup <- function(fn=NULL, df=NULL, oc="scr",
                    resptype="cr", transf=NULL,
                    subids=NULL) {
  # Regroup the participants according to perfomance in reversal index
  #
  # Args:
  # Returns:
  r            <- readdat(fn=fn, df=df, resptype=resptype, subids=subids)
  #r$outcome    <- transfy(r[, oc], transf) 
  rr           <- subset(r, trial==2, select=c(id, group, revlearn))

  #meanrl       <- mean(subset(r, trial==2, select=revlearn)$revlearn,
  #                     na.rm=TRUE)
  
  library('Hmisc')
  rr$newgroup   <- as.numeric(cut2(rr$revlearn, g=3))
  rr$newgroup[rr$newgroup==1] <- "WeakRI"
  rr$newgroup[rr$newgroup==2] <- "MedRI"
  rr$newgroup[rr$newgroup==3] <- "StrongRI"
  df <- merge(r, rr, all=TRUE)
  #df <- rr
  df$newgroup <- factor(df$newgroup,
                        levels=c("WeakRI", "MedRI", "StrongRI"))
  return(df)
}

mergegroups <- function(df=NULL, act="merge") {
  # Merges VCC and VPTSD
  #df$group <- factor(df$group, levels=c("Aware", "Unaware", "All")) 
  switch(act,
         merge={
           df$oldgroup <- df$group
           df$group <- "All"
         },
         unmerge={
           df$group <- df$oldgroup
         })
  df$group <- factor(df$group)
  return(df)
}


plotbar <- function(matmean, mateb, col, xlim, ylim,
                    spc=c(0, 0.2), axes=TRUE,
                    xpd=TRUE, abl=NULL, main=NULL,
                    yaxt=NULL, ylb=NULL, cxf=2,
                    cxn=1.2) {
  # Bar plots with error bars
  #
  #   Args:
  #
  #   Returns:

  # create means and se
  #mxy         <- max(df$m)*1.45
  #tabbedmeans <- tapply(df$m, factors, function(x) c(x=x))
  #print(tabbedmeans)
  #tabbedeb    <- tapply(df$eb, factors, function(x) c(x=x))
  barcenters  <- barplot(height=matmean,
                         beside=TRUE,
                         ylim=ylim,
                         col=col,
                         space=spc,
                         main=main,
                         yaxt="n",
                         cex.names=cxn,
                         font=cxf)

  if(is.null(yaxt)) {
    axis(2)
  }
        
                        #las=1, axes=axes, col=cols, space=c(0, 0.4),
                        #cex.names=0.8, ylim=ylim, xpd=xpd)
    s           <- segments(barcenters, matmean-mateb,
                            barcenters, matmean+mateb,
                            lwd=1.0)
    a           <- arrows  (barcenters, matmean-mateb,
                            barcenters, matmean+mateb,
                            lwd=1.0, angle=90, code=3, length=0.02)
  #axis(1, at=1.5, labels="")
  if (!is.null(abl)) abline(h=abl, lwd=2)
  mtext(ylb, 2, 3, font=2)
  #lines(0:6, rep(-0.01, 7), lwd=2)
  #box(bty="l")
  #box(bty="7", col="white")
}

plotgroupacqrev <- function(df=NULL,
                            xlim=c(-0.3, 0.3),
                            ylim=c(-0.3, 0.3),
                            col=c("white", "gray82", "gray42"),
                            pch=c(21, 21, 21), cex=2, tcol="gray94",
                            tlbs=c("A-R-", "A-R+", "A+R-", "A+R+"),
                            tc=matrix(c( 0.15, -0.15,
                                         0.15,  0.2,
                                        -0.15,  0.2,
                                        -0.15, -0.15),
                                      nrow=4, ncol=2,
                                      byrow=TRUE),
                            tcex=3,
                            xlab="Differential SCR in Acquisition",
                            ylab="Differential SCR in Revesal") {
  # Plot reversal as a function of acquisition
  #
  #   Args:
  #   Returns:
  
  par(mar=c(2, 2, 2, 2), oma=c(14, 4, 4, 4))
  layout(matrix(1))

  #layout(matrix(c(1:100), byrow=TRUE, ncol=10))
  dfa <- subset(df, trial==2)

  # set xlim and ylim
  xlim <- c(-max(dfa$acq, na.rm=TRUE), max(dfa$acq, na.rm=TRUE))
  ylim <- c(-max(dfa$rev, na.rm=TRUE), max(dfa$rev, na.rm=TRUE))

  # set coords for background text
  tc <- matrix(c(xlim[1]/2, ylim[1]/2,
                 xlim[1]/2, ylim[2]/2,
                 xlim[2]/2, ylim[1]/2,
                 xlim[2]/2, ylim[2]/2),
               nrow=4,
               ncol=2,
               byrow=TRUE)

  # create plot
  plot(dfa$acq[dfa$group=="HCC"],
       dfa$rev[dfa$group=="HCC"],
       xlim=xlim, ylim=ylim,
       pch=1, xlab="", ylab="", cex=2, col="white")

  for (i in 1:nrow(tc)) {
    # plot conditions
    text(tc[i, 1], tc[i, 2],
         labels=tlbs[i], cex=3, col=tcol, font=2)
  }

  # plot axis labels
  mtext(xlab, 1, 3, font=2, cex=1.5) 
  mtext(ylab, 2, 3, font=2, cex=1.5) 

  # plot points for each group
  for (i in 1:nlevels(dfa$group)) {
    points(dfa$acq[dfa$group==levels(dfa$group)[i]],
           dfa$rev[dfa$group==levels(dfa$group)[i]], pch=pch[i],
           col="black", bg=col[i],
           cex=2)
  }

  # plot category lines
  abline(h=0, lty=3)
  abline(v=0, lty=3)

  # plot legend
  legend("topleft", legend=levels(dfa$group),
         col=c("black", "black", "black"),
         pt.bg=col,
         pch=pch,
         bty="n", pt.cex=2)
}

regroupcaps <- function(fn=NULL, df=NULL, oc="scr",
                        resptype="cr", transf=NULL,
                        subids=NULL, cutoff=50) {
  # Regroup the participants according to cutoff in CAPS
  #
  # Thanks Philipp! Collapsing to the VSS and PTSD make sense clinically
  # to me as the PTSD cutoff is arbitrary in many respects. When we use
  # the VCC/PTSD combined group we should not do a group comparison to
  # the HCC. We can use more rigid PTSD diagnosis cutoff and include in
  # the PTSD group only those veterans with CAPS =>50 and mover the <50
  # to VCC and then try to re-run the between groups analyses (HCC, VCC,
  # PTSD).
  # Args:
  # Returns:
  r            <- readdat(fn=fn, df=df, resptype=resptype, subids=subids)
  #r$outcome    <- transfy(r[, oc], transf) 
  rr           <- subset(r, trial==2, select=c(id, group, caps))
  rr$newgroup[rr$caps>=cutoff]    <- "VPTSD"
  rr$newgroup[rr$caps<cutoff]     <- "VCC"
  rr$newgroup[rr$group=="HCC"]    <- "HCC"
  df <- merge(r, rr, all=TRUE)
  #df <- rr
  df$newgroup <- factor(df$newgroup,
                        levels=c("HCC", "VCC", "VPTSD"))
  df  <- df[order(df$id, df$trial), ]
  return(df)
}


doplotbar <- function(df=NULL, resptype="cr",
                      oc="scr", transf=NULL, subids=NULL,
                      col=c("white", "gray"), ylim=NULL,
                      main=NULL, ylab=NULL) {
  # Summarize by group, stage, stim with bar plots
  #
  # Args:
  #
  df         <- readdat(df=df, resptype=resptype, subids=subids)
  df$outcome <- transfy(df[, oc], transf) 

  par(mar=c(4, 1, 2, 0), oma=c(5, 7, 5, 1))
  l <- layout(matrix(c(1, 1, 2, 2, c(3:6)), nrow=2, ncol=4,
                     byrow=TRUE))
  #layout.show(l)

  
  dfm <- ddply(subset(df, !stim %in% c("CSplusUS", "CSminusUS")),
               c("group", "stage", "stim", "id"), summarize,
               m=mean(outcome, na.rm=TRUE), sd=sd(outcome, na.rm=TRUE),
               n=sum(!is.na(outcome)), eb=sd/sqrt(n))


  dfms <- ddply(dfm, c("group", "stage", "stim"), summarize,
                mm=mean(m, na.rm=TRUE), sd=sd(m, na.rm=TRUE),
                n=sum(!is.na(m)), eb=sd/sqrt(n))
  dfms$stage <- factor(dfms$stage)
  dfms$stim <- factor(dfms$stim)
  if(is.null(ylim)) {
    ylim <- c(0, max(dfms$mm) + 3 * max(dfms$eb))
  }
  

  print(dfms)

  for (i in 1:nlevels(dfms$group)) {
    print(i)
    if (i == 1) {
      yaxt <- NULL
    } else {
      yaxt = "n"
    }
    
    dfmm <- subset(dfms, group==levels(dfms$group)[i])
    m <- matrix(ncol=2, nrow=2, dfmm$mm[1:4])
    e <- matrix(ncol=2, nrow=2, dfmm$eb[1:4])
    colnames(m) <- levels(dfmm$stage)
    rownames(m) <- levels(dfmm$stim)
    colnames(e) <- levels(dfmm$stage)
    rownames(e) <- levels(dfmm$stim)
    plotbar(m, e, spc=c(0, 0.7), col=col, ylim=ylim, yaxt=yaxt,
            main=levels(dfms$group)[i])
  }


  dfm <- ddply(subset(df, !stim %in% c("CSplusUS", "CSminusUS")),
               c("group", "stage", "phase", "stim", "id"), summarize,
               m=mean(outcome, na.rm=TRUE), sd=sd(outcome, na.rm=TRUE),
               n=sum(!is.na(outcome)), eb=sd/sqrt(n))


  dfms <- ddply(dfm, c("group", "stage", "phase", "stim"), summarize,
                mm=mean(m, na.rm=TRUE), sd=sd(m, na.rm=TRUE),
                n=sum(!is.na(m)), eb=sd/sqrt(n))
  dfms$stage <- factor(dfms$stage)
  dfms$stim <- factor(dfms$stim)

  # now by phases
  for (i in 1:nlevels(dfms$group)) {
    for (j in 1:nlevels(dfms$stage)) {
      print(i)
      if (i == 1 && j == 1) {
        yaxt <- NULL
      } else {
        yaxt = "n"
      }
      
      dfmm <- subset(dfms, group==levels(dfms$group)[i]&
                           stage==levels(dfms$stage)[j])
      m <- matrix(ncol=2, nrow=2, dfmm$mm[1:4])
      e <- matrix(ncol=2, nrow=2, dfmm$eb[1:4])
      colnames(m) <- levels(dfmm$phase)
      colnames(m) <- c("E", "L") 
      rownames(m) <- levels(dfmm$stim)
      colnames(e) <- levels(dfmm$phase)
      colnames(e) <- c("E", "L") 
      rownames(e) <- levels(dfmm$stim)
      #plotbar(m, e, col=col, ylim=ylim, yaxt=yaxt, spc=c(0, 0),
      #        main=levels(dfms$stage)[j])

      plotbar(m, e, col=col, ylim=ylim, yaxt=yaxt, spc=c(0, 0))
    }
  }
  mtext(ylab, 2, 3, font=2, cex=1.5, outer=TRUE)

}

plotcorrplot <- function(df=NULL) {
  # plot correlation matrix of caps factors
  library('corrplot')
  dfav <- subset(df, trial==2&!group=="HCC"&!is.na(caps))
  m    <- dfav[, c("factor1caps", "factor2caps", "factor3caps",
                   "factor4caps", "factor5caps")]
  cm   <- cor(m)
  pdf("../output/figures/cfs_corrcaps.pdf")
  corrplot(cm, method="color", tl.col="black")
  dev.off()
}

calcsubsetmai <- function(df=NULL, subset="odd") {
  # Calculates the mean awareness index for a subset of trials
  #
  # Args:
  #   df: data frame
  #   subset: string (one of "odd" or "even")
  # Returns:
  #   dfp: data frame
  switch(subset,
         # restrict to CR
         odd = {
           dfp <- subset(df, trial %% 2 == 1)
         },
         even = {
           dfp <- subset(df, trial %% 2 == 0)
         })
  
         d             <- aggregate(dfp$bleedlevel, by=list(dfp$id),
                                    FUN=mean,
                                    na.rm=TRUE)
         d$id          <- d$Group.1
         d$id          <- factor(d$id)
         d$meanpercsub <- d$x
         d             <- d[c(3, 4)]
         d$id          <- factor(d$id)
         #dfp           <- merge(dfp, d, all=TRUE)
         return (d)
}

createprobflag <- function(df=NULL, ollim=6) {
  # Flags subjects that might have been tracking the spider
  #
  # Args:
  #  df: data frame
  #  ollim: overlap limit (6)
  # Returns:
  #  dfp: processed data frame
  # set overlap limit

  df         <- subset(df, group=="Unaware"&!is.na(id))
  tmp        <- subset(df, select=c(id, visrespnum, spidernum))
  tmp$overlap = !xor(tmp$spidernum, tmp$visrespnum)
  tmp[, 2:3] <- NULL
  
  hasproblem <- vector(mode="character", length=nlevels(tmp$id))
  id         <- vector(mode="character", length=nlevels(tmp$id))
  for (i in 1:nlevels(tmp$id)) {
    #print(levels(tmp$id)[i])
    s <- subset(tmp, id==levels(tmp$id)[i]&!is.na(overlap))

    # skip if s is empty
    if (nrow(s) == 0) {
      hasproblem[i] <- FALSE
    } else {
      hasproblem[i] <- (max(diff(cumsum(!s$overlap), lag=ollim),
                            na.rm=TRUE) >= ollim) |
        (max(diff(cumsum( s$overlap), lag=ollim), na.rm=TRUE) >= ollim)
    }
    id[i] <- levels(tmp$id)[i]

  }
  return(data.frame(id=id, hasproblem=hasproblem))
}

plotextrapol <- function(df=NULL, main=NULL, line=TRUE) {
  # Plots the regression of awareness on reversal learning
  #
  # The idea is to plot the strength of reversal learning against the
  # mean awareness index. If we hadn't achieved full suppression, we
  # could still extrapolate the regression line so that it intersects
  # the y axis. If the intercept is greater than zero, we can conclude
  # that even under full suppression the indirect measure is non-zero.
  # However, the usual caveats against prediction outside the known
  # range apply (see for instance:
  # http://stats.stackexchange.com/questions/219579/)
  #
  # References: Hannula et al. 2005, Nat Rev Neurosci
  #
  # Args:
  # Returns:
  r <- df
  plot(x=r$meanperc, y=r$revlearn, ylab="", xlab="", pch=19)
  mtext("Awareness index", 1, 2.2, font=2)
  mtext("Reversal learning index", 2, 2.2, font=2)
  mtext(main, 3, 0.2, font=2)
  l <- lm(revlearn ~ meanperc, data=r)

  if (line==TRUE) abline(l, lwd=2)

  print(summary(l))
  print(confint(l, level=0.95))
}


plot_intercepts <- function(df=NULL, oc="logdcm2") {
  # Create regression plots for subsamples (even and uneven trials)
  #
  # Args:
  # Returns:

  dfa <- subset(df, trial==2)
  dfa$meanpercsub <- dfa$meanperc
  df$outcome    <- df[, oc]
  dfs <- subset(df,
                select=c(id, group, stage, stim, trial, ctrial,
                         outcome, stait, stais, fsq, meanperc,
                         randomguess)) 
  dfp <- subset(df, trial==2, select=c(id, hasproblem))
  #dfs <- readdat(resptype="crur")
  
  par(mar=c(3, 3, 3, 3), oma=c(1, 1, 1, 1))
  layout(matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE))


  # even trials
  #dffe <- createfulldf(df=subset(df, trial %% 2 == 1), oc=ocscr[h],
  #                     resptype="cr", ocsur=ocsur[h])

  out             <- calcstrengthint(df=subset(df, trial %% 2 == 1),
                                     fcts=c("id", "stage", "stim"),
                                     invf=-1)
  dffe           <- out[[2]]
  dffe[, "revindex"]  <- dffe$m
  dffe$m         <- NULL
  dffe$stage     <- NULL
  dfee <- merge(dfs, dffe)
  dfee <- merge(dfee, dfp)
  deven <- calcsubsetmai(df=df, subset="even")
  df1 <- merge(dfee, deven, all=TRUE)
  dfe <- subset(df1, trial==2&group=="Unaware")


  # odd trials
  #dffo <- createfulldf(df=subset(df, trial %% 2 == 0), oc=ocscr[h],
  #                     resptype="cr", ocsur=ocsur[h])

  out             <- calcstrengthint(df=subset(df, trial %% 2 == 0),
                                     fcts=c("id", "stage", "stim"),
                                     invf=-1)
  dffo           <- out[[2]]
  dffo[, "revindex"]  <- dffo$m
  dffo$m         <- NULL
  dffo$stage     <- NULL
  dfoo <- merge(dfs, dffo)
  dfoo <- merge(dfoo, dfp)
  dodd  <- calcsubsetmai(df=df, subset="odd")
  df2 <- merge(dfoo, dodd, all=TRUE)
  dfo <- subset(df2, trial==2&group=="Unaware")

  coefs <- vector(mode="numeric", 12)
  ci <- vector(mode="numeric", 12)
  tstat <- vector(mode="numeric", 12)
  pval <- vector(mode="numeric", 12)

  # linear regressions
  ds <- list(dfa, dfe, dfo)
  fm <- c(
    "revindex ~ meanpercsub",
    "revindex ~ stait + meanpercsub",
    "revindex ~ stait + hasproblem + meanpercsub",
    "revindex ~ stait + hasproblem + randomguess + meanpercsub"
    )
  c <-  0
  ti <- c(5, 7, 9, 11)
  pi <- c(7, 10, 13, 16)
  for (i in 1:length(ds)) {
    for (j in 1:length(fm)) {
      c <- c + 1
      l <- lm(fm[j], data=ds[[i]])
      print(summary(l))
      coefs[c] <- l$coef[1]
      ci[c] <- - confint(l)[1] + coefs[c] 
      tstat[c] <- summary(l)$coef[ti[j]]
      pval[c] <- summary(l)$coef[pi[j]]
    }
  }


  dfe$evenodd   <- "even"
  dfo$evenodd   <- "odd"
  # merge the data sets for paired tests
  dfeo <- rbind(subset(dfe, select=c(id, meanpercsub, revindex, evenodd)),
                subset(dfo, select=c(id, meanpercsub, revindex, evenodd)))

  dfe$meanperc <- dfe$meanpercsub
  d <- prepntile(subset(dfe, meanperc>0), 3)
  d$ntile <- d$ntile + 1
  dfee <- merge(dfe, d, all=TRUE)
  dfee$ntile[dfee$meanperc==0] <- 1
   symsize <- ddply(dfee, "ntile", summarize, symsize=sum(!is.na(ntile)))
  print(subset(dfee, select=c(id, revindex, meanperc, ntile)))
  dfee$x <- dfee$meanperc
  dfee$y <- dfee$revindex
  dfoe <- merge(dfo, d, all=TRUE)
  dfoe$ntile[dfee$meanperc==0] <- 1
  dfoe$meanperc <- dfoe$meanpercsub
  dfoe$x <- dfoe$meanperc
  dfoe$y <- dfoe$revindex

  # calculate correlation between mean perc for even/odd
  print(cor.test(dfe$meanpercsub, dfo$meanpercsub))

  # find N per ntile
  nt <- ddply(dfee, "ntile", summarize, nntile=sum(!is.na(ntile))) 

  # plot the scatter using the even trial colorcode
  ntileplot(data.frame(id=dfee$id, x=dfee$meanperc, y=dfoe$meanperc,
                          ntile=dfee$ntile),
            xlab="Awareness index (even)",
            ylab="Awareness index (odd)",
            ntilec=c(paste("UA (N=", nt$nntile[1], ")", sep=""),
                     paste("T1 (N=", nt$nntile[2], ")", sep=""),
                     paste("T2 (N=", nt$nntile[3], ")", sep=""),
                     paste("T3 (N=", nt$nntile[4], ")", sep="")),
            sy=FALSE, reg=TRUE,
            xlim=c(0, 1),
            ylim=c(0, 1))
  mtext("a", 3, 2, adj=0, at=-0.3, font=2, cex=1.2)

  ntileplot(dfee, main="Even trials",
            ntilec=c(paste("UA (N=", nt$nntile[1], ")", sep=""),
                     paste("T1 (N=", nt$nntile[2], ")", sep=""),
                     paste("T2 (N=", nt$nntile[3], ")", sep=""),
                     paste("T3 (N=", nt$nntile[4], ")", sep="")),
            xlim=c(0.0, 1.0), ylim=c(-2, 8.2), sy=FALSE)
  mtext("b", 3, 2, adj=0, at=-0.3, font=2, cex=1.2)

  ntileplot(dfoe, main="Odd trials", 
            ntilec=c(paste("UA (N=", nt$nntile[1], ")", sep=""),
                     paste("T1 (N=", nt$nntile[2], ")", sep=""),
                     paste("T2 (N=", nt$nntile[3], ")", sep=""),
                     paste("T3 (N=", nt$nntile[4], ")", sep="")),
            xlim=c(0.0, 1.0), ylim=c(-2, 8.2), sy=FALSE)
  mtext("c", 3, 2, adj=0, at=-0.3, font=2, cex=1.2)
             
  # plot the regression
  dfe$meanperc <- dfe$meanpercsub
  dfo$meanperc <- dfo$meanpercsub
  plotextrapol(dfa, main="All trials", line=TRUE)
  mtext("d", 3, 2, adj=0, at=-0.3, font=2, cex=1.2)

  #plotextrapol(dfo, main="Odd trials")
  
  # use a dot plot for visualization
  par(mar=c(8, 5, 6, 5), oma=c(6, 4, 6, 4))
  layout(matrix(c(1:1), nrow=1, ncol=1, byrow=TRUE))
  
  labs <- c("All trials", "adj. STAI-T", "adj. STAI-T + overlaps",
            "adj. STAIT + overlaps + predictability", "Even trials",
            "adj. STAI-T", "adj. STAI-T + overlaps",
            "adj. STAIT + overlaps + predictability", "Odd trials",
            "adj. STAI-T", "adj. STAI-T + overlaps",
            "adj. STAIT + overlaps + predictability")
  dftmp <- data.frame(ctrial=c(1:12), m=coefs, eb=ci, labels=labs,
                      pval=pval)
  print(dftmp)
  mxy <- max(dftmp$m) * 1.95
  plot(dftmp$ctrial, dftmp$m, col="white", xlab="", ylab="",
       ylim=c(-mxy*0.5, mxy), xaxt="n", bty="n")
  plotscr(df=dftmp, l=FALSE, pch=19, lwd=1.2, pcex=0.8)
  print("here")
  abline(h=0, lty=3)
  axis(1, at=c(1:length(dftmp$labels)), labels=dftmp$labels, las=2, font=2)
  mtext("Intercept with 95% CI", 2, 3, font=2)


  # plot catgorical
  dfee$trials <- "Even"
  dfoe$trials <- "Odd"
  dd <- rbind(dfee, dfoe)
  dd$trialsn[dd$trials=="Even"] <- 0
  dd$trialsn[dd$trials=="Odd"] <- 1
  ntileplot(data.frame(id=dfee$id, x=dd$trialsn, y=dd$meanperc,
                          ntile=dd$ntile),
            xlab="Trials",
            ylab="Awareness index",
            ntilec=c("T1", "T2", "T3"), pl=TRUE)

  ddd <- data.frame(id=dd$id, x=dd$trials, y=dd$meanperc,
                    ntile=dd$ntile, trials=c(dfee$trials, dfoe$trials))
  #plot(rep(c(0, 1), rep(nrow(ddd)/2, 2)), c(ddd$y[ddd$trials=="Even"],
  #                                          ddd$y[ddd$trials=="Odd"]))
  #segments(rep(0, nrow(ddd)/2), ddd$y[ddd$trials=="Even"],
  #         rep(1, nrow(ddd)/2), ddd$y[ddd$trials=="Odd"])

  plot(dd$id[dd$trials=="Even"],
       dd$meanperc[dd$trials=="Even"],
       pch=19,
       cex=0.7)
  points(dd$id[dd$trials=="Odd"],
         dd$meanperc[dd$trials=="Odd"],
         pch=23,
         cex=0.7)

}

sumstats <- function(df=NULL) {
  # produce summary stats
  library(tidyr)
  library(dplyr)

  # filter and select data
  df$gend <- as.numeric(df$gender)
  df$gr <- as.numeric(df$group)
  dfa <- filter(df, trial==2, !manipulation=="")
  #dfa <- select(dfa, manipulation, gend, gr, age, stais, stait, fsq,
  #              meanperc, meanresp, meanconf, correctresp, meanur, acq,
  #              rev, revindex) 

  dfa <- subset(dfa, select=c(manipulation, gend, gr, age, stais, stait, fsq,
                meanperc, meanresp, meanconf, correctresp, meanur, acq,
                rev, revindex)) 

  ## reshape2 still does its thing:
  library(reshape2)
  melted <- melt(dfa, id.vars=c("manipulation"))
  head(melted)
  
  library(dplyr)
  grouped <- group_by(melted, manipulation, variable)
  s <- summarise_each(grouped, funs(n=sum(!is.na(value)),
                                    missing=sum(is.na(value)),
                                    mean=round(mean(value, na.rm=TRUE), 2),
                                    sd=round(sd(value, na.rm=TRUE), 2),
                                    se=round(sd/sqrt(n), 2)))
  sw <- cbind(filter(s, manipulation=="noCFS"),
              filter(s, manipulation=="CFS"))
  colnames(sw) <- c(paste(colnames(s), "NoCFS", sep=""),
                    paste(colnames(s), "CFS", sep=""))
  
  # use dplyr to get wide format
  #su <- unite(s, variable, n, r)
  #sw <- spread(su, n, missing, mean, sd, se)

  # pooled sd
  #sp <- sqrt(((sw$nNoCFS-1)*sw$sdNoCFS^2 + (sw$nCFS-1)*sw$sdCFS^2)/(sw$nNoCFS + sw$nCFS - 2))
  
  sw$tstat <- (sw$meanNoCFS - sw$meanCFS)/sqrt(sw$sdNoCFS^2/sw$nNoCFS + sw$sdCFS^2/sw$nCFS)
  #sw$tstat <- (sw$meanNoCFS - sw$meanCFS)/sqrt(sw$seNoCFS + sw$seCFS)
  #sw$tstat <- (sw$meanNoCFS - sw$meanCFS)/(sp * sqrt(1/sw$nNoCFS + 1/sw$nCFS))
  sw$Pval  <- round(2*pt(abs(sw$tstat), (sw$nNoCFS + sw$nCFS -2), lower=FALSE), 4)
  sw$tstat <- round(sw$tstat, 2)
  sw$df    <- (sw$nNoCFS + sw$nCFS -2)
  colnames(sw) <- c(paste(colnames(s)), paste(colnames(s)), "t", "P", "df")
  return(sw)

}

ntileplot <- function(df=NULL,
                      col=c("dodgerblue3", "firebrick3",
                            "goldenrod3", "darkseagreen3"),
                      main=NULL, xlab="Awareness index",
                      ylab="Reversal learning index",
                      ntilec=c("Q1", "Q2", "Q3", "Q4"),
                      pl=FALSE,
                      sy=FALSE,
                      xlim=c(0, 1),
                      ylim=c(0, 1),
                      reg=FALSE) {
  # Plot regression using ntile color coding as in Shanks 2016, Fig
  # 3a

  plot(x=df$x, y=df$y, col="white",
       xlab="", ylab="", xlim=xlim, ylim=ylim)
  mtext(xlab, 1, 2.2, font=2)
  mtext(ylab, 2, 2.2, font=2)
  mtext(main, 3, 0.2, font=2)
  for (i in 1:max(df$ntile)) {
    p <- subset(df, df$ntile==i)
    points(p$x, p$y, col=col[i],
           bg=col[i], pch=19, cex=1.0)
    #print(col[i])
    if (pl == TRUE) {
      segments(p$x[p$x==0], p$y[p$x==0],
               p$x[p$x==1], p$y[p$x==1], lwd=0.4, lty=3)
      print(p)
    }
    if (sy == TRUE) {
      x <- mean(p$x)
      y <- mean(p$y)
      points(x, y, cex=nrow(p)*0.05, bg=col[i], pch=21,
             col="black")
    }
    
  }
  if (reg == TRUE) {
    l <- lm(y ~ x, data=df)
    abline(l, lty=1)
    rsq <- summary(l)$r.squared
    text(0.1, 0.8, paste("r=", round(sqrt(rsq), 2), sep=""), cex=0.85)
  }
  
  legend("top", col=col, pch=rep(19, max(df$ntile)),
         cex=0.5, bty="n", ncol=4,
         legend=ntilec)
}

prepntile <- function(df=NULL, nt=4) {
  # Prepare ntiles
  df$ntile <- ntile(df$meanperc, nt)
  d <- subset(df, select=c(id, ntile))
  return(d)
}

revsims <- function(nsims=10^5) {
  # Simulate a few scenarios of reversal learning
  #
  # Args:
  #   nsims: # of simulations

  par(mar=c(7, 5, 5, 5), oma=c(1, 2, 1, 2))
  l        <- layout(matrix(c(1:6), byrow=TRUE, ncol=3))

  scenarios <- c("Acq sig, Rev sig\nRev Ind sig\nStage diff n.s.",
                 "Acq sig, Rev n.s.\nRev Ind sig\nStage diff n.s.",
                 "Acq sig, Rev sig.\nRev Ind sig\nStage diff n.s.",
                 "Acq sig, Rev n.s.\nRev Ind sig\nStage diff n.s.",
                 "Acq sig, Rev n.s.\nRev Ind sig\nStage diff n.s.",
                 "Acq sig, Rev n.s.\nRev Ind sig\nStage diff n.s.")

  scenarios <- c("Scenario 1\nNon-zero Reversal Index\nZero Stage Difference",
                 "Scenario 2\nNon-zero Reversal Index\nNon-zero Stage Difference",
                 "Scenario 3\nZero Reversal Index\nNon-Zero Stage Difference",
                 "Scenario 4\nZero Reversal Index\nZero Stage Difference")

  
  da <- c( 1.2,  1.8,  0.5,  0.25)
  dr <- c(-0.7, -0.5,  0.2, -0.25)
  sigma <- c(7.5, 2.5, 10, 10)
  r <- c(0.5, 0.3, 0.3, 0.3)
  set.seed(432)

  for (i in 1:length(scenarios)) {
    lbs    <- c("Acquisition",
                "Reversal",
                "Reversal Index",
                "Stage Difference")
    hccn   <- 40
    X0     <- mvrnorm(n=hccn, mu=c(da[i], dr[i]),
                      Sigma=rbind(c(sigma[i], r[i]),
                                  c(r[i], sigma[i])))
    hccacq <- X0[,1] 
    hccrev <- X0[,2]
    hccri  <- X0[,1] - X0[,2]
    hccdc  <- X0[,1] + X0[,2]
    df     <- data.frame(m      =c(mean(hccacq),
                                   mean(hccrev),
                                   mean(hccri),
                                   mean(hccdc)),
                         eb     =c(2*(sd(hccacq)/sqrt(hccn)),
                                   2*(sd(hccrev)/sqrt(hccn)),
                                   2*(sd(hccri)/sqrt(hccn)),
                                   2*(sd(hccdc)/sqrt(hccn))),
                         ctrial =c(1:4),
                         group  =c("Acq", "Rev"))
    plot(c(1:4), df$m, col="white",
         ylim=c(-3.5, 3.5),
         xlim=c(1, 4),
         axes=FALSE, xlab="", ylab="",
         main=scenarios[i], cex.main=0.75)
    #axis(1, at=c(1:4), labels=c("", "", "", ""))
    axis(1, at=c(1:4), labels=lbs, las=2, font=2)
    if (i == 1 || i == 3) {
      #axis(2)
    }
    axis(2)
    plotscr(df, l=FALSE, pch=19, lwd=1.2)
    abline(h=0, lty=3)

    # run t test as control
    print(t.test(hccacq, hccrev, paired=TRUE))
    print(t.test(hccri))
    print(t.test(hccdc))

    if (i == 2 || i == 4) {
      plot.new()
    }
  }
  mtext("Cohen's d with 95% CI", 2, -1, outer=TRUE, font=2, cex=1.5)
  text(-0.9, 3.075, "Indices:", xpd=NA, adj=0, cex=0.7, font=2)
  text(-0.9, 3.0, "Acquisition: Face A - Face B", xpd=NA, adj=0, cex=0.7)
  text(-0.9, 2.925, "Reversal: Face A - Face B", xpd=NA, adj=0, cex=0.7)
  text(-0.9, 2.85, "Reversal Index: (Face A - Face B) - (Face A - Face B)",
       xpd=NA, adj=0, cex=0.7)
  text(-0.9, 2.775, "Stage Difference: (Face A - Face B) - (Face B - Face A)",
       xpd=NA, adj=0, cex=0.7)
}


rw <- function(v0=0.5, alpha=0.05, ntrials=32,
               stim=NULL, lambda=NULL, shock=NULL) {
  v <- vector("numeric", ntrials)
  for (i in 1:ntrials) {
    v[i] <- v0
    v0 <- v0 + alpha * (lambda[i] - v0) 
  }
  return(v)
}


  
simtrialseq <- function(stim=NULL, ntrials=32,
                        lambda=NULL, alpha=0.05,
                        seed=20170719, rr=0.33) {
  #
  #
  # Run the trial sequence simulations

  # The problem is this: assuming a trial-switch expectancy that updates
  # from trial to trial, how do we map this on to the expected value?
  # For reasonable trial sequences, an alternating 0-1 pattern seems to
  # be helpful; but what about an imaginery trial sequence of just CS-?
  # Here, the trial-switch expectancy would be zero after just a few
  # trials (depending on the learning rate), and so should be the
  # expected value.
  # Thus, the following model might not be enough:
  #
  # V = Vts * s + (1 - s) * (1 -Vts)
  #
  # With zero switch expectancy after a few CS- trials, this would
  # result in V = 1.
  #
  # V = Vts * s + (1 - s) * (1 - Vts)
  # if Vts > 0.5 & s == 1: snew = 0
  # if Vts > 0.5 & s == 0: snew = 1
  # else snew = s
  #
  # V = Vts * s + (1 - s) * (1 - Vts) iff Vts >= 0.5
  # V = Vts * s iff Vts < 0.5
  #
  # The model assumes a trial switch when the expectancy is > 0.5:
  #
  # s = 1 iff Vts > 0.5
  # s = 0 other wise
  #
  # V = Vts * lambda  + (1 -lambda) * (1 - Vts) iff s = 1
  # V = Vts otherwise

  set.seed(seed)
  library(data.table)

  if (is.null(stim)) {
    # create stim sequence
    # do not allow more than 2 equal stimuli in a row
    stim = sample(c(rep(0, 16), rep(1, 16)))
    done <- FALSE
    while(!done) {
      stim = sample(c(rep(0, 16), rep(1, 16)))
      stimi = abs(stim-1)
      stim[17] <- 1
      done <- max(diff(cumsum(stim), lag = 3)) < 3 && 
        max(diff(cumsum(stimi), lag = 3)) < 3 &&
        sum(stim[1:ntrials/2]) == ntrials/2
    }
    print(stim)
  }
  
    
  # create lambda
  if (is.null(lambda)) {
    done <- FALSE
    while(!done) {
      lambda <- sample(c(rep(1, round(ntrials * rr)),
                         rep(0, ntrials-round(ntrials*rr))))
      lambda[17] <- 1
      done <- (max(diff(cumsum(lambda), lag = 3)) < 3) &&
        sum(stim[lambda==1] == 1) == sum(lambda==1)
    }


    # switch sign in second half of stim
    stim[(ntrials/2+1):ntrials] <- abs(stim[(ntrials/2+1):ntrials] - 1)
  }
  
    

  # RW-model
  v <- vector("numeric", ntrials)
  vm <- vector("numeric", ntrials/2)
  vp <- vector("numeric", ntrials/2)
  vm <- rw(alpha=alpha, lambda=lambda[stim==0],
           ntrials=ntrials/2, v0=0.5)
  vp <- rw(alpha=alpha, lambda=lambda[stim==1],
           ntrials=ntrials/2, v0=0.5)
  v[stim==0] <- vm
  v[stim==1] <- vp
  trials <- c(1:ntrials)

  # layout
  layout(matrix(c(1:4), nrow=2, byrow=TRUE))
  par(mar=c(3, 3, 3, 2), oma=c(3, 3, 3, 0))
  
  plot(c(1:ntrials), v, lty=3, "o", pch=19, ylim=c(0, 1.2),
       cex=0.5, xlab="", ylab="", yaxt="n", yaxt="n")
  points(trials[stim==1], v[stim==1], lty=3, pch=19, col="red",
         cex=0.5)
  points(trials[stim==0], v[stim==0], lty=3, pch=19, col="blue",
         cex=0.5)
  abline(v=16.5, lty=2, col="gray")
  legend("top", legend=c("CS-", "CS+"), pch=c(19, 19),
         col=c("blue", "red"), bty="n", ncol=2, cex=0.5)
  axis(2, at=seq(0, 1, 0.2))
  mtext("Conventional RW-model", 3, 0, font=2, cex=0.75)
  mtext("Expected value", 2, 2.5, font=2, cex=0.75)
  plot.new()
  mtext(paste("Reinforcement rate: ", rr*100, "%\nLearning rate: ",
              alpha, sep=""), 3, -3, font=2)


  # now the model
  # lambdabin tracking trial shifts
  stimbin <- rep(0, ntrials)
  stimbin[lambda==1] <- 1
  lambdabin <- abs((stimbin == shift(stimbin, 1))-1)
  lambdabin[1] <- 0
  vst <- rw(alpha=alpha, lambda=lambdabin, ntrials=ntrials)
  plot(c(1:ntrials), vst, lty=3, "o", pch=19, ylim=c(0, 1.2),
       cex=0.5, xlab="", ylab="", yaxt="n")
  points(trials[stimbin==1], vst[stimbin==1], lty=3, pch=19,
         col="indianred", cex=0.5)
  points(trials[stimbin==0], vst[stimbin==0], lty=3, pch=19,
         col="lightblue", cex=0.5)
  legend("top", legend=c("Non-reinforced", "Reinforced"), pch=c(19, 19),
         col=c("lightblue", "indianred"), bty="n", ncol=2, cex=0.5)
  axis(2, at=seq(0, 1, 0.2))
  #abline(h=0.5, lty=3)
  mtext("Trial-switch model", 3, 0, font=2, cex=0.75)
  mtext("Trial", 1, 0, outer=TRUE, font=2)
  mtext("Expected trial switch", 2, 2.5, font=2, cex=0.75)

  
  # mapping on to expected value
  vste <- vector("numeric", ntrials)
  for (t in 1:ntrials) {
    if (vst[t] <= 0.5) {
      # no trial switch expected
      vste[t] <- vst[t]
    } else {
      # trials switch expected
      s <- abs(lambda[t-1] - 1)
      vste[t] <- vst[t] * s + (1 - s) * (1 - vst[t])
    }
  }
  plot(c(1:ntrials), vste, lty=3, "o", pch=19, ylim=c(0, 1.2),
       cex=0.5, xlab="", ylab="", yaxt="n")
  points(trials[stim==0], vste[stim==0], lty=3, pch=19, col="blue",
         cex=0.5)
  points(trials[stim==1], vste[stim==1], lty=3, pch=19, col="red",
         cex=0.5)
  legend("top", legend=c("CS-", "CS+"), pch=c(19, 19),
         col=c("blue", "red"), bty="n", ncol=2, cex=0.5)
  axis(2, at=seq(0, 1, 0.2))
  abline(v=16.5, lty=2, col="gray")
  mtext("Predictability model", 3, 0, font=2, cex=0.75)
  mtext("Expected value", 2, 2.5, font=2, cex=0.75)
}

pdf()
alpha <- c(0.1, 0.3, 0.5)
for (a in 1:length(alpha)) {
  stim <-   c(0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0,
              0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0)

  # 100 % reinforcement
  lambda <- c(0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0,
              1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1)
  simtrialseq(seed=3,  rr=1, alpha=alpha[a], stim=stim, lambda=lambda)

  # 75 % reinforcement
  lambda <- c(0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0,
              1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0)
  simtrialseq(seed=3, rr=0.75, alpha=alpha[a], stim=stim, lambda=lambda)
  
  # 50 % reinforcement
  lambda <- c(0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0,
              1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0)
  simtrialseq(seed=3, rr=0.5, alpha=alpha[a], stim=stim, lambda=lambda)

  # 38 % reinforcement
  lambda <- c(0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0,
              1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0)
  simtrialseq(seed=3, rr=0.38, alpha=alpha[a], stim=stim, lambda=lambda)
}
dev.off()
  
  
 
summarize_models <- function () {
  # summarize the model fits
  #
  #
  library(R.matlab)
  m1 <- readMat("../data/cfs_phtdindiv_rw_logdcm2_allbounded.mat")
  m2 <- readMat("../data/cfs_phtdindiv_trialswitch_logdcm2_allbounded.mat")
  df1 <- data.frame(id=rep(unlist(m1$idlist), each=32),
                    trial=rep(c(1:32), 98),
                    model="RW",
                    stim=unlist(m1$CSUSs),
                    v=unlist(m1$VFits),
                    alpha=rep(m1$minsum[, 1], each=32),
                    bic=rep(unlist(m1$BICs), each=32),
                    ll=rep(unlist(m1$llss), each=32),
                    order=unlist(m1$Orders))

  df2 <- data.frame(id=rep(unlist(m2$idlist), each=32),
                    trial=rep(c(1:32), 98),
                    model="Sequence",
                    stim=unlist(m2$CSUSs),
                    v=unlist(m2$VFits),
                    alpha=rep(m2$minsum[, 1], each=32),
                    bic=rep(unlist(m2$BICs), each=32),
                    ll=rep(unlist(m2$llss), each=32),
                    order=unlist(m2$Orders))


  dff <- read.csv("../data/cfs_rwinput.csv")
  dff <- subset(dff, select=c(id, trial, group,
                              meanperc, logdcm2))
  dff1 <- merge(dff, df1)
  dff2 <- merge(dff, df2)

  dfs <- rbind(dff1, dff2) 
  dfm <- ddply(dfs, c("group", "model", "order", "trial"),
               summarize,
               alpha=mean(alpha),
               bic=mean(bic),
               meanv=mean(v),
               sdv=sd(v),
               nv=sum(!is.na(v)),
               sev=sdv/sqrt(nv),
               civ=2*sev,
               means=mean(logdcm2),
               sds=sd(logdcm2),
               ns=sum(!is.na(logdcm2)),
               ses=sds/sqrt(ns),
               cis=2*ses)

  # orders
  # sort dfs by order, trial
  dfs <- dfs[order(dfs$order, dfs$trial),]
    
  dda <- unique(subset(dfs, order=="A", select=c(order, trial, stim)))
  ddb <- unique(subset(dfs, order=="B"&id==9525,
                       select=c(order, trial, stim)))
  ddc <- unique(subset(dfs, order=="C", select=c(order, trial, stim)))
  ddd <- unique(subset(dfs, order=="D", select=c(order, trial, stim)))

  dd <- rbind(dda, ddb, ddc, ddd)

  dfm2 <- merge(dfm, dd)
  dfm <- dfm2

  dfm <- dfm[order(dfm$group, dfm$model, dfm$order, dfm$trial),]


  print(paste("Unaware BIC1 =",
              round(sum(subset(dfs,
                               trial==1& group==
                               "Unaware"& model=="RW")$bic), 2)))

  print(paste("Unaware BIC2 =",
              round(sum(subset(dfs,
                               trial==1& group==
                               "Unaware"& model==
                               "Sequence")$bic), 2)))

  print(paste("All BIC1 =",
              round(sum(subset(dfs, trial==1& 
                                    model=="RW")$bic), 2)))

  print(paste("All BIC2 =",
              round(sum(subset(dfs, trial==1& 
                                    model=="Sequence")$bic), 2)))

  print(paste("Meanperc==0 BIC1 =",
              round(sum(subset(dfs, trial==1& 
                                    meanperc==0&
                                    model=="RW")$bic), 2)))

  print(paste("Meanperc==0 BIC2 =",
              round(sum(subset(dfs, trial==1&
                                    meanperc==0&
                                    model=="Sequence")$bic), 2)))


  # write results to disk
  bics <- c(
    sum(subset(dfs, trial==1&group=="Unaware"&model=="RW")$bic),
    sum(subset(dfs, trial==1&group=="Unaware"&model=="Sequence")$bic),
    sum(subset(dfs, trial==1&model=="RW")$bic),
    sum(subset(dfs, trial==1&model=="Sequence")$bic),
    sum(subset(dfs, trial==1&meanperc==0&model=="RW")$bic),
    sum(subset(dfs, trial==1&meanperc==0&model=="Sequence")$bic)
    )
            
  bicsdf <- data.frame(bic=bics, group=c("CFS", "CFS", "All", "All",
                                         "CFS", "CFS"),
                       model=rep(c("RW", "Sequence"), 3))
  write.csv(bicsdf, "../output/tables/cfs_bic.csv")
  #print(dfm$bic)

  #print(dfm$bic)

  # plot by group, model, and order
  for (i in 1:nlevels(dfm$group)) {
    for (j in 1:nlevels(dfm$model)) {
      for (k in 1:nlevels(dfm$order)) {
        dfms <- subset(dfm, group==levels(dfm$group)[i]&
                            model==levels(dfm$model)[j]&
                            order==levels(dfm$order)[k])

        if (dfms$group[1] == "Aware") {
          manipulation <- "noCFS"
        } else {
          manipulation <- "CFS"
        }
        #main=paste("Group=", dfms$group[1], ", Model=", dfms$model[1],
        #           ", Order=", dfms$order[1], ", alpha=",
        #           round(dfms$alpha[1], 2),
        #           sep="")

        main=paste("alpha=", round(dfms$alpha[1], 2), sep="")
        #main=paste("Order:", dfms$order[1])
        plot(dfms$trial, dfms$means, "o", pch=19, ylim=c(-0.2, 2.0),
             xlab="", ylab="", yaxt="n", xaxt="n", cex=0.5, lty=3)
        #mtext("Trial", 1, 2.5, font=2)
        #mtext("Log SNA", 2, 2.5, font=2)
        mtext(main, 3, 0, font=2, cex=1.0, adj=1)

        # pearson correlation
        r <- sqrt(summary(lm(dfms$means ~ dfms$meanv))$r.squared)
        p <- summary(lm(dfms$means ~ dfms$meanv))$coef[2, 4]
        as <- starsfromp(p) 
        lb <- paste("r=", round(r, 2), as, sep="")
        mtext(lb, 3, 0, font=2, cex=1.0, adj=0)
        points(dfms$trial, dfms$meanv, "o", pch=19, lty=1, cex=0.5,
               lwd=0.75, col="lightgray")
        points(dfms$trial[dfms$stim=="CSminus"],
               dfms$means[dfms$stim=="CSminus"], 
               pch=19, col="blue", cex=0.5)
        points(dfms$trial[dfms$stim=="CSplus"],
               dfms$means[dfms$stim=="CSplus"], 
               pch=19, col="red", cex=0.5)
        #legend("top", legend=c("CSminus", "CSplus", "Expected value"),
        #       pch=c(19, 19, 19),
        #       col=c("blue", "red", "lightgray"),
        #       bty="n", ncol=3, cex=0.5)

        # plot orders
        if (i == 1 & j == 1) {
          mtext(paste("Order", dfms$order[1]), 3, 2, font=2)
        }

        # plot x axis 
        if (i == 2 & j == 2) {
          axis(1, cex.axis=1.5)
        }

        # plot x axis label
        if (i == 2 & j == 2 & k == 3) {
          mtext("Trial", 1, 2, font=2, outer=TRUE, cex=1.5)
          #legend("bottom", legend=c("CSminus", "CSplus", "Expected value"),
          #       pch=c(19, 19, 19), col=c("blue", "red", "lightgray"),
          #       bty="n", ncol=1, cex=1, xpd=NA, inset=-0.8)
        }

        # plot legend
        if (i == 1 & j == 1 & k == 2) {
          legend(x=2, y=4, legend=c("CSminus", "CSplus", "Expected value"),
                 pch=c(19, 19, 19), col=c("blue", "red", "lightgray"),
                 bty="n", ncol=3, cex=1.5, xpd=NA, inset=-0.8)
        }



        # plot y axis
        if (k == 1) {
          axis(2, at=seq(0, 1, 0.5), las=1, cex.axis=1.5)
        }
        

        # plot y axis label
        if (i == 2 & j == 2 & k == 1) {
          mtext("SNA amplitude", 2, 3, font=2, outer=TRUE,
                cex=1.5)
        }
        

        # plot model
        if (k == 4) {
          #mtext(dfms$model[1], 4, 2, font=2, las=3)
          text(36, 0.9, labels=dfms$model[1], font=2, 
               srt = -90, xpd=NA, cex=1.5)
        }

        # plot group
        if (j == 1 & k == 4) {
          #mtext(dfms$model[1], 4, 2, font=2, las=3)
          text(41, -0.75, labels=manipulation, font=2, 
               srt = -90, xpd=NA, cex=2.0)
        }

        #legend("top", inset=0.05, legend=c("SCR", "fitted"),
        #       lty=c(1, 3), ncol=2, bty="n", cex=0.5)
      }
    }
  }
}

# Recursion for consecutive hits
f <- function(j, m) {
  k <- 8
  p <- 1/2
  if ((k - j) > m) {
    return(0)
  }
  if (j == k & m >= 0) {
    return(1)
  }
  return(p * f(j+1, m-1) + (1 - p) * f(0, m-1))
}

# Extract legend
g_legend<-function(agplot){
  tmp <- ggplot_gtable(ggplot_build(agplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
